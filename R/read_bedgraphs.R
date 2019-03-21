#' Versatile BedGraph reader.
#' @details Reads BedGraph files and generates methylation and coverage matrices. Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files bedgraph files.
#' @param pipeline Default NULL. Can be "Bismark" or "MethylDeckal". If not known use idx arguments for manual column assignments.
#' @param ideal If all bedgraphs contain same number of CpG's and in same order. If TRUE instead of merging tables, matrices are crated by simply cbinding which is significantly faster.
#' @param genome BSgenome object or name of the installed BSgenome package. Example: BSgenome.Hsapiens.UCSC.hg19
#' @param contigs contigs to restrict genomic CpGs to. Default all autosomes and allosomes - ignoring extra contigs.
#' @param zero_based Are bedgraph regions zero based ? Default TRUE
#' @param vect To use vectorized code. Default TRUE, memory intese. Set to FALSE if you have large number of BedGraph files.
#' @param vect_batch_size Default NULL. Process samples in batches. Applicable only when vect = TRUE
#' @param coldata An optional DataFrame describing the samples. Row names, if present, become the column names of the matrix. If NULL, then a DataFrame will be created with basename of files used as the row names.
#' @param chr_idx column index for chromosome in bedgraph files
#' @param start_idx column index for start position in bedgraph files
#' @param end_idx column index for end position in bedgraph files
#' @param beta_idx column index for beta values in bedgraph files
#' @param M_idx column index for read counts supporting Methylation in bedgraph files
#' @param U_idx column index for read counts supporting Un-methylation in bedgraph files
#' @param strand_idx column index for strand information in bedgraph files
#' @param cov_idx column index for total-coverage in bedgraph files
#' @param h5 Should the coverage and methylation matrices be stored as "HDF5Array"
#' @param h5_dir directory to store H5 based object
#' @param verbose Be little chatty ? Default TRUE.
#' @param bored Tell me a Chuck Norris joke while I wait for my CpG extraction. Default TRUE. jokes can be very explicit! you have been warned!
#' @export
#' @import data.table
#' @import parallel
#' @import SummarizedExperiment
#'
#'

#one has to check if it works correctly with Bismark and MethylDackel data.

read_bedgraphs = function(files = NULL, pipeline = NULL, zero_based = TRUE, genome_ob=NULL, genome_name = NULL, contigs = NULL, vect = T, vect_batch_size = NULL, coldata = NULL, chr_idx = NULL,
                          start_idx = NULL, end_idx = NULL, beta_idx = NULL, stranded = T, h5temp=NULL,
                          M_idx = NULL, U_idx = NULL, strand_idx = NULL, cov_idx = NULL, h5 = FALSE, h5_dir = NULL, verbose = TRUE, bored = TRUE){

  if(is.null(files)){
    stop("Missing input files.", call. = FALSE)
  }
  if(is.null(genome)){
    stop("Missing genome. Please provide a valid genome name", .call = F)
  }
  if (vect && is.null(vect_batch_size)){
    vect_batch_size <- length(files)
  }



  #Extract CpG's
  if(is.null(genome_ob)){
    genome = extract_CPGs(ref_genome = genome_name, bored = bored)
  } else {
    message(paste0("Using the provided genome object with ", nrow(genome), "CpG sites."))
      genome <- copy(genome_ob)
      rm(genome_ob)
    }

    if(zero_based){
      genome[, start := start - 1][, end := end - 1]
    }
    #check it with the strand column
    if (stranded){
      genome_plus <- copy(genome)
      genome_plus[, strand := "+"]
      genome[, start := start + 1]
      genome[, strand := "-"]
      genome <- rbindlist(list(genome, genome_plus), use.names = T)
      setkeyv(genome, cols=c("chr", "start"))
      message(paste0("Splitted into  ", nrow(genome), " CpGs with strand information"))
      rm(genome_plus)
    }

    if(is.null(contigs)){
      #Work with only main contrigs (either with chr prefix - UCSC style, or ensemble style)
      ####it should work more generally
      contigs = c(paste0("chr", c(1:22, "X", "Y", "M")), 1:22, "X", "Y", "MT")
    }

    genome_contigs = genome[,.N,chr][,chr]
    genome = genome[chr %in% as.character(contigs)]

    if(nrow(genome) == 0){
      message("No more CpG's left after subsetting for contigs. It appears provided contig names do not match to the BSgenome.")
      message("Contigs provied:")
      print(contigs)
      message("Contigs from BSGenome:")
      print(genome_contigs)
      stop(call. = FALSE)
    }
    if(verbose){

      message(paste0("Retained ", nrow(genome), " CpGs after filtering for contigs"))
    }



  #Set colData
  if(is.null(coldata)){
    coldata = data.frame(row.names = unlist(data.table::tstrsplit(x = basename(files), split = '\\.', keep = 1)),
                         stringsAsFactors = FALSE)
  }else{
    if(length(files) != nrow(coldata)){
      stop("Number of samples in coldata does not match the number of input files.")
    }
  }

  start_proc_time = proc.time()

  #Final aim is to bring input data to the following order:
  ##chr start end beta cov starnd <rest..>
  if(is.null(pipeline)){
    col_idx = .parse_source_idx(chr = chr_idx, start = start_idx, end = end_idx, beta = beta_idx,
                                cov = cov_idx, strand = strand_idx, n_meth = M_idx, n_unmeth = U_idx, verbose = verbose)
    col_idx$col_classes = NULL
  }else{
    pipeline = match.arg(arg = pipeline, choices = c("Bismark_cov", "MethylDackel"))
    if(verbose){
      message(paste0("Using ", pipeline, " as a preset.."))
    }
    col_idx = .get_source_idx()
  }

  #Summarize bedgraphs and create a matrix

  if(vect){
    mat_list = .vect_code_batch(files = files, col_idx = col_idx, batch_size = vect_batch_size, col_data = coldata,  genome = genome)
  } else {
    mat_list = .non_vect_code(files = files, col_idx = col_idx, coldata = coldata, verbose = verbose,  genome = genome, h5 = h5, h5temp = h5temp)
  }

  if(nrow(mat_list$beta_matrix) != nrow(mat_list$cov_matrix)){
    stop("Discrepancies in dimensions of coverage and beta value matrices.")
  }

  if(vect && h5){
    #se = methrix(SE = SummarizedExperiment::SummarizedExperiment(assays = list(methylation_matrix = DelayedArray::DelayedArray(mat_list$beta_matrix[,which(!(colnames(mat_list$beta_matrix) %in% c("chr", "start")))]),
    #                                                                           coverage_matrix = DelayedArray::DelayedArray(mat_list$cov_matrix[,which(!(colnames(mat_list$cov_matrix) %in% c("chr", "start")))]))),
    #             CpGs = nrow(mat_list$cov_matrix), samples = nrow(coldata), h5 = as.logical(h5), genome = genome)
    se <-  methrix(beta=  as(mat_list$beta_matrix, "HDF5Array"),
                    cov= as(mat_list$cov_matrix, "HDF5Array"),
                    position=genome[,.(chr, start, strand)],
                    is.HDF5 = T, genome=genome_name, sample_annotation=coldata)
    rm(mat_list)
    gc()
    if(!is.null(h5_dir)){
      tryCatch(HDF5Array::saveHDF5SummarizedExperiment(x = se, dir = h5_dir, replace = TRUE),  error = function(e) message("The dataset is not saved."))
    }
  }else if (!h5) {
    #se = methrix(SE = SummarizedExperiment::SummarizedExperiment(assays = list(methylation_matrix = mat_list$beta_matrix[,which(!(colnames(mat_list$beta_matrix) %in% c("chr", "start")))],
    #                                                                           coverage_matrix = mat_list$cov_matrix[,which(!(colnames(mat_list$cov_matrix) %in% c("chr", "start")))]), rowData = mat_list$beta_matrix[,1:2]),
    #             CpGs = nrow(mat_list$cov_matrix), samples = nrow(coldata), h5 = as.logical(h5), genome = genome)

    se <-  methrix(beta= mat_list$beta_matrix,
                    cov= mat_list$cov_matrix,
                    position=genome[,.(chr, start, strand)],
                    is.HDF5 = F, genome=genome_name, sample_annotation=coldata)
  } else if (!vect && h5){

    se <-  methrix(beta= mat_list$beta_matrix,
                    cov= mat_list$cov_matrix,
                    position=genome[,.(chr, start, strand)],
                    is.HDF5 = T, genome=genome_name, sample_annotation=coldata)

    if(!is.null(h5_dir)){
      tryCatch(
      HDF5Array::saveHDF5SummarizedExperiment(x = se, dir = h5_dir, replace = TRUE),  error = function(e) message("The dataset is not saved."))
    }
  }

  cat(data.table:::timetaken(started.at = start_proc_time), sep = "\n")

  return(se)

}


#Bismark and methyldackel have same output format
.get_source_idx = function(){
  return(list(col_idx = c(chr = 1, start = 2, end = 3, beta = 4, M = 5, U = 6),
              col_classes = c("character", "numeric", "numeric", "numeric", "integer", "integer"),
              fix_missing = c("cov := M+U", "strand := '.'")))
}

#Parse custom indices and return missing info
.parse_source_idx = function(chr = NULL, start = NULL, end = NULL, strand = NULL,
                             beta = NULL, n_meth = NULL, n_unmeth = NULL, cov = NULL, beta_fract = FALSE, verbose = TRUE){

  #mandatory chr and start field
  if(is.null(chr) | is.null(start)){
    stop("missing chromosome/start indices\nUse pipeline argument if the files are from Bismark or MethyDeckal", call. = FALSE)
  }

  #See if any indices are duplicated
  if(length(which(duplicated(c(chr, start, end, strand, beta, n_meth, n_unmeth, cov)))) > 0){
    stop("Duplicated indices.", call. = FALSE)
  }

  #Check maximum betavalues (Can be 1 or 100)
  if(beta_fract){
    cov_scale = 1
  }else{
    cov_scale = 100
  }
  fix_missing = vector()
  if(is.null(strand)){
    fix_missing = "strand := '*'"
  }

  #If beta and cov are missing
  if(all(is.null(beta), is.null(cov))){
    if(is.null(n_meth) | is.null(n_unmeth)){
      stop("Missing beta or coverage values.\nU and M are not available either!", call. = FALSE)
    }else{
      if(verbose){
        message("Missing beta and coverage info.\nEstimating them from M and U values.")
      }

      return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                              beta = beta, M = n_meth, U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "cov := M+U", "beta := M/(M+U)")))
    }
  }else if(is.null(beta) & !is.null(cov)){ #If anyone of them is present (case-1: coverage available, estimate beta)
    if(all(is.null(n_meth), is.null(n_unmeth))){
      stop("Missing beta values but coverage info available.\nEither U or M are required for estimating beta values!", call. = FALSE)
    }else if(all(!is.null(n_meth), !is.null(n_unmeth))){
      message("Estimating beta values from M and U..")
      return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                              beta = beta, M = n_meth, U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "beta := M/(M+U")))
    }else if(!is.null(n_meth)){ #M available
      message("Estimating beta values from M and coverage..")
      return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                              beta = beta, M = n_meth, U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "beta := M/cov")))
    }else if(!is.null(n_unmeth)){ #U available
      message("Estimating beta values from U and coverage..")
      return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                              beta = beta, M = n_meth, U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, paste0("beta := 1- (U/cov)"))))
    }
  } else if(!is.null(beta) & is.null(cov)){ #If anyone of them is present (case-2: beta available, estimate coverage)
    if(all(is.null(n_meth), is.null(n_unmeth))){
      stop("Missing coverage info but beta values are available.\nU and M are required for estimating coverage values!", call. = FALSE)
    }else{
      if(verbose){
        message("Estimating coverage from M and U..")
      }

      return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                              beta = beta, M = n_meth, U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "cov := M+U")))
    }
  }else{
    if(verbose){
      message("All fields are present. Nice.")
    }

    return(list(col_idx = c(chr = chr, start = start, end = end, strand = strand,
                            beta = beta, cov = cov),
                fix_missing = NULL))
  }
}

#Read bedgraphs, and add missing info
.read_bdg = function(bdg, col_list = NULL, genome = NULL, verbose = TRUE){

  bdg_dat = data.table::fread(file = bdg, sep = "\t", colClasses = col_list$col_classes)
  colnames(bdg_dat)[col_list$col_idx] = names(col_list$col_idx)

  if(!is.null(col_list$fix_missing)){
    for(cmd in col_list$fix_missing){
      bdg_dat[,eval(parse(text=cmd))]
    }
  }

  bdg_dat = bdg_dat[,.(chr, start, beta, cov, strand)]
  bdg_dat[, chr := as.character(chr)][, start := as.integer(start)]
  data.table::setkey(x = bdg_dat, "chr", "start")
#browser()
  if(!is.null(genome)){
    if(!grepl("chr", bdg_dat[1,chr])){
      #warning("Dynamically determined that input bedgraph do not contain `chr` prefix. Will add it for convenience.", immediate. = TRUE)
      #Not very efficient way to do it! Might cause an issue. Check later and keep and eye.
      #Since data.table changes values in place wihout making copies, running gsub second time would be unnecessary. Check before gsub
      bdg_dat[, chr := paste0("chr", chr)]

      if(!grepl("chr", bdg_dat[1,chr])){
        bdg_dat[, chr := paste0("chr", chr)]
      }
    }

    missing_cpgs = genome[!bdg_dat[,list(chr, start)], on = c("chr", "start")]
    if(verbose){
      message(paste0("Missing ", nrow(missing_cpgs), " from: ", basename(bdg)))
    }
    if(nrow(missing_cpgs)>0){
      missing_cpgs[, width := NULL][, end := NULL][, beta := NA][, cov := 0]
      bdg_dat = data.table::rbindlist(list(bdg_dat, missing_cpgs), use.names = TRUE)
    }
    data.table::setkey(x = bdg_dat, "chr", "start")
  }
if(identical(bdg_dat[,list(chr, start)], genome[,list(chr, start)])){
  return(bdg_dat)}
  else{
    stop("Something went wrong with filling up the not covered CpG sites.")
  }
}


#Process all samples in one go (ideal for few number of samples)
# .vect_code = function(files, col_idx, col_data, ideal = FALSE, genome = NULL){
#   bdgs = lapply(files, .read_bdg, col_list = col_idx, genome = genome)
#   names(bdgs) = rownames(col_data)
#   gc()
#   cov_mat = data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE)
#   beta_mat = data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE)
#   colnames(cov_mat) = colnames(beta_mat) = rownames(col_data)
#   cov_mat = cbind(bdgs[[1]][,.(chr, start)], cov_mat)
#   beta_mat = cbind(bdgs[[1]][,.(chr, start)], beta_mat)
#   return(list(beta_matrix = beta_mat, cov_matrix = cov_mat))
# }

#Process samples in batches. Batches are processed in vectorized manner (ideal for large number of samples)
.vect_code_batch = function(files, col_idx, batch_size,  col_data = NULL, genome = NULL){
  batches = split(files, ceiling(seq_along(files)/batch_size))
  batches_samp_names = split(rownames(col_data), ceiling(seq_along(rownames(col_data))/batch_size))

  beta_mat_final = data.table::data.table()
  cov_mat_final = data.table::data.table()

  for(i in seq_along(batches)){
    #browser()
    message(paste0("Processing batch ",  i , " of ", length(batches)))
    batch_files = batches[[i]]
    samp_names = batches_samp_names[[i]]
    bdgs = lapply(batch_files, .read_bdg, col_list = col_idx, genome = genome)
    names(bdgs) = samp_names

    if(i == 1){
      cov_mat_final = data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE)
      beta_mat_final = data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE)
      colnames(cov_mat_final) = colnames(beta_mat_final) = samp_names
    }else{
      cov_mat = data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE)
      beta_mat = data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE)
      colnames(cov_mat) = colnames(beta_mat) = samp_names
      cov_mat_final = cbind(cov_mat_final, cov_mat)
      beta_mat_final = cbind(beta_mat_final, beta_mat)
      rm(cov_mat)
      rm(beta_mat)
      gc()
    }
  }

  return(list(beta_matrix = data.table::setDT(beta_mat_final), cov_matrix = data.table::setDT(cov_mat_final)))
}

#Use for loop for sample-by-sample processing, memory efficient, uses HDF5Array
.non_vect_code = function(files, col_idx, coldata, verbose = TRUE,  genome = NULL, h5temp=NULL, h5=NULL){

  if (h5){
    if(is.null(h5temp)){
      h5temp <- tempdir()}
    sink_counter <- 1
    while(any(c(paste0("M_sink_",sink_counter, ".h5"),
             paste0("cov_sink_",sink_counter, ".h5")) %in%  dir(h5temp))){
      sink_counter <- sink_counter+1

    }
    grid <- DelayedArray::RegularArrayGrid(
      refdim = c(nrow(genome), length(files)),
      spacings = c(nrow(genome), 1L))

    M_sink <- HDF5Array::HDF5RealizationSink(
      dim = c(nrow(genome), length(files)),
      dimnames = NULL,
      type = "double",
      filepath = file.path(h5temp, paste0("M_sink_",sink_counter, ".h5")),
      name = "M", chunkdim = HDF5Array::getHDF5DumpChunkDim(c(nrow(genome), length(files)), "double"),
      level = 6)
    cov_sink <- HDF5Array::HDF5RealizationSink(
      dim = c(nrow(genome), length(files)),
      dimnames = NULL,
      type = "integer",
      filepath = file.path(h5temp, paste0("cov_sink_",sink_counter, ".h5")),
      name = "cov", chunkdim = HDF5Array::getHDF5DumpChunkDim(c(nrow(genome), length(files)), "integer"),
      level = 6)
  } else {
    beta_mat = data.table::data.table()
    cov_mat = data.table::data.table()
  }
  if(h5){
    #browser()
    for(i in seq_along(files)){
      if(verbose){
        message("Processing: ", files[i])
      }
      b = .read_bdg(bdg = files[i], col_list = col_idx, genome = genome)
      DelayedArray::write_block_to_sink(block=as.matrix(b[, .(beta)]), viewport = grid[[i]], sink = M_sink)
      DelayedArray::write_block_to_sink(block=as.matrix(b[, .(cov)]), viewport = grid[[i]], sink = cov_sink)
      rm(b)
      gc()
    }
    return(list(beta_matrix = as(M_sink, "HDF5Array"), cov_matrix = as(cov_sink, "HDF5Array")))
  } else {
    for(i in seq_along(files)){
      if(verbose){
        message("Processing: ", files[i])
      }
      if(i == 1){
        b = .read_bdg(bdg = files[i], col_list = col_idx, genome = genome)
        beta_mat = b[,.(chr, start, beta)]
        cov_mat = b[,.(chr, start, cov)]
      }else{
        b = .read_bdg(bdg = files[i], col_list = col_idx, genome = genome)
        beta_mat = cbind(beta_mat, b[,.(beta)])
        cov_mat = cbind(cov_mat, b[,.(cov)])
      }
      colnames(beta_mat)[ncol(beta_mat)] = colnames(cov_mat)[ncol(cov_mat)] = rownames(coldata)[i]
    }
    return(list(beta_matrix = beta_mat, cov_matrix = cov_mat))

  }
}

