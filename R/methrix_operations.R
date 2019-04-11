#' Extracts methylation of regions of interest
#' @details Takes \code{\link{methrix}} object and summarizes regions
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param type data that should be summarized. Coule be methylation or coverage.
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min
#' @export
get_region_summary = function(m, regions = NULL, type = NULL, how = NULL){
  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }

  if(is.null(type)){
    stop("Please specify if you want to summarize methylation (meth) or coverage (cov).")
  }else if (type == "M") {
    message("Summarize methylation.")
    dat = get_matrix(m = m, type = "M", add_loci = TRUE)
  }else if (type == "C") {
    message("Summarize coverage.")
    dat = get_matrix(m = m, type = "C", add_loci = TRUE)
  }else{
      stop("Invalid input for summarization of regions.")
    }

  if(!is.null(regions)){
    message("Subsetting by genomic regions..")

    if(class(regions)[1] == "GRanges"){
      regions = as.data.frame(regions)
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions$chr = as.character(regions$chr)
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
      regions = regions[,.(chr, start, end)]
    }else if(class(regions)[1] == "data.table"){
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions = regions[,.(chr, start, end)]
      regions[, chr := as.character(chr)]
      regions[, start := as.numeric(start)]
      regions[, end := as.numeric(end)]
      regions$name <- paste0("chr",regions$chr, ".",regions$start, ".",regions$end)
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
    }else{
      stop("Invalid input class for regions. Must be a data.table or GRanges object")
    }

    dat[,end := start+1]
    overlap = data.table::foverlaps(x = dat, y = regions, type = "within", nomatch = NULL)

    if(nrow(overlap) == 0){
      stop("Subsetting resulted in zero entries")
    }
    if(is.null(how)) {
      stop("Please specify how regions should be specified.")
    }else if (how=="mean") {
      message("Summarizing by average..")
      output = overlap[, lapply(.SD, mean, na.rm=TRUE), by=name, .SDcols=rownames(colData(m))]
    }else if (how=="max") {
      message("Summarizing by maximum..")
      output = overlap[, lapply(.SD, max, na.rm=TRUE), by=name, .SDcols=rownames(colData(m))]
    }else if (how=="min") {
      message("Summarizing by minimum..")
      output = overlap[, lapply(.SD, min, na.rm=TRUE), by=name, .SDcols=rownames(colData(m))]
    }else if (how=="sum") {
      message("Summarizing by sum..")
      output = overlap[, lapply(.SD, sum, na.rm=TRUE), by=name, .SDcols=rownames(colData(m))]
    }else{
      stop("Invalid input for region summarization.")
    }
    
  output = cbind(regions[,c("chr", "start", "end"), with=F], mean[,rownames(colData(m)), with=F])

  return(output)
}
}

#' Order mathrix object by SD
#' @details Takes \code{\link{methrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{methrix}} object
#' @export
order_by_sd = function(m){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }else{
    row_order = order(matrixStats::rowSds(x = get_matrix(m = m, type = "M"), na.rm = TRUE), decreasing = TRUE)
    assay(m, i = 1) = assay(m, i = 1)[row_order,, drop = FALSE]
    assay(m, i = 2) = assay(m, i = 2)[row_order,, drop = FALSE]
    rowData(x = m) = S4Vectors::DataFrame(as.data.frame(x = rowData(x = m))[row_order,, drop = FALSE])
  }
  m
}

#--------------------------------------------------------------------------------------------------------------------------

#' Subsets \code{\link{methrix}} object based on given conditions.
#' @details Takes \code{\link{methrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{methrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param contigs chromosome names to subset by
#' @param samples sample names to subset by
#' @export
subset_methrix = function(m, regions = NULL, contigs = NULL, samples = NULL){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }

  m_dat = get_matrix(m = m, type = "M", add_loci = TRUE)
  c_dat = get_matrix(m = m, type = "C", add_loci = TRUE)

  if(!is.null(regions)){
    message("Subsetting by genomic regions..")

    if(class(regions)[1] == "GRanges"){
      regions = as.data.frame(regions)
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions$chr = as.character(regions$chr)
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
      regions = regions[,.(chr, start, end)]
    }else if(class(regions)[1] == "data.table"){
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions = regions[,.(chr, start, end)]
      regions[, chr := as.character(chr)]
      regions[, start := as.numeric(start)]
      regions[, end := as.numeric(end)]
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
    }else{
      stop("Invalid input class for regions. Must be a data.table or GRanges object")
    }

    m_dat[,end := start+1]
    m_dat = data.table::foverlaps(x = m_dat, y = regions, type = "within", nomatch = NULL)

    if(nrow(m_dat) == 0){
      stop("Subsetting resulted in zero entries")
    }

    m_dat[,i.end := NULL]
    m_dat[,i.start := NULL]

    c_dat[,end := start+1]
    c_dat = data.table::foverlaps(x = c_dat, y = regions, type = "within", nomatch = NULL)
    c_dat[,i.end := NULL]
    c_dat[,i.start := NULL]
  }

  if(!is.null(contigs)){

    message("Subsetting by contigs..")

    m_dat = m_dat[chr %in% contigs]
    if(nrow(m_dat) == 0){
      stop("Subsetting resulted in zero entries")
    }
    c_dat = c_dat[chr %in% contigs]
  }

  if(!is.null(samples)){

    message("Subsetting by samples..")

    samples = samples[samples %in% colnames(m_dat)]

    if(length(samples) == 0){
      stop("None of the samples are present in the object")
    }

    m_dat = m_dat[,c("chr", "start", "strand", samples), with = FALSE]
    c_dat = c_dat[,c("chr", "start", "strand", samples), with = FALSE]
  }


  range_idx = colnames(m_dat)[colnames(m_dat) %in% c("chr", "start", "end", "strand")]
  samp_names = colnames(m_dat)[!colnames(m_dat) %in% c("chr", "start", "end", "strand")]

  m = create_methrix(
    beta_mat = m_dat[,samp_names, with = FALSE],
    cov_mat = c_dat[,samp_names, with = FALSE],
    cpg_loci = m_dat[,range_idx, with = FALSE],
    is_hdf5 = is_h5(m),
    genome_name = m@metadata$genome,
    col_data = colData(m)[samp_names,, drop = FALSE],
    h5_dir = NULL, ref_cpg_dt = m@metadata$ref_CpG)

  return(m)
}

#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by coverage
#' @details Takes \code{\link{methrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{methrix}} object
#' @param cov_thr minimum coverage required to call a loci covered
#' @param min_samples At-least these many samples should have a loci with coverage >= \code{cov_thr}
#' @param n_threads number of threads to use. Default 4.
#' @export
coverage_filter = function(m, cov_thr = 1, min_samples = 1, n_threads = 4){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }else{

    cov_dat = get_matrix(m = m, type = "C")

    row_idx = parallel::mclapply(X = seq_len(nrow(cov_dat)), function(i){
      x = cov_dat[i,]
      length(x[x > cov_thr])
    }, mc.cores = n_threads)

    row_idx = unlist(row_idx)
    row_idx = row_idx >= min_samples

    message(paste0("Retained ", format(length(row_idx[row_idx]), big.mark = ","), " of ", format(nrow(cov_dat), big.mark = ","), " sites"))

    rm(cov_dat)
    gc()

    m = create_methrix(beta_mat = assay(m, i = 1)[row_idx,, drop = FALSE],
                       cov_mat = assay(m, i = 2)[row_idx,, drop = FALSE],
                       cpg_loci = data.table::as.data.table(as.data.frame(x = rowData(x = m))[row_idx,, drop = FALSE]),
                       is_hdf5 = is_h5(m), genome_name = m@metadata$genome,
                       col_data = colData(m), h5_dir = NULL, ref_cpg_dt = m@metadata$ref_CpG)

    return(m)
  }
}

#--------------------------------------------------------------------------------------------------------------------------

#' Extract methylation or coverage matrices
#' @details Takes \code{\link{methrix}} object and returns user specified \code{methylation} or \code{coverage} matrix
#' @param m \code{\link{methrix}} object
#' @param type can be \code{M} or \code{C}. Default "M"
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @export
#'
get_matrix= function(m, type = "M", add_loci = FALSE){

  type = match.arg(arg = type, choices = c("M", "C"))

  if(type == "M"){
    d = SummarizedExperiment::assay(x = m, i = 1)
  }else{
    d = SummarizedExperiment::assay(x = m, i = 2)
  }

  if(add_loci){
    d = as.data.frame(cbind(SummarizedExperiment::rowData(x = m), d))
    data.table::setDT(x = d)
  }

  d
}


#--------------------------------------------------------------------------------------------------------------------------

#' Convert methrix to bsseq object
#' @details Takes \code{\link{methrix}} object and returns a \code{\link{BSseq}} object
#' @param m \code{\link{methrix}} object
#' @export
#'
methrix2bsseq = function(m){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }

  n_samps = as.numeric(m@metadata$summary[1, Summary])
  warning("BSseq does not allow any uncovered loci (NA's). This will filter out all uncovered loci ", immediate. = TRUE)
  m_clean = methrix::coverage_filter(m = m, cov_thr = 1, min_samples = n_samps)


  #Thanks to Maxi for pointing out the bug related to M estimation
  #To-do: Find solution to avoid matrix multiplication (for small datasets it shouldn't affect)
  b = bsseq::BSseq(M = get_matrix(m_clean) * get_matrix(m_clean, type = "C"), Cov = get_matrix(m_clean, type = "C"),
                   pData = colData(x = m_clean), pos = rowData(x = m_clean)[,"start"], chr = rowData(x = m_clean)[,"chr"],
                   sampleNames = rownames(m_clean@colData))
  b
}

#--------------------------------------------------------------------------------------------------------------------------

#' Remove loci that are uncovered across all samples
#' @details Takes \code{\link{methrix}} object and removes loci that are uncovered across all samples
#' @param m \code{\link{methrix}} object
#' @export
#'
remove_uncovered = function(m){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }

  cov_dat = get_matrix(m = m, type = "C")

  row_idx = which(matrixStats::rowSums2(x = cov_dat) == 0)

  message(paste0("Removed ", format(length(row_idx), big.mark = ","),
                 " [", round(length(row_idx)/nrow(cov_dat) * 100, digits = 2), "%] uncovered loci of ",
                 format(nrow(cov_dat), big.mark = ","), " sites"))

  rm(cov_dat)
  gc()

  m = create_methrix(beta_mat = assay(m, i = 1)[-row_idx,, drop = FALSE],
                     cov_mat = assay(m, i = 2)[-row_idx,, drop = FALSE],
                     cpg_loci = data.table::as.data.table(as.data.frame(x = rowData(x = m))[-row_idx,, drop = FALSE]),
                     is_hdf5 = is_h5(m), genome_name = m@metadata$genome,
                     col_data = colData(m), h5_dir = NULL, ref_cpg_dt = m@metadata$ref_CpG)
  m
}

