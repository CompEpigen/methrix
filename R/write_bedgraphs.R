#' Writes bedGraphs from methrix object
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved. If \code{NULL} creats a \code{tempdir}
#' @param n_thr Default 4. Only used if \code{plot_beta_dist} is TRUE
#' @param make_bigWig Default FALSE. Whether to convert bedGraph to BigWig
#' @param bedtools_path Default NULL. Path to bedtools. Default looks in $PATH. Only applicable if \code{make_bigWig} = TRUE
#' @param bedGraphToBigWig_path Default NULL. Path to bedGraphToBigWig Default looks in $PATH. Only applicable if \code{make_bigWig} = TRUE
#' @param bedClip_path Default NULL. Path to bedClip Default looks in $PATH. Only applicable if \code{make_bigWig} = TRUE
#' @export

write_bedgraphs = function(m, output_dir = NULL, rm_NA = TRUE, force = FALSE, n_thr = 4, make_BigWig = FALSE,
                           bedtools_path = NULL, bedGraphToBigWig_path = NULL, bedClip_path = NULL){


  if(!dir.exists(output_dir)){
    dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  mat = methrix::get_matrix(m = m, type = "M", add_loci = TRUE)
  mat[,end := start + 1]
  mat = mat[,c("chr", "start", "end", "strand", rownames(colData(x = m))), with = FALSE]

  if(is.null(output_dir)){
    output_dir = getwd()
  }

  cat("----------------------\n")
  cat("*Writing bedGraphs:\n")
  op_bdgs = lapply(seq_len(nrow(colData(m))), function(i){
                mat_i = mat[,c(1:3, i+4), with = FALSE]
                if(rm_NA){
                  mat_i = mat_i[complete.cases(mat_i),,]
                }
                op_bdg = paste0(output_dir, "/", rownames(colData(m))[i], ".bedGraph")
                if(file.exists(op_bdg)){
                  if(force){
                    cat(paste0("**Writing ", rownames(colData(m))[i], "\n"))
                    colnames(mat_i) = paste0("V", 1:ncol(mat_i))
                    #https://github.com/Rdatatable/data.table/issues/2020
                    #Remove below two lines once #2020 is fixed
                    mat_i$V2 = format(as.integer(mat_i$V2), scientific = FALSE)
                    mat_i$V3 = format(as.integer(mat_i$V3), scientific = FALSE)
                    data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", col.names = FALSE, nThread = n_thr)
                    paste0(output_dir, "/", rownames(colData(m))[i], ".bedGraph")
                  }else{
                    cat(paste0("**File ", basename(op_bdg), " already exists. Skipped re-writing\n"))
                  }
                }else{
                  cat(paste0("**Writing ", rownames(colData(m))[i], "\n"))
                  colnames(mat_i) = paste0("V", 1:ncol(mat_i))
                  #https://github.com/Rdatatable/data.table/issues/2020
                  #Remove below two lines once #2020 is fixed
                  mat_i$V2 = format(as.integer(mat_i$V2), scientific = FALSE)
                  mat_i$V3 = format(as.integer(mat_i$V3), scientific = FALSE)
                  data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", col.names = FALSE, nThread = n_thr)
                  paste0(output_dir, "/", rownames(colData(m))[i], ".bedGraph")
                }
                op_bdg
              })
  cat("----------------------\n")

  if(make_BigWig){
    if(is.na(m@metadata$chrom_sizes)){
      stop("meta data contains no chrom_sizes slot.")
    }
    cat("*Converting bedGraphs to bigWig:\n")
    cat("**Checking for bed* binaries\n")
    bt = check_bedtools(path = bedtools_path)
    bb = check_bedGraphToBigWig(path = bedGraphToBigWig_path)
    bc = check_bedClip(path = bedClip_path)
    cat("**All good\n")

    chr_lens = data.table::copy(x = m@metadata$chrom_sizes)
    #chr_lens[,start := 1]
    chrom_sizes = paste0(output_dir, "/", "chrom.sizes")
    data.table::fwrite(x = chr_lens, file = chrom_sizes, sep = "\t", col.names = FALSE)

    bdg2bw = system.file('extdata', 'bdg2bw.sh', package = 'methrix')
    lapply(op_bdgs, function(b){
      cat(paste0("***Converting ", basename(b), "\n"))
      cmd = paste("bash", bdg2bw, b, chrom_sizes, bt, bc, bb)
      system(command = cmd, intern = TRUE)
    })
    cat("----------------------\n")
  }
}
