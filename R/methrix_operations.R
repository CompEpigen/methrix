#' Order mathrix object by SD
#' @details Takes \code{\link{methrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{methrix}} object
#' @export
order_by_sd = function(m){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }else{
    row_order = order(matrixStats::rowSds(x = assay(x, 1)), decreasing = TRUE)
    assay(m, i = 1) = assay(m, i = 1)[row_order,, drop = FALSE]
    assay(m, i = 2) = assay(m, i = 2)[row_order,, drop = FALSE]
    rowData(x = m) = S4Vectors::DataFrame(as.data.frame(x = rowData(x = x))[row_order,, drop = FALSE])
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

  m_dat = get_matrix(m = m, type = "M", add_lcoi = TRUE)
  c_dat = get_matrix(m = m, type = "C", add_lcoi = TRUE)

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
    h5_dir = NULL
  )

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
filter_methrix = function(m, cov_thr = 1, min_samples = 1, n_threads = 4){

  if(is_h5(m)){
    stop("This function only supports non HDF5 matrices for now.")
  }else{
    cov_dat = assay(x = m, i = 2)

    row_idx = parallel::mclapply(X = seq_len(nrow(cov_dat)), function(i){
      x = cov_dat[i,]
      length(x[x > cov_thr])
    }, mc.cores = n_threads)

    row_idx = unlist(row_idx)
    row_idx = row_idx >= min_samples

    message(paste0("Retained ", length(row_idx[row_idx]), " of ", nrow(cov_dat), " sites"))

    rm(cov_dat)

    m = create_methrix(beta_mat = assay(m, i = 1)[row_idx,, drop = FALSE],
                       cov_mat = assay(m, i = 2)[row_idx,, drop = FALSE],
                       cpg_loci = S4Vectors::DataFrame(as.data.frame(x = rowData(x = m))[row_idx,, drop = FALSE]),
                       is_hdf5 = is_h5(m), genome_name = m@metadata$genome,
                       col_data = colData(m), h5_dir = NULL)

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
get_matrix= function(m, type = "M", add_lcoi = FALSE){

  type = match.arg(arg = type, choices = c("M", "C"))

  if(type == "M"){
    d = SummarizedExperiment::assay(x = m, i = 1)
  }else{
    d = SummarizedExperiment::assay(x = m, i = 2)
  }

  if(add_lcoi){
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

  b = bsseq::BSseq(M = get_matrix(m), Cov = get_matrix(m, type = "C"),
                   pData = colData(x = m), pos = rowData(x = m)[,"start"], chr = rowData(x = m)[,"chr"],
                   sampleNames = rownames(m@colData))
  b
}
