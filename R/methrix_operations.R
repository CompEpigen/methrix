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
