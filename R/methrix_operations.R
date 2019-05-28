#' Extract and summarize methylation or coverage info by regions of interest
#' @details Takes \code{\link{methrix}} object and summarizes regions
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param type matrix which needs to be summarized. Coule be `M`, `MR`(Methylated Reads) or `C`. Default "M"
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default "mean"
#' @param na_rm Remove NA's ? Default \code{TRUE}
#' @return a coverage or methylation matrix
#' @examples
#' data("mm9_bsmap")
#' get_region_summary(m = mm9_bsmap, regions = data.table(chr = "chr1", start = 3020220, end = 3209171), type = "M", how = "mean")
#' @export
get_region_summary = function(m, regions = NULL, type = "M", how = "mean", na_rm = TRUE){

  type = match.arg(arg = type, choices = c('M', 'C','MR'))
  how = match.arg(arg = how, choices = c('mean', 'sum', 'max', 'min'))

  if(is(regions[1], "GRanges")){
    regions = as.data.frame(regions)
    colnames(regions)[1:3] = c("chr", "start", "end")
    regions$chr = as.character(regions$chr)
    regions_work = data.table::copy(regions)
    data.table::setDT(x = regions_work, key = c("chr", "start", "end"))
    regions_work = regions_work[,.(chr, start, end)]
  }else if(is(regions[1], 'data.table')){
    regions_work = data.table::copy(regions)
    data.table::setDT(x = regions_work)
    colnames(regions_work)[1:3] = c("chr", "start", "end")
    regions_work = regions_work[,.(chr, start, end)]
    regions_work[, chr := as.character(chr)]
    regions_work[, start := as.numeric(start)]
    regions_work[, end := as.numeric(end)]
  }else{
    stop("Invalid input class for regions. Must be a data.table, data.frame or GRanges object")
  }

  regions_work[, id := paste0(chr, ":", start, "-", end)]
  data.table::setDT(x = regions_work, key = c("chr", "start", "end"))

  if (type == "M") {
    dat = get_matrix(m = m, type = "M", add_loci = TRUE)
  }else if (type == "C") {
    dat = get_matrix(m = m, type = "C", add_loci = TRUE)
  }else if (type == "MR") {
    dat_C = get_matrix(m = m, type = "C", add_loci = FALSE)
    dat_M = get_matrix(m = m, type = "M", add_loci = FALSE)
    reg = get_matrix(m = m, type = "M", add_loci = TRUE)[,1:3]
    dat = cbind(reg,dat_C*dat_M)
    }

  dat[,end := start+1]
  #region_overlap <- unique(data.table::foverlaps(x = dat, y = regions_work, type = "any", nomatch = NULL, which=T)$yid)
  overlap = data.table::foverlaps(x = dat, y = regions_work, type = "any", nomatch = NULL)


  if(nrow(overlap) == 0){
    stop("Subsetting resulted in zero entries")
  }

  if(how == "mean") {
    cat("-Summarizing by average\n")
    output = overlap[, lapply(.SD, mean, na.rm = na_rm), by = id, .SDcols = rownames(colData(m))]
  }else if (how == "max") {
    cat("-Summarizing by maximum\n")
    output = overlap[, lapply(.SD, max, na.rm = na_rm), by = id, .SDcols = rownames(colData(m))]
  }else if (how == "min") {
    cat("-Summarizing by minimum\n")
    output = overlap[, lapply(.SD, min, na.rm = na_rm), by = id, .SDcols = rownames(colData(m))]
  }else if (how == "sum") {
    cat("-Summarizing by sum\n")
    output = overlap[, lapply(.SD, sum, na.rm = na_rm), by = id, .SDcols = rownames(colData(m))]
  }

  output <- merge(regions_work[,list(id)],output, by="id", all=T, sort=F)
  output <- output[order(order(regions$chr, regions$start, regions$end)),]
  #output = cbind(regions[,c("chr", "start", "end"), with = FALSE], mean[,rownames(colData(m)), with = FALSE])
  return(output)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Order mathrix object by SD
#' @details Takes \code{\link{methrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data("mm9_bsmap")
#' order_by_sd(m = mm9_bsmap)
#' @export
order_by_sd = function(m){

  if(is_h5(m)){
    row_order = order(DelayedMatrixStats::rowSds(x = get_matrix(m = m, type = "M"), na.rm = TRUE), decreasing = TRUE)
  }else{
    row_order = order(matrixStats::rowSds(x = get_matrix(m = m, type = "M"), na.rm = TRUE), decreasing = TRUE)
    m = m[row_order,]
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
#' @examples
#' data("mm9_bsmap")
#' #Subset to chromosome 1
#' subset_methrix(mm9_bsmap, contigs = "chr1")
#' @return An object of class \code{\link{methrix}}
#' @export
subset_methrix = function(m, regions = NULL, contigs = NULL, samples = NULL){



  #m_dat = get_matrix(m = m, type = "M", add_loci = TRUE)
  #c_dat = get_matrix(m = m, type = "C", add_loci = TRUE)
  r_dat <- as.data.table(m@elementMetadata)
  if(!is.null(regions)){
    cat("-Subsetting by genomic regions\n")

    if(is(regions[1], "GRanges")){
      regions = as.data.frame(regions)
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions$chr = as.character(regions$chr)
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
      regions = regions[,.(chr, start, end)]
    }else if(is(regions[1], "data.table")){
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions = regions[,.(chr, start, end)]
      regions[, chr := as.character(chr)]
      regions[, start := as.numeric(start)]
      regions[, end := as.numeric(end)]
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
    }else{
      stop("Invalid input class for regions. Must be a data.table or GRanges object")
    }

    r_dat[, end := start + 1]
    data.table::setDT(x = r_dat, key = c("chr", "start", "end"))
    overlaps = data.table::foverlaps(
      x = r_dat,
      y = regions,
      type = "within",
      nomatch = NULL,
      which = TRUE
    )
    if (nrow(overlaps) == 0) {
      stop("Subsetting resulted in zero entries")
    }

    m <- m[overlaps$xid, ]
  }

  if(!is.null(contigs)) {
    cat("-Subsetting by contigs\n")
    selected_rows <- which(r_dat$chr %in% contigs)

    if (length(selected_rows) == 0) {
      stop("Subsetting resulted in zero entries")
    }
    m <- m[selected_rows, ]
  }

  if(!is.null(samples)) {
    cat("Subsetting by samples\n")

    samples = which(rownames(m@colData) %in% samples)
    if (length(samples) == 0) {
      stop("None of the samples are present in the object")
    }

    m <- m[, samples]
  }


  if(is_h5(m)){
    n_non_covered = length(which(DelayedMatrixStats::rowSums2(x = m@assays[["cov"]]) == 0))
  } else {
    n_non_covered = length(which(matrixStats::rowSums2(x = m@assays[["cov"]]) == 0))}


  se_summary = data.table::data.table(ID = c("n_samples", "n_CpGs", "n_uncovered", "n_chromosomes", "Reference_Build", "is_H5"),
                                      Summary = c(ncol(m), format(nrow(m), big.mark = ","),
                                                  n_non_covered, length(unique(m@elementMetadata$chr)), m@metadata$genome, m@metadata$is_h5))
  m@metadata$summary <- se_summary


  return(m)
}

#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by coverage
#' @details Takes \code{\link{methrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{methrix}} object
#' @param cov_thr minimum coverage required to call a loci covered
#' @param min_samples At-least these many samples should have a loci with coverage >= \code{cov_thr}
#' @importFrom methods is as new
#' @examples
#' data("mm9_bsmap")
#' #keep only CpGs which are covered by at-least 1 read across 3 samples
#' coverage_filter(m = mm9_bsmap, cov_thr = 1, min_samples = 3, n_threads = 1)
#' @return An object of class \code{\link{methrix}}
#' @export
coverage_filter = function(m, cov_thr = 1, min_samples = 1){

    res <- as.data.table(which(cov_dat <= cov_thr, arr.ind = T))
    res <- res[,.(Count=(.N)), by=V1]
    row_idx = res$V1[res$Count > min_samples]


    cat(paste0("-Retained ", format(length(row_idx[row_idx]), big.mark = ","), " of ", format(nrow(cov_dat), big.mark = ","), " sites\n"))

    rm(cov_dat)
    gc()

    m = m[row_idx,]

    if(is_h5(m)){
      n_non_covered = length(which(DelayedMatrixStats::rowSums2(x = m@assays[["cov"]]) == 0))
    } else {
      n_non_covered = length(which(matrixStats::rowSums2(x = m@assays[["cov"]]) == 0))}


    m@metadata$summary = data.table::data.table(ID = c("n_samples", "n_CpGs", "n_uncovered", "n_chromosomes", "Reference_Build", "is_H5"),
                                                Summary = c(ncol(m), format(nrow(m), big.mark = ","),
                                                            n_non_covered, length(unique(m@elementMetadata$chr)), m@metadata$genome, m@metadata$is_h5))

    return(m)

}

#--------------------------------------------------------------------------------------------------------------------------

#' Extract methylation or coverage matrices
#' @details Takes \code{\link{methrix}} object and returns user specified \code{methylation} or \code{coverage} matrix
#' @param m \code{\link{methrix}} object
#' @param type can be \code{M} or \code{C}. Default "M"
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @return Coverage or Methylation matrix
#' @examples
#' data("mm9_bsmap")
#' #Get methylation matrix
#' get_matrix(m = mm9_bsmap, type = "M")
#' #Get methylation matrix along with loci
#' get_matrix(m = mm9_bsmap, type = "M", add_loci = TRUE)
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
    if (is_h5(m)){
      d = as.data.frame(cbind(SummarizedExperiment::rowData(x = m), as.data.frame(d)))
      } else {
    d = as.data.frame(cbind(SummarizedExperiment::rowData(x = m), d))
    }
    data.table::setDT(x = d)
  }
  d
}

#--------------------------------------------------------------------------------------------------------------------------

#' Convert methrix to bsseq object
#' @details Takes \code{\link{methrix}} object and returns a \code{\link{BSseq}} object
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{\link{BSseq}}
#' @examples
#' \dontrun{
#' data("mm9_bsmap")
#' methrix2bsseq(m = mm9_bsmap)
#' }
#' @export
#'
methrix2bsseq = function(m){

  n_samps = nrow(SummarizedExperiment::colData(x = m))
  M_clean <- get_matrix(m) * get_matrix(m, type = "C")
  M_clean[is.na(M_clean)] <- 0

  #Thanks to Maxi for pointing out the bug related to M estimation
  #To-do: Find solution to avoid matrix multiplication (for small datasets it shouldn't affect)
  b = bsseq::BSseq(M = M_clean,
                   Cov = get_matrix(m, type = "C"),
                   pData = colData(x = m),
                   pos = rowData(x = m)[,"start"],
                   chr = rowData(x = m)[,"chr"],
                   sampleNames = rownames(m@colData))
  b
}

#--------------------------------------------------------------------------------------------------------------------------

#' Remove loci that are uncovered across all samples
#' @details Takes \code{\link{methrix}} object and removes loci that are uncovered across all samples
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data("mm9_bsmap")
#' remove_uncovered(m = mm9_bsmap)
#' @export
#'
remove_uncovered = function(m){

  if(is_h5(m)){
    cov_dat = m@assays[["cov"]]

    row_idx = which(DelayedMatrixStats::rowSums2(x = cov_dat) == 0)

  } else {

    cov_dat = get_matrix(m = m, type = "C")

    row_idx = which(matrixStats::rowSums2(x = cov_dat, na.rm = TRUE) == 0)
  }


  cat(paste0("-Removed ", format(length(row_idx), big.mark = ","),
                 " [", round(length(row_idx)/nrow(cov_dat) * 100, digits = 2), "%] uncovered loci of ",
                 format(nrow(cov_dat), big.mark = ","), " sites\n"))

  rm(cov_dat)
  gc()

  m = m[-row_idx,]

  if(is_h5(m)){
    n_non_covered = length(which(DelayedMatrixStats::rowSums2(x = m@assays[["cov"]]) == 0))
  } else {
    n_non_covered = length(which(matrixStats::rowSums2(x = m@assays[["cov"]]) == 0))
  }


  m@metadata$summary = data.table::data.table(ID = c("n_samples", "n_CpGs", "n_uncovered", "n_chromosomes", "Reference_Build", "is_H5"),
                                              Summary = c(ncol(m), format(nrow(m), big.mark = ","),
                                                          n_non_covered, length(unique(m@elementMetadata$chr)), m@metadata$genome, m@metadata$is_h5))

  m
}

#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by region
#' @details Takes \code{\link{methrix}} object and filters CpGs based on supplied regions in data.table or GRanges format
#' @param m \code{\link{methrix}} object
#' @param regions genomic regions to filter-out. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data("mm9_bsmap")
#' region_filter(m = mm9_bsmap, regions = data.table(chr = "chr1", start = 3020220, end = 3209171))
#' @export
region_filter = function(m, regions){


    if(is(regions[1], "GRanges")){
      regions = as.data.frame(regions)
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions$chr = as.character(regions$chr)
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
      regions = regions[,.(chr, start, end)]
    }else if(is(regions[1], "data.table")){
      colnames(regions)[1:3] = c("chr", "start", "end")
      regions = regions[,.(chr, start, end)]
      regions[, chr := as.character(chr)]
      regions[, start := as.numeric(start)]
      regions[, end := as.numeric(end)]
      data.table::setDT(x = regions, key = c("chr", "start", "end"))
    }else{
      stop("Invalid input class for regions. Must be a data.table or GRanges object")
    }

    current_regions <-  as.data.table(m@elementMetadata)
    current_regions[,end := start+1]
    data.table::setDT(x = current_regions, key = c("chr", "start", "end"))
    overlap = data.table::foverlaps(x = current_regions, y = regions, type = "within", nomatch = NULL, which = TRUE)

    if(nrow(overlap) == 0){
      stop("No CpGs found within the query intervals. Nothing to remove.")
    }

    cat(paste0("-Removed ", format(nrow(overlap), big.mark = ','), " CpGs\n"))

    m <- m[-overlap$xid,]
    if(is_h5(m)){
      n_non_covered = length(which(DelayedMatrixStats::rowSums2(x = m@assays[["cov"]]) == 0))
    } else {
      n_non_covered = length(which(matrixStats::rowSums2(x = m@assays[["cov"]]) == 0))}

    se_summary = data.table::data.table(ID = c("n_samples", "n_CpGs", "n_uncovered", "n_chromosomes", "Reference_Build", "is_H5"),
                                        Summary = c(ncol(m), format(nrow(m), big.mark = ","),
                                                    n_non_covered, nrow(current_regions[overlap$xid,.N,chr]), m@metadata$genome, m@metadata$is_h5))
    m@metadata$summary <- se_summary
    return(m)

}

#--------------------------------------------------------------------------------------------------------------------------
