#' Extract and summarize methylation or coverage info by regions of interest
#' @details Takes \code{\link{methrix}} object and summarizes regions
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param type matrix which needs to be summarized. Coule be `M`, `C`. Default "M"
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default "mean"
#' @param na_rm Remove NA's ? Default \code{TRUE}
#' @param verbose Default TRUE
#' @return a coverage or methylation matrix
#' @examples
#' data("methrix_data")
#' get_region_summary2(m = methrix_data, regions = data.table(chr = "chr21", start = 27867971, end =  27868103), type = "M", how = "mean")
#' @export
get_region_summary = function(m, regions = NULL, type = "M", how = "mean", na_rm = TRUE, verbose = TRUE){

  type = match.arg(arg = type, choices = c('M', 'C'))
  how = match.arg(arg = how, choices = c('mean', 'median', 'max', 'min', 'sum'))

  start_proc_time = proc.time()

  target_regions = cast_ranges(regions)
  #Add a unique id for every target range (i.e, rows)
  target_regions[, rid := paste0("rid_", 1:nrow(target_regions))]

  r_dat = data.table::as.data.table(rowData(x = m))
  r_dat[, chr := as.character(chr)]
  r_dat[, end := start + 1]
  data.table::setDT(x = r_dat, key = c("chr", "start", "end"))

  if(verbose){
    cat("-Checking for overlaps..\n")
  }

  overlap_indices = data.table::foverlaps(x = r_dat, y = target_regions, type = "any", nomatch = NULL, which = TRUE)

  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }

  overlap_indices[,yid := paste0("rid_", yid)]
  n_overlap_cpgs = overlap_indices[,.N,yid]
  colnames(n_overlap_cpgs) = c('rid', 'n_overlap_CpGs')

  #overlap_indices = split(overlap_indices, as.factor(as.character(overlap_indices$yid)))

  if (type == "M") {
    dat = get_matrix(m = m[overlap_indices$xid,], type = "M", add_loci = TRUE)
  }else if (type == "C") {
    dat = get_matrix(m = m[overlap_indices$xid,], type = "C", add_loci = TRUE)
  }

  if(nrow(overlap_indices) != nrow(dat)){
    stop("Something went wrong")
  }

  dat = cbind(overlap_indices, dat)

 #cat("-Summarizing overlaps..\n")
 if(how == "mean") {
   cat("-Summarizing by average\n")
   output = dat[, lapply(.SD, mean, na.rm = na_rm), by = yid, .SDcols = rownames(colData(m))]
 }else if (how == "median") {
   cat("-Summarizing by median\n")
   output = dat[, lapply(.SD, median, na.rm = na_rm), by = yid, .SDcols = rownames(colData(m))]
 }else if (how == "max") {
   cat("-Summarizing by maximum\n")
   output = dat[, lapply(.SD, max, na.rm = na_rm), by = yid, .SDcols = rownames(colData(m))]
 }else if (how == "min") {
   cat("-Summarizing by minimum\n")
   output = dat[, lapply(.SD, min, na.rm = na_rm), by = yid, .SDcols = rownames(colData(m))]
 }else if (how == "sum") {
   cat("-Summarizing by sum\n")
   output = dat[, lapply(.SD, sum, na.rm = na_rm), by = yid, .SDcols = rownames(colData(m))]
 }

  output = merge(target_regions, output, by.x = 'rid', by.y = 'yid', all.x = TRUE)
  output = merge(n_overlap_cpgs, output, by = 'rid')
  output[,rid := NULL]


  if(verbose){
    cat("-Done! Finished in:",data.table::timetaken(start_proc_time),"\n")
  }

  return(output)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Order mathrix object by SD
#' @details Takes \code{\link{methrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{methrix}} objectw
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data("methrix_data")
#' order_by_sd(m = methrix_data)
#' @export
order_by_sd = function(m){

  if(is_h5(m)){
    sds = DelayedMatrixStats::rowSds(x = get_matrix(m, type = "M"))
  }else{
    sds = matrixStats::rowSds(x = get_matrix(m, type = "M"))
  }
  row_order = order(sds, na.last = TRUE, decreasing = TRUE)
  m = m[row_order,]
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
#' data("methrix_data")
#' #Subset to chromosome 1
#' subset_methrix(methrix_data, contigs = "chr21")
#' @return An object of class \code{\link{methrix}}
#' @export
subset_methrix = function(m, regions = NULL, contigs = NULL, samples = NULL){

  r_dat <- data.table::as.data.table(rowData(m))

  if(!is.null(regions)){
    cat("-Subsetting by genomic regions\n")

    target_regions = cast_ranges(regions)

    r_dat[, end := start + 1]
    data.table::setDT(x = r_dat, key = c("chr", "start", "end"))
    overlaps = data.table::foverlaps(
      x = r_dat,
      y = target_regions,
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

    samples = which(rownames(colData(m)) %in% samples)
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
#' data("methrix_data")
#' #keep only CpGs which are covered by at-least 1 read across 3 samples
#' coverage_filter(m = methrix_data, cov_thr = 1, min_samples = 3)
#' @return An object of class \code{\link{methrix}}
#' @export
coverage_filter = function(m, cov_thr = 1, min_samples = 1){

  cov_dat = get_matrix(m = m, type = "C")

  res <- data.table::as.data.table(which(cov_dat >= cov_thr, arr.ind = TRUE))

  if (is_h5(m)){
  res <- res[,.(Count=(.N)), by=V1]
  row_idx = res$V1[res$Count >= min_samples]} else {
    res <- res[,.(Count=(.N)), by=row]
    row_idx = res$row[res$Count >= min_samples]
  }



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
#' data("methrix_data")
#' #Get methylation matrix
#' get_matrix(m = methrix_data, type = "M")
#' #Get methylation matrix along with loci
#' get_matrix(m = methrix_data, type = "M", add_loci = TRUE)
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
#' data("methrix_data")
#' methrix2bsseq(m = methrix_data)
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
#' data("methrix_data")
#' remove_uncovered(m = methrix_data)
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
#' data("methrix_data")
#' region_filter(m = methrix_data, regions = data.table(chr = "chr21", start = 27867971, end =  27868103))
#' @export
region_filter = function(m, regions){


  target_regions = cast_ranges(regions)

  current_regions <-  data.table::as.data.table(colData(m))
  current_regions[,end := start+1]
  data.table::setDT(x = current_regions, key = c("chr", "start", "end"))
  overlap = data.table::foverlaps(x = current_regions, y = target_regions, type = "within", nomatch = NULL, which = TRUE)

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
#' Combine methrix objects
#' @details Takes two \code{\link{methrix}} objects and combines them row- or column-wise
#' @param m1 \code{\link{methrix}} object
#' @param m1 \code{\link{methrix}} object
#' @param by The direction of combine. "column" (cbind) combines samples with same regions, "row" combines different regions,
#' e.g. different chromosomes.
#' @return An object of class \code{\link{methrix}}
#' @export
#'
combine_methrix = function(m1, m2, by = c("row", "col")){

  by = match.arg(arg = by, choices = c("row", "col"), several.ok = FALSE)

  if (by=="row"){
    if (!(all(rownames(m1@colData)==rownames(m2@colData)))){
      stop("You have different samples in your dataset. You need the same samples in your datasets. ")
    } else {
      m <- rbind(m1, m2)
    }
    if (any(duplicated((as.data.table(m@elementMetadata))))){
      stop("There are overlapping regions in your datasets. This function only takes distinct objects. ")
    }
  }
  if (by=="col"){
    if (any(rownames(m1@colData) %in% rownames(m2@colData))){
      stop("You have the same samples in your datasets. You need different samples for this merging.  ")
    } else if (!identical(m1@elementMetadata, m2@elementMetadata)){
      stop("You have to have the same regions in your datasets. ")
    } else {
      m <- cbind(m1, m2)
    }
  }
  gc()
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
#' Estimate descriptive statistics
#' @details Calculate descriptive statistics
#' @param m \code{\link{methrix}} object
#' @param per_chr Estimate stats per chromosome. Default FALSE
#' @param skip_cov Default \code{TRUE}
#' @seealso \code{\link{plot_stats}}
#' @examples
#' data("methrix_data")
#' get_stats(methrix_data)
#' @export
get_stats = function(m, per_chr = FALSE, skip_cov = TRUE){

  #m_sub = remove_uncovered(m = m)

  if(per_chr){
    row_df = data.frame(rowData(x = m))
    row_chrs = names(table(row_df$chr))
    cat("-Processing", length(row_chrs), "chromosomes")

    stats = lapply(row_chrs, function(x){
      if(is_h5(m)){
        x_idx = which(row_df$chr == x)
        mean_meth = DelayedMatrixStats::colMeans2(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        median_meth = DelayedMatrixStats::colMedians(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        sd_meth = DelayedMatrixStats::colSds(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        gc(verbose = FALSE)

        mean_cov = median_cov = sd_cov = NA
        if(!skip_cov){
          cov_mat = get_matrix(m = m[x_idx,], type = "C")
          #Set uncovered loci to NA (zeros can affect mean and median)
          cov_mat[cov_mat == 0] = NA
          mean_cov = DelayedMatrixStats::colMeans2(x = cov_mat, na.rm = TRUE)
          median_cov = DelayedMatrixStats::colMedians(cov_mat, na.rm = TRUE)
          sd_cov = DelayedMatrixStats::colSds(cov_mat, na.rm = TRUE)

          rm(cov_mat)
          gc(verbose = FALSE)
        }

      }else{
        x_idx = which(row_df$chr == x)
        mean_meth = matrixStats::colMeans2(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        median_meth = matrixStats::colMedians(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        sd_meth = matrixStats::colSds(x = get_matrix(m = m[x_idx,], type = "M"), na.rm = TRUE)
        gc(verbose = FALSE)

        mean_cov = median_cov = sc_cov = NA
        if(!skip_cov){
          cov_mat = get_matrix(m = m[x_idx,], type = "C")
          #Set uncovered loci to NA (zeros can affect mean and median)
          cov_mat[cov_mat == 0] = NA
          mean_cov = matrixStats::colMeans2(x = cov_mat, na.rm = TRUE)
          median_cov = matrixStats::colMedians(cov_mat, na.rm = TRUE)
          sd_cov = matrixStats::colSds(cov_mat, na.rm = TRUE)

          rm(cov_mat)
          gc(verbose = FALSE)
        }
      }

      chr_stat = data.table(
        Sample_Name = rownames(colData(x = m)),
        mean_meth = mean_meth,
        median_meth = median_meth,
        sd_meth = sd_meth,
        mean_cov = mean_cov,
        median_cov = median_cov,
        sd_cov = sd_cov,
        stringsAsFactors = FALSE
      )

      chr_stat
    })

    names(stats) = row_chrs
    stats = data.table::rbindlist(l = stats, use.names = TRUE, fill = TRUE, idcol = "Chromosome")

  }else{

    if(is_h5(m)){
      mean_meth = DelayedMatrixStats::colMeans2(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      median_meth = DelayedMatrixStats::colMeans2(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      sd_meth = DelayedMatrixStats::colSds(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      gc(verbose = FALSE)

      if(!skip_cov){
        cov_mat = get_matrix(m = m, type = "C")
        #Set uncovered loci to NA (zeros can affect mean and median)
        cov_mat[cov_mat == 0] = NA
        mean_cov = DelayedMatrixStats::colMeans2(x = cov_mat, na.rm = TRUE)
        median_cov = DelayedMatrixStats::colMedians(cov_mat, na.rm = TRUE)
        sd_cov = DelayedMatrixStats::colSds(cov_mat, na.rm = TRUE)
        rm(cov_mat)
        gc(verbose = FALSE)
      }
    }else{
      mean_meth = matrixStats::colMeans2(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      median_meth = matrixStats::colMedians(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      sd_meth = matrixStats::colSds(x = get_matrix(m = m, type = "M"), na.rm = TRUE)
      gc(verbose = FALSE)

      mean_cov = median_cov = sc_cov = NA
      if(!skip_cov){
        cov_mat = get_matrix(m = m, type = "C")
        #Set uncovered loci to NA (zeros can affect mean and median)
        cov_mat[cov_mat == 0] = NA
        mean_cov = matrixStats::colMeans2(x = cov_mat, na.rm = TRUE)
        median_cov = matrixStats::colMedians(cov_mat, na.rm = TRUE)
        sd_cov = matrixStats::colSds(cov_mat, na.rm = TRUE)
        rm(cov_mat)
        gc(verbose = FALSE)
      }
    }

    stats = data.table(
      Sample_Name = rownames(colData(x = m)),
      mean_meth = mean_meth,
      median_meth = median_meth,
      sd_meth = sd_meth,
      mean_cov = mean_cov,
      median_cov = median_cov,
      sd_cov = sd_cov,
      stringsAsFactors = FALSE
    )
  }

  stats
}

#--------------------------------------------------------------------------------------------------------------------------
#' Extract CpGs covered per chromosome information
#' @details Takes \code{\link{methrix}} object and returns a table of CpGs covered per chromosome
#' @param m \code{\link{methrix}} object
#' @examples
#' data("methrix_data")
#' get_chr_summary(methrix_data)
#' @export
get_chr_summary = function(m = NULL){
  chr_tbl = table(rowData(x = m)[,"chr"])
  chr_tbl = data.table::data.table(chr_tbl)
  names(chr_tbl) = c("chr", "n_CpG")
  chr_tbl = merge(chr_tbl, m@metadata$ref_CpG, by = 'chr')
  names(chr_tbl) = c('chr', 'n_CpG', 'n_ref_CpG')
  chr_tbl[,percent_ref_CpG_covered := round(n_CpG/n_ref_CpG, digits = 4) * 100]
  chr_tbl
}


#--------------------------------------------------------------------------------------------------------------------------
#' Saves HDF5 methrix object
#' @details Takes \code{\link{methrix}} object and saves it
#' @param m \code{\link{methrix}} object
#' @param dir The directory to use. Created, if not existing.
#' @param replace Should it overwrite the pre-existing data? FALSE by default.
#' @param ... Parameters to pass to saveHDF5SummarizedExperiment
#' @examples
#' data("methrix_data")
#' save_HDF5_methrix(methrix_data, dir="/my_methrix_folder", replace=T)
#' @export
save_HDF5_methrix = function(m=NULL, dir="", replace=FALSE, ...){

  if (class(m)=="methrix" && is_h5(m)){
    HDF5Array::saveHDF5SummarizedExperiment(x=m, dir=dir, replace = replace, ...)
  } else {
    stop("The object is not a mthrix object or not in an HDF5 format. ")
  }

}

#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 methrix object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{methrix}} object and loads it
#' @param dir The directory to read in from.
#' @param ... Parameters to pass to loadHDF5SummarizedExperiment
#' @return An object of class \code{\link{methrix}}
#' @examples
#' load_HDF5_methrix(dir="/my_methrix_folder")
#' @export
load_HDF5_methrix = function(dir="", ...){

  m <- HDF5Array::loadHDF5SummarizedExperiment(dir=dir,  ...)
  m <- as(m, "methrix")
  return(m)
}
