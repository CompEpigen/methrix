#' Extract and summarize methylation or coverage info by regions of interest
#' @details Takes \code{\link{methrix}} object and summarizes regions
#' @param m \code{\link{methrix}} object
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param type matrix which needs to be summarized. Coule be `M`, `C`. Default 'M'
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default 'mean'
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description,
#' see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' @param na_rm Remove NA's? Default \code{TRUE}
#' @param elementMetadata.col columns in \code{\link{methrix}}@elementMetadata which needs to be summarised. Default = NULL.
#' @param n_chunks Number of chunks to split the \code{\link{methrix}} object in case it is very large. Default = 1.
#' @param n_cores Number of parallel instances. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @param verbose Default TRUE
#' @return a coverage or methylation matrix
#' @examples
#' data('methrix_data')
#' get_region_summary(m = methrix_data,
#' regions = data.table(chr = 'chr21', start = 27867971, end =  27868103),
#' type = 'M', how = 'mean')
#' @export
get_region_summary = function(m, regions = NULL, type = "M", how = "mean", overlap_type = "within", 
                              na_rm = TRUE, elementMetadata.col = NULL, verbose = TRUE, n_chunks=1, n_cores=1){

    yid  <- NULL
    if (n_cores > n_chunks){
        n_chunks <- n_cores
        message("n_cores should be set to be less than or equal to n_chunks.","\n","n_chunks has been set to be equal to n_cores = ", n_cores)
    }


    type = match.arg(arg = type, choices = c('M', 'C'))
    how = match.arg(arg = how, choices = c('mean', 'median', 'max', 'min', 'sum', 'sd'))

    start_proc_time = proc.time()


    if("data.table" %in% class(regions)){
        regions <- GenomicRanges::makeGRangesFromDataFrame(df = regions)
    }
    

    target_regions = (regions)
    #Add a unique id for every target range (i.e, rows)
    target_regions@elementMetadata$rid <- paste0("rid_", 1:length(target_regions))

    r_dat = as.data.frame(rowData(x = m))
    r_dat$seqnames <- as.character(r_dat$chr)
    r_dat$chr <- NULL
    
    if(is.null(r_dat$end)){
        r_dat$end <- r_dat$start + 1
    }
    
    r_dat <- GenomicRanges::makeGRangesFromDataFrame(r_dat, keep.extra.columns = FALSE)


    if(!all(elementMetadata.col %in% colnames(m@elementMetadata))) stop("variables provided to elementMetadata.col not correct")

    if(verbose){
        message("-Checking for overlaps..\n")
    }

    overlap_indices = as.data.table(findOverlaps(r_dat, target_regions, type = overlap_type))



    colnames(overlap_indices) <- c("xid", "yid")

    if(nrow(overlap_indices) == 0){
        stop("No overlaps detected")
    }


    overlap_indices[,yid := paste0("rid_", yid)]
    n_overlap_cpgs = overlap_indices[, .N, yid]
    colnames(n_overlap_cpgs) = c('rid', 'n_overlap_CpGs')

    #overlap_indices = split(overlap_indices, as.factor(as.character(overlap_indices$yid)))

    if(n_chunks==1){
        if (type == "M") {
            dat = get_matrix(m = m[overlap_indices$xid,], type = "M", add_loci = TRUE)
        } else if (type == "C") {
            dat = get_matrix(m = m[overlap_indices$xid,], type = "C", add_loci = TRUE)
        }
    } else {
        #To handle the situation where the number of overlapping sites are lower than the n_chunks
        if(nrow(overlap_indices) < n_chunks){
            n_chunks <- nrow(overlap_indices)
            if (n_cores > n_chunks){
                n_cores <- n_chunks 
                #message("n_cores should be set to be less than or equal to n_chunks.","\n","n_chunks has been set to be equal to n_cores = ", n_cores)
            }
        }
        
        if (type == "M") {
            dat = tryCatch(do.call("rbind",mclapply(mc.cores=n_cores, 
                                           split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid)/ceiling(length(overlap_indices$xid)/n_chunks))), 
                                           function(i) {
                get_matrix(m[i,], type = "M", add_loci = TRUE)
            })), 
            error = function(e){
                message( e, "\n Try to use less cores. ")
            })
        } else if (type == "C") {
            dat = tryCatch(do.call("rbind",mclapply(mc.cores=n_cores, split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid)/ceiling(length(overlap_indices$xid)/n_chunks))), function(i) {
                get_matrix(m[i,], type = "C", add_loci = TRUE)
            })), 
            error = function(e){
                message( e, "\n Try to use less cores. ")
            })
        }

    }



    if(nrow(overlap_indices) != nrow(dat)){
        stop("Something went wrong")
    }

    dat = cbind(overlap_indices, dat)

    #message("-Summarizing overlaps..\n")
    if(how == "mean") {
        message("-Summarizing by average\n")
        output = dat[, lapply(.SD, mean, na.rm = na_rm), by = yid, .SDcols = c(elementMetadata.col,rownames(colData(m)))]
    } else if (how == "median") {
        message("-Summarizing by median\n")
        output = dat[, lapply(.SD, median, na.rm = na_rm), by = yid, .SDcols = c(elementMetadata.col,rownames(colData(m)))]
    } else if (how == "max") {
        message("-Summarizing by maximum\n")
        output = dat[, lapply(.SD, max, na.rm = na_rm), by = yid, .SDcols = c(elementMetadata.col,rownames(colData(m)))]
    } else if (how == "min") {
        message("-Summarizing by minimum\n")
        output = dat[, lapply(.SD, min, na.rm = na_rm), by = yid, .SDcols = c(elementMetadata.col,rownames(colData(m)))]
    } else if (how == "sum") {
        message("-Summarizing by sum\n")
        output = dat[, lapply(.SD, sum, na.rm = na_rm), by = yid, .SDcols = c(elementMetadata.col,rownames(colData(m)))]
    }

    output = merge(target_regions, output, by.x = 'rid', by.y = 'yid', all.x = TRUE)
    output = merge(n_overlap_cpgs, output, by = 'rid')
    output$rid <- as.numeric(gsub("rid_","",output$rid))

    output <- output[order(output$rid),]
    setnames(output, "seqnames", "chr")
    keep <- c("chr", "start", "end", "n_overlap_CpGs", "rid", elementMetadata.col, colnames(m))
    output <- output[, ..keep]


    if(verbose){
        message("-Done! Finished in:",data.table::timetaken(start_proc_time),"\n")
    }

    return(output)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Order mathrix object by SD
#' @details Takes \code{\link{methrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' order_by_sd(m = methrix_data)
#' @export
order_by_sd <- function(m) {

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    if (is_h5(m)) {
        sds <- DelayedMatrixStats::rowSds(x = get_matrix(m, type = "M"))
    } else {
        sds <- matrixStats::rowSds(x = get_matrix(m, type = "M"))
    }
    row_order <- order(sds, na.last = TRUE, decreasing = TRUE)
    m <- m[row_order, ]
    m
}

#--------------------------------------------------------------------------------------------------------------------------

#' Subsets \code{\link{methrix}} object based on given conditions.
#' @details Takes \code{\link{methrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{methrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs chromosome names to subset by
#' @param samples sample names to subset by
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description,
#' see the \code{foverlaps} function of the \code{\link{data.table}} package.
#' @examples
#' data('methrix_data')
#' #Subset to chromosome 1
#' subset_methrix(methrix_data, contigs = 'chr21')
#' @return An object of class \code{\link{methrix}}
#' @export
subset_methrix <- function(m, regions = NULL, contigs = NULL, samples = NULL, overlap_type="within") {

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    r_dat <- data.table::as.data.table(rowData(m))
    if (!is.null(regions)) {
        message("-Subsetting by genomic regions")

        target_regions <- cast_ranges(regions)

        r_dat[, `:=`(end, start + 1)]
        #data.table::setDT(x = r_dat, key = c("chr", "start", "end"))
        overlaps <- data.table::foverlaps(x = r_dat, y = target_regions,
            type = overlap_type, nomatch = NULL, which = TRUE)
        if (nrow(overlaps) == 0) {
            stop("Subsetting resulted in zero entries")
        }

        m <- m[overlaps$xid, ]
    }

    if (!is.null(contigs)) {
        message("-Subsetting by contigs")
        selected_rows <- which(r_dat$chr %in% contigs)

        if (length(selected_rows) == 0) {
            stop("Subsetting resulted in zero entries")
        }
        m <- m[selected_rows, ]
    }

    if (!is.null(samples)) {
        message("Subsetting by samples")

        samples <- which(rownames(colData(m)) %in% samples)
        if (length(samples) == 0) {
            stop("None of the samples are present in the object")
        }

        m <- m[, samples]
    }

    return(m)
}

#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by coverage
#' @details Takes \code{\link{methrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{methrix}} object
#' @param cov_thr minimum coverage required to call a loci covered
#' @param min_samples Minimum number of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param prop_samples Minimum proportion of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be
#' tested group-wise.
#' @param n_chunks Number of chunks to split the \code{\link{methrix}} object in case it is very large. Default = 1.
#' @param n_cores Number of parallel instances. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @importFrom methods is as new
#' @examples
#' data('methrix_data')
#' #keep only CpGs which are covered by at-least 1 read across 3 samples
#' coverage_filter(m = methrix_data, cov_thr = 1, min_samples = 3)
#' @return An object of class \code{\link{methrix}}
#' @export
coverage_filter <- function(m, cov_thr = 1, min_samples = 1, prop_samples=0, group = NULL, n_chunks=1, n_cores=1) {

    start_proc_time <- proc.time()
    V1 <- . <- col2 <- Count2 <- i.to <- NULL
    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    if (!is.numeric(cov_thr)){
        stop("cov_thr is not numeric.")
    }

    if (!(is.numeric(min_samples) & is.numeric(prop_samples))){
        stop("min_samples and prop_samples variables are not numeric.")
    }

    if (!is.null(group) && !(group %in% colnames(m@colData))){
        stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."))
    }


    if (n_cores > n_chunks){
        n_chunks <- n_cores
        message("n_cores should be set to be less than or equal to n_chunks.", "\n", "n_chunks has been set to be equal to n_cores = ", n_cores)
    }



    if (is_h5(m)) {
        if (n_chunks == 1) {
            cov_dat = get_matrix(m = m, type = "C")
            if (!is.null(group)) {
                row_idx <- sapply(unique(m@colData[, group]), function(c) {
                    res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[, 
                                                                            group] == c] >= cov_thr, na.rm = TRUE)
                    row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
                                                                    sum(m@colData[, group] == c))))
                })
                row_idx <- DelayedMatrixStats::rowAlls(row_idx)
            } else {
                res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
                row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
                                                                ncol(cov_dat))))
            }
        } else {
            if (!is.null(group)) {
                row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
                                           function(i) {
                                               cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
                                                                    type = "C")
                                               row_idx <- sapply(unique(m@colData[, group]), function(c) {
                                                   res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[,group] == c] >= cov_thr, na.rm = TRUE)
                                                   row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
                                               })
                                               row_idx <- DelayedMatrixStats::rowAlls(row_idx)
                                           }))
            } else {
                row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
                                           function(i) {
                                               cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
                                                                    type = "C")
                                               res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
                                               row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
                                           }))
            }
        }
    } else {
        cov_dat = get_matrix(m = m, type = "C")
        if (!is.null(group)) {
            row_idx <- sapply(unique(m@colData[, group]), function(c) {
                res <- matrixStats::rowSums2(cov_dat[, m@colData[, group] == 
                                                         c] >= cov_thr, na.rm = TRUE)
                row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
            })
            row_idx <- matrixStats::rowAlls(row_idx)
        } else {
            res <- matrixStats::rowSums2(cov_dat >= cov_thr, na.rm = T)
            row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
        }
    }

    gc()
    message(paste0("-Retained ", format(sum(row_idx), big.mark = ","),
                   " of ", format(nrow(m), big.mark = ","), " sites"))
    message("-Finished in:  ", data.table::timetaken(start_proc_time))


    return(m[row_idx, ])
}


#--------------------------------------------------------------------------------------------------------------------------

#' Extract methylation or coverage matrices
#' @details Takes \code{\link{methrix}} object and returns user specified \code{methylation} or \code{coverage} matrix
#' @param m \code{\link{methrix}} object
#' @param type can be \code{M} or \code{C}. Default 'M'
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{GRanges}?
#' @return Coverage or Methylation matrix
#' @examples
#' data('methrix_data')
#' #Get methylation matrix
#' get_matrix(m = methrix_data, type = 'M')
#' #Get methylation matrix along with loci
#' get_matrix(m = methrix_data, type = 'M', add_loci = TRUE)
#' #' #Get methylation data as a GRanges object
#' get_matrix(m = methrix_data, type = 'M', add_loci = TRUE, in_granges=TRUE)
#' @export
get_matrix <- function(m, type = "M", add_loci = FALSE, in_granges=FALSE) {

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }


    type <- match.arg(arg = type, choices = c("M", "C"))
    if (add_loci == FALSE & in_granges == TRUE) {
        warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
                "the output will be a data.table object. ")
    }

    if (type == "M") {
        d <- SummarizedExperiment::assay(x = m, i = 1)
    } else {
        d <- SummarizedExperiment::assay(x = m, i = 2)
    }

    if (add_loci) {
        if (is_h5(m)) {
            d <- as.data.frame(cbind(SummarizedExperiment::rowData(x = m),
                as.data.frame(d)))
        } else {
            d <- as.data.frame(cbind(SummarizedExperiment::rowData(x = m),
                d))
        }
        if (in_granges){
            d$end <- d$start +1
            d <- GenomicRanges::makeGRangesFromDataFrame(d, keep.extra.columns = TRUE)
        } else {
        data.table::setDT(x = d)
        }
    }
    d
}

#--------------------------------------------------------------------------------------------------------------------------

#' Convert \code{\link{methrix}} to \code{bsseq} object
#' @details Takes \code{\link{methrix}} object and returns a \code{bsseq} object
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{bsseq}
#' @examples
#' \dontrun{
#' data('methrix_data')
#' methrix2bsseq(m = methrix_data)
#' }
#' @export
#'
methrix2bsseq <- function(m) {

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    n_samps <- nrow(SummarizedExperiment::colData(x = m))
    M_clean <- get_matrix(m) * get_matrix(m, type = "C")
    M_clean[is.na(M_clean)] <- 0
    assays(m)[[2]][is.na(assays(m)[[2]])] <- 0

    b <- bsseq::BSseq(M = M_clean, Cov = get_matrix(m, type = "C"), pData = colData(x = m),
        pos = rowData(x = m)[, "start"], chr = rowData(x = m)[, "chr"],
        sampleNames = rownames(m@colData))
    b
}

#--------------------------------------------------------------------------------------------------------------------------

#' Remove loci that are uncovered across all samples
#' @details Takes \code{\link{methrix}} object and removes loci that are uncovered across all samples
#' @param m \code{\link{methrix}} object
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' remove_uncovered(m = methrix_data)
#' @export
#'
remove_uncovered <- function(m) {

    start_proc_time <- proc.time()
    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }


    row_idx<-rowSums(!is.na(assays(m)[[2]]))==0

    message(paste0("-Removed ", format(sum(row_idx), big.mark = ","),
                   " [", round(sum(row_idx)/nrow(m) * 100, digits = 2), "%] uncovered loci of ",
                   format(nrow(m), big.mark = ","), " sites"))

    gc()

    message("-Finished in:  ", data.table::timetaken(start_proc_time))

    if (sum(row_idx) == 0){
        m
    } else {
        m[!row_idx, ]
    }


}

#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by region
#' @details Takes \code{\link{methrix}} object and filters CpGs based on supplied regions in data.table or GRanges format
#' @param m \code{\link{methrix}} object
#' @param regions genomic regions to filter-out. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param type defines the type of the overlap of the CpG sites with the target regions. Default value is `within`. For detailed description,
#' see the \code{foverlaps} function of the \code{\link{data.table}} package.
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' region_filter(m = methrix_data,
#' regions = data.table(chr = 'chr21', start = 27867971, end =  27868103))
#' @export
region_filter <- function(m, regions, type = "within") {

    start_proc_time <- proc.time()
    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    target_regions <- cast_ranges(regions)

    current_regions <- data.table::as.data.table(rowData(m))
    current_regions[, `:=`(end, start + 1)]
    data.table::setDT(x = current_regions, key = c("chr", "start", "end"))
    overlap <- data.table::foverlaps(x = current_regions, y = target_regions,
        type = type, nomatch = NULL, which = TRUE)

    if (nrow(overlap) == 0) {
        stop("No CpGs found within the query intervals. Nothing to remove.")
    }

    message(paste0("-Removed ", format(nrow(overlap), big.mark = ","),
        " CpGs"))
    message("-Finished in:  ", data.table::timetaken(start_proc_time))

    m[-overlap$xid, ]
}

#--------------------------------------------------------------------------------------------------------------------------

#' Masks too high or too low coverage
#' @details Takes \code{\link{methrix}} object and masks sites with too high or too low coverage
#'  by putting NA for coverage and beta value. The sites will remain in the object.
#' @param m \code{\link{methrix}} object
#' @param low_count The minimal coverage allowed. Everything below, will get masked. Default = NULL, nothing gets masked.
#' @param high_quantile The quantile limit of coverage. Quantiles are calculated for each sample and everything that belongs to a
#' higher quantile than the defined will be masked. Default = 0.99.
#' @param n_cores Number of parallel instances. Can only be used if \code{\link{methrix}} is in HDF5 format. Default = 1.
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' mask_methrix(m = methrix_data, low_count = 5, high_quantile = 0.99 )
#' @export
mask_methrix <- function(m, low_count = NULL, high_quantile = 0.99, n_cores=1) {

    start_proc_time <- proc.time()
    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    if (!is_h5(m) & n_cores != 1) {
        stop("Parallel processing not supported for a non-HDF5 methrix object due to probable high memory usage. \nNumber of cores (n_cores) needs to be 1.")
    }


    if (!is.null(low_count)) {

        if(!is.numeric(low_count)){
            stop("low_count must be a numeric value.")
        }


        message("-Masking coverage lower than ", low_count)

        row_idx1 <- (assays(m)[[2]] < low_count)
        assays(m)[[1]][row_idx1] <- as.double(NA)
        assays(m)[[2]][row_idx1] <- as.integer(NA)


        if (is_h5(m)) {
            if (n_cores == 1){
                n <- DelayedMatrixStats::colSums2(row_idx1, na.rm = T)}
            else {
                    n <- simplify2array(mclapply(mc.cores = n_cores, 1:ncol(row_idx1), function(i) sum(row_idx1[, i], na.rm = T)))}
        } else {
            n <- colSums(row_idx1, na.rm = T)
        }

        for (i in seq_along(colnames(m))) {
            message(paste0("-Masked ", n[i], " CpGs due to too low coverage in sample ",  colnames(row_idx1)[i], "."))

        }

    }


    if (!is.null(high_quantile)) {
        if (high_quantile >= 1 | high_quantile <= 0) {
            stop("High quantile should be between 0 and 1. ")
        }

        message("\n\n-Masking coverage higher than ", high_quantile*100," percentile")


        
        if (is_h5(m)) {
            
            if (n_cores == 1){
                quantiles <- DelayedMatrixStats::colQuantiles(assays(m)[[2]], probs = high_quantile, na.rm = TRUE, drop = F)}
            else {
                quantiles <- simplify2array(mclapply(mc.cores = n_cores, 1:ncol(assays(m)[[2]]), 
                                                     function(i) quantile(assays(m)[[2]][, i], probs = high_quantile, na.rm = TRUE)))}
            quantiles <- as.vector(quantiles)
            names(quantiles) <- rownames(m@colData)
        } else {
            quantiles <- matrixStats::colQuantiles(assays(m)[[2]], probs = high_quantile, na.rm = TRUE)
            quantiles <- as.vector(quantiles)
            names(quantiles) <- rownames(m@colData)
        }
        

        row_idx2 <- t(t((assays(m)[[2]])) > quantiles)
        assays(m)[[1]][row_idx2] <- as.double(NA)
        assays(m)[[2]][row_idx2] <- as.integer(NA)


        if (is_h5(m)) {
            if (n_cores == 1){
                n <- DelayedMatrixStats::colSums2(row_idx2, na.rm = T)}
            else {
                n <- simplify2array(mclapply(mc.cores = n_cores, 1:ncol(row_idx2), 
                                             function(i) sum(row_idx2[, i], na.rm = T)))}
        } else {
            n <- colSums(row_idx2, na.rm = T)
        }


        for (i in seq_along(colnames(m))) {
            message(paste0("-Masked ", n[i], " CpGs due to too high coverage in sample ",
                           colnames(row_idx2)[i], "."))

        }


    }
    message("-Finished in:  ", data.table::timetaken(start_proc_time))
    return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Combine methrix objects
#' @details Takes two \code{\link{methrix}} objects and combines them row- or column-wise
#' @param m1 Frist \code{\link{methrix}} object
#' @param m2 Second \code{\link{methrix}} object
#' @param by The direction of combine. 'column' (cbind) combines samples with same regions, 'row' combines different regions,
#' e.g. different chromosomes.
#' @return An object of class \code{\link{methrix}}
#' @export
#'
combine_methrix <- function(m1, m2, by = c("row", "col")) {

    if (!is(m1, "methrix") || !is(m2, "methrix") ){
        stop("Valid methrix objects need to be supplied.")
    }

    by <- match.arg(arg = by, choices = c("row", "col"), several.ok = FALSE)

    if (by == "row") {
        if (nrow(colData(m1))!= nrow(colData(m2))  || !(all(rownames(m1@colData) == rownames(m2@colData)))) {
            stop("You have different samples in your dataset. You need the same samples in your datasets. ")
        } else {
            m <- rbind(m1, m2)
        }
        if (any(duplicated((as.data.table(m@elementMetadata))))) {
            stop("There are overlapping regions in your datasets. This function only takes distinct objects. ")
        }
    }
    if (by == "col") {
        if (any(rownames(m1@colData) %in% rownames(m2@colData))) {
            stop("You have the same samples in your datasets. You need different samples for this merging.  ")
        } else if (!identical(m1@elementMetadata, m2@elementMetadata)) {
            stop("You have to have the same regions in your datasets. ")
        } else {
            m <- cbind(m1, m2)
        }
    }
    gc()
    return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Estimate descriptive statistics
#' @details Calculate descriptive statistics
#' @param m \code{\link{methrix}} object
#' @param per_chr Estimate stats per chromosome. Default TRUE
#' @seealso \code{\link{plot_stats}}
#' @examples
#' data('methrix_data')
#' get_stats(methrix_data)
#' @return data.table of summary stats
#' @export
get_stats <- function(m, per_chr = TRUE) {

    median <- . <- sd <- chr <- NULL
    start_proc_time <- proc.time()
    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }


    #Counting NA's seems to be not-ideal (especially in cases of masked/coverage filtered cases.)
    # row_idx <- data.table::as.data.table(which(is.na(get_matrix(m = m,
    #     type = "C")), arr.ind = TRUE))
    # colnames(row_idx) <- c("row", "col")
    # row_idx <- split(row_idx, as.factor(as.character(row_idx$col)))

    if (per_chr) {
        cov_stat <- lapply(1:ncol(m), function(i){
            get_matrix(m = m[,i], type = "C", add_loci = TRUE)[, c(1, 4), with = FALSE][, .(
                mean_cov = lapply(.SD,
                                  matrixStats::mean2, na.rm = TRUE),
                median_cov = lapply(.SD,
                                    median, na.rm = TRUE),
                sd_cov = lapply(.SD, sd, na.rm = TRUE)
            ),
            by = chr]
        })

        meth_stat <- lapply(1:ncol(m), function(i){
            get_matrix(m = m[,i], type = "M", add_loci = TRUE)[, c(1, 4), with = FALSE][, .(
                mean_meth = lapply(.SD,
                                  matrixStats::mean2, na.rm = TRUE),
                median_meth = lapply(.SD,
                                    median, na.rm = TRUE),
                sd_meth = lapply(.SD, sd, na.rm = TRUE)
            ),
            by = chr]
        })

        names(meth_stat) <-  colnames(m)
        names(cov_stat) <- colnames(m)

        cov_stat <- data.table::rbindlist(l = cov_stat, use.names = TRUE,
                                          idcol = "Sample_Name")
        meth_stat <- data.table::rbindlist(l = meth_stat, use.names = TRUE,
                                           idcol = "Sample_Name")
        stats <- merge(meth_stat, cov_stat, by = c("chr", "Sample_Name"))
        colnames(stats)[1] <- "Chromosome"
        stats$Chromosome <- factor(x = stats$Chromosome, levels = m@metadata$chrom_sizes$contig)
    } else {
        if (is_h5(m)) {
            cov_stat <- data.table::data.table(Sample_Name = colnames(m),
                mean_cov = DelayedMatrixStats::colMeans2(get_matrix(m = m, "C"), na.rm = TRUE),
                median_cov = DelayedMatrixStats::colMedians(get_matrix(m = m, "C"), na.rm = TRUE),
                sd_cov = DelayedMatrixStats::colSds(get_matrix(m = m, "C"), na.rm = TRUE)
            )

            meth_stat <- data.table::data.table(Sample_Name = colnames(m),
                mean_meth = DelayedMatrixStats::colMeans2(get_matrix(m = m, "M"), na.rm = TRUE),
                median_meth = DelayedMatrixStats::colMedians(get_matrix(m = m, "M"), na.rm = TRUE),
                sd_meth = DelayedMatrixStats::colSds(get_matrix(m = m, "M"), na.rm = TRUE)
            )

        } else {
            cov_stat <- data.table::data.table(Sample_Name = colnames(m),
                mean_cov = matrixStats::colMeans2(get_matrix(m = m, "C"), na.rm = TRUE),
                median_cov = matrixStats::colMedians(get_matrix(m = m, "C"), na.rm = TRUE),
                sd_cov = matrixStats::colSds(get_matrix(m = m, "C"), na.rm = TRUE)
            )

            meth_stat <- data.table::data.table(Sample_Name = colnames(m),
                mean_meth = matrixStats::colMeans2(get_matrix(m = m, "M"), na.rm = TRUE),
                median_meth = matrixStats::colMedians(get_matrix(m = m, "M"), na.rm = TRUE),
                sd_meth = matrixStats::colSds(get_matrix(m = m, "M"), na.rm = TRUE)
            )
        }


        stats <- merge(meth_stat, cov_stat, by = c("Sample_Name"))
    }

    gc()
    message("-Finished in:  ", data.table::timetaken(start_proc_time))

    stats
}

#--------------------------------------------------------------------------------------------------------------------------
#' Saves HDF5 methrix object
#' @details Takes \code{\link{methrix}} object and saves it
#' @param m \code{\link{methrix}} object
#' @param dir The directory to use. Created, if not existing. Default NULL
#' @param replace Should it overwrite the pre-existing data? FALSE by default.
#' @param ... Parameters to pass to saveHDF5SummarizedExperiment
#' @examples
#' data('methrix_data')
#' methrix_data_h5 <- convert_methrix(m=methrix_data)
#' target_dir = paste0(getwd(), '/temp/')
#' save_HDF5_methrix(methrix_data_h5, dir = target_dir, replace = TRUE)
#' @return Nothing
#' @export
save_HDF5_methrix <- function(m = NULL, dir = NULL, replace = FALSE, ...) {


    if (is.null(dir)) {
        stop("Please provide an target directory to save the results")
    }

    if (is(m, "methrix") && is_h5(m)) {
        HDF5Array::saveHDF5SummarizedExperiment(x = m, dir = dir, replace = replace,
            ...)
    } else {
        stop("The object is not a methrix object or not in an HDF5 format. ")
    }

}

#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 methrix object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{methrix}} object and loads it
#' @param dir The directory to read in from. Default NULL
#' @param ... Parameters to pass to loadHDF5SummarizedExperiment
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' methrix_data_h5 <- convert_methrix(m=methrix_data)
#' target_dir = paste0(getwd(), '/temp1/')
#' save_HDF5_methrix(methrix_data_h5, dir = target_dir, replace = TRUE)
#' load_HDF5_methrix(target_dir)
#' @export
load_HDF5_methrix <- function(dir = NULL, ...) {

    if (is.null(dir)) {
        stop("Please provide the target directory containing ")
    }

    m <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
    m <- as(m, "methrix")
    return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts HDF5 methrix object to standard in-memory object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with in-memory assay slots.
#' @param m An object of class \code{\link{methrix}}, HDF5 format
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' m <- convert_HDF5_methrix(m=m2)
#' @export
convert_HDF5_methrix <- function(m = NULL) {

    if (is.null(m) | !is(m, "methrix")) {
        stop("No or not valid input data provided.")
    }
    if (!is_h5(m)) {
        stop("The input data is not in HDF5 format. No conversion happened.")
    }

    assays(m)[[1]] <- as.matrix(assays(m)[[1]])
    assays(m)[[2]] <- as.matrix(assays(m)[[2]])
    m@metadata$is_h5 <- FALSE
    return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts an in-memory object to an on-disk HDF5 object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @param m An object of class \code{\link{methrix}}
#' @return An object of class \code{\link{methrix}}, HDF5 format
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' @export
convert_methrix <- function(m = NULL) {

    if (is.null(m) | !is(m, "methrix")) {
        stop("No or not valid input data provided.")
    }
    if (is_h5(m)) {
        stop("The input data is already in HDF5 format. No conversion happened.")
    }

    m <- create_methrix(beta_mat = assays(m)[[1]], cov_mat = assays(m)[[2]],
        cpg_loci = m@elementMetadata, is_hdf5 = TRUE, genome_name = m@metadata$genome,
        col_data = m@colData, chrom_sizes = m@metadata$chrom_sizes, ref_cpg_dt = m@metadata$ref_CpG,
        desc = m@metadata$descriptive_stats)
    return(m)
}
