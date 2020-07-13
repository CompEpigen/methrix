#' Removes CpG sites from the object if they overlap with common SNPs
#' @details Takes \code{\link{methrix}} object and removes common SNPs. SNPs overlapping with a CpG site and have a minor allele frequency (MAF)
#' above a threshold in any of the populations used will be selected and the corresponding CpG sites will be removed from the \code{\link{methrix}} object.
#' With the reduce_filtering option, SNPs with MAP < 0.1 will be further evaluated. If they show low variance in the dataset, there is probably no genotype
#' variability in the population, therefore the corresponding CpG site won't be removed. Please keep in mind that variance thresholds are
#' @param m \code{\link{methrix}} object
#' @param populations Populations to use. Default is all.
#' @param maf_threshold The frequency threshold, above which the SNPs will be removed. Default is 0.01
#' @param reduce_filtering If TRUE, the SNPs with a MAF < 0.1 will be evaluated and only the highly variable ones will be removed.
#' Default FALSE.
#' @param forced the reduce_filtering is not recommended with less than 10 samples, but can be forced. Default is FALSE.
#' @param keep Do you want to keep the sites that were filtered out? In this case, the function will return with a list of wo methrix objects.
#' @param n_chunks Number of chunks to split the \code{\link{methrix}} object in case it is very large. Can only be used if input data is in HDF5 format. Default = 1.
#' @param n_cores Number of parallel instances. Can only be used if input data is in HDF5 format. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @return methrix object or a list of methrix objects
#' @importFrom BSgenome score
#' @examples
#' data('methrix_data')
#' remove_snps(m = methrix_data, maf_threshold=0.01)
#' @export
remove_snps <- function(m, populations = NULL, maf_threshold = 0.01, reduce_filtering = FALSE,    forced = FALSE, keep = FALSE, n_chunks=1, n_cores=1) {

    start_proc_time <- proc.time()

    genome <- m@metadata$genome
    chr <- m2 <- NULL



    if (!is_h5(m) & n_cores!=1){
        stop("Parallel processing not supported for a non-HDF5 methrix object due to probable high memory usage. \nNumber of cores (n_cores) needs to be 1.")
    }

    if (!is_h5(m) & n_chunks!=1){
        stop("Splitting into chunks not supported for a non-HDF5 methrix object. \nNumber of chunks (n_chunks) needs to be 1.")
    }

    if (n_cores>n_chunks){
        n_chunks<-n_cores
        message("n_cores should be set to be less than or equal to n_chunks.","\n","n_chunks has been set to be equal to n_cores = ",n_cores)
    }


    if (grepl("hg19|GRCh37|Hs37|hs37", genome)) {
        if (requireNamespace("MafDb.1Kgenomes.phase3.hs37d5", quietly = TRUE) &
            requireNamespace("GenomicScores", quietly = TRUE)) {
            message("Used SNP database: MafDb.1Kgenomes.phase3.hs37d5. \n")
            mafdb <- MafDb.1Kgenomes.phase3.hs37d5::MafDb.1Kgenomes.phase3.hs37d5
        } else {
            stop("Packages MafDb.1Kgenomes.phase3.hs37d5 and/or GenomicScores not found.
                 Please install them before proceed.")
        }
    } else if (grepl("hg38|GRCh38|Hs38|hs38", genome)) {
        if (requireNamespace("MafDb.1Kgenomes.phase3.GRCh38", quietly = TRUE) &
            requireNamespace("GenomicScores", quietly = TRUE)) {
            message("Used SNP database: MafDb.1Kgenomes.phase3.hs38. \n")
            mafdb <- MafDb.1Kgenomes.phase3.GRCh38::MafDb.1Kgenomes.phase3.GRCh38
        } else {
            stop("Packages MafDb.1Kgenomes.phase3.GRCh38 and/or GenomicScores not found.
                 Please install them before proceed.")
        }
    } else {
        stop("Only hg19 and hg38 genomes are currently supported..\n")
    }
    if (is.null(populations)) {
        populations <- GenomicScores::populations(mafdb)
    } else if (!all(populations %in% GenomicScores::populations(mafdb))) {
        stop("The population selected is not available. Please select from ",
             paste(GenomicScores::populations(mafdb), collapse = ", "), ". \n")
    }


    regions <- gr.nochr(GenomicRanges::makeGRangesFromDataFrame(elementMetadata(m), start.field = "start", end.field = "start"))


    if(n_cores==1) {
        snp_rows<-which((rowSums(as.matrix(score(mafdb, regions, pop = populations))>= maf_threshold,na.rm=T)>0) |
                            (rowSums(as.matrix(score(mafdb,  IRanges::shift(regions, 1), pop = populations))>= maf_threshold,na.rm=T)>0))
    } else {
        snp_rows<-which(unlist(mclapply(mc.cores=n_cores,1:n_chunks, function(i){
            sub_dat = regions[((i-1)*ceiling(length(regions)/n_chunks)+1):min(i*ceiling(length(regions)/n_chunks),length(regions))]
            (rowSums(as.matrix(score(mafdb, sub_dat, pop = populations))>= maf_threshold,na.rm=T)>0) |
                (rowSums(as.matrix(score(mafdb,  IRanges::shift(sub_dat, 1), pop = populations))>= maf_threshold,na.rm=T)>0)
        })))
    }

    snp_rows <- snp_rows[order(snp_rows)]



    if (reduce_filtering) {

        message("Keep in mind that the filtering is in experimental state, the cut-off is arbitrary.")

        if (ncol(m) < 10 & forced == FALSE) {
            stop("The reduce filtering option is not recommended for sample number below 10.
                 Use the forced option if you still want to do it.")
        } else if (ncol(m) < 10 & forced == TRUE) {
            message("The reduce filtering option is not recommended for sample number below 10.")
        }




        if(n_chunks==1) {
            snp_test<-which((rowSums(as.matrix(score(mafdb, regions, pop = populations))<=0.1,na.rm=T)>0) |
                                (rowSums(as.matrix(score(mafdb,  IRanges::shift(regions, 1), pop = populations))<=0.1,na.rm=T)>0))
        } else {
            snp_test<-which(unlist(mclapply(mc.cores=n_cores,1:n_chunks, function(i){
                sub_dat = regions[((i-1)*ceiling(length(regions)/n_chunks)+1):min(i*ceiling(length(regions)/n_chunks),length(regions))]
                (rowSums(as.matrix(score(mafdb, sub_dat, pop = populations))<=0.1,na.rm=T)>0) |
                    (rowSums(as.matrix(score(mafdb,  IRanges::shift(sub_dat, 1), pop = populations))<=0.1,na.rm=T)>0)
            })))
        }




        snp_test <- snp_test[order(snp_test)]
        snp_test <- intersect(snp_rows, snp_test)
        snp_rows <- snp_rows[!(snp_rows %in% snp_test)]
        m_mat <- get_matrix(m[snp_test, ])
        if (is_h5(m)) {

            if(n_chunks==1)  {
                vars <- DelayedMatrixStats::rowVars(m_mat, na.rm = TRUE)
            } else {
                vars<-unlist(mclapply(mc.cores=n_cores,1:n_chunks, function(i){
                    m_mat_dat = m_mat[((i-1)*ceiling(nrow(m_mat)/n_chunks)+1):min(i*ceiling(nrow(m_mat)/n_chunks),nrow(m_mat)),]
                    vars<-DelayedMatrixStats::rowVars(m_mat_dat, na.rm = TRUE)
                }))
            }



        } else {
            vars <- matrixStats::rowVars(m_mat, na.rm = TRUE)
        }
        snp_test <- snp_test[which(vars > quantile(vars[complete.cases(vars)],0.25))]


        snp_rows <- c(snp_rows, snp_test)
    }

    removed_snps <- data.table::as.data.table(m@elementMetadata)[snp_rows,    ]

    message("Number of SNPs removed:")
    print(removed_snps[, .N, by = chr])
    message("Sum: ")
    print(removed_snps[, .N])

    if (keep){
        m2 <- m[snp_rows,]
    }

    m <- m[-snp_rows, ]


    gc()

    message("-Finished in:  ", data.table::timetaken(start_proc_time))


    if (keep){
        return(list("snp_filtered" = m, "removed_snps" = m2))
    } else {
        return(m)
    }

}

gr.nochr <- function(gr) {

    if (grepl("^chr", GenomeInfoDb::seqlevels(gr)[1])) {
        GenomeInfoDb::seqlevels(gr) <- gsub("^chr", "", GenomeInfoDb::seqlevels(gr))
    }
    return(gr)
}
