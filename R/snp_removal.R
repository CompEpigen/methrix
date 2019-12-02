#' Removes CpG sites from the object if they overlap with common SNPs
#' @details Takes \code{\link{methrix}} object and removes common SNPs. SNPs overlapping with a CpG site and have a minor allele frequency (MAF)
#' above a threshold in any of the populations used will be selected and the corresponding CpG sites will be removed from the \code{\link{methrix}} object.
#' With the reduce_filtering option, SNPs with MAP < 0.1 will be further evaluated. If they show low variance in the dataset, there is probably no genotype
#' variability in the population, therefore the corresponding CpG site won't be removed. Please keep in mind that variance thresholds are
#' @param m \code{\link{methrix}} object
#' @param populations Populations to use. Default is all.
#' @param maf_threshold The frequency threshold, above which the SNPs will be removed. Accepted thresholds: 0.01 and 0.05. Default is 0.05
#' @param reduce_filtering If TRUE, the SNPs with a MAF < 0.1 will be evaluated and only the highly variable ones will be removed.
#' Default FALSE.
#' @param forced the reduce_filtering is not recommended with less than 10 samples, but can be forced. Default is FALSE.
#' @param keep Do you want to keep the sites that were filtered out? In this case, the function will return with a list of wo methrix objects. 
#' @return methrix object or a list of methrix objects
#' @importFrom BSgenome score
#' @examples
#' data('methrix_data')
#' remove_snps(m = methrix_data, maf_threshold=0.01)
#' @export


remove_snps <- function(m, populations = NULL, maf_threshold = 0.01, reduce_filtering = FALSE,
    forced = FALSE, keep = FALSE) {
    genome <- m@metadata$genome
    chr <- m2 <- NULL

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
    } else if (!all(populations %in% populations(mafdb))) {
        stop("The population selected is not available. Please select from ",
            paste(populations(mafdb), collapse = ", "), ". \n")
    }

    if (!(maf_threshold %in% c(0.01, 0.05))) {
        stop("The maf_threshold should be 0.01 or 0.05. \n")
    }

    regions <- gr.nochr(GenomicRanges::makeGRangesFromDataFrame(elementMetadata(m),
        start.field = "start", end.field = "start"))

    snp_rows <- unique(c(unique(which(as.data.table(score(mafdb, regions,
        pop = populations)) >= maf_threshold, arr.ind = TRUE)[, 1]),
        unique(which(as.data.table(score(mafdb,
        shift(regions, 1), pop = populations)) >= maf_threshold, arr.ind = TRUE)[, 1])))

    snp_rows <- snp_rows[order(snp_rows)]

    if (reduce_filtering) {

        message("Keep in mind that the filtering is in experimental state, the cut-off is arbitrary.")

        if (ncol(m) < 10 & forced == FALSE) {
            stop("The reduce filtering option is not recommended for sample number below 10.
                 Use the forced option if you still want to do it.")
        } else if (ncol(m) < 10 & forced == TRUE) {
            message("The reduce filtering option is not recommended for sample number below 10.")
        }

        snp_test <- unique(c(unique(which(as.data.table(score(mafdb, regions,
            pop = populations)) <= 0.1, arr.ind = TRUE)[, 1]),
            unique(which(as.data.table(score(mafdb,
            shift(regions, 1), pop = populations)) <= 0.1, arr.ind = TRUE)[, 1])))
        snp_test <- snp_test[order(snp_test)]
        snp_test <- intersect(snp_rows, snp_test)
        snp_rows <- snp_rows[!(snp_rows %in% snp_test)]
        m_mat <- get_matrix(m[snp_test, ])
        if (is_h5(m)) {
            vars <- DelayedMatrixStats::rowVars(m_mat, na.rm = TRUE)
        } else {
            vars <- matrixStats::rowVars(m_mat, na.rm = TRUE)
        }
        snp_test <- snp_test[which(vars > quantile(vars[complete.cases(vars)],
            0.25))]


        snp_rows <- c(snp_rows, snp_test)
    }

    removed_snps <- data.table::as.data.table(m@elementMetadata)[snp_rows,
        ]

    message("Number of SNPs removed:")
    print(removed_snps[, .N, by = chr])
    message("Sum: ")
    print(removed_snps[, .N])
    
    if (keep){
        m2 <- m[snp_rows,]
    }
    
    m <- m[-snp_rows, ]


    gc()
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
