#' Creates a detailed interative html summary report from Methrix object
#' @description Creates a detailed interative html summary report from Methrix object.
#' If the directory contains required files (from previous run), it directly proceeds to generate html report.
#' @param meth \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved. If \code{NULL} creats a \code{tempdir}
#' @param recal_stats Whether summary statistics should be recalculated? If you are using subsetted methrix object set this to TRUE.
#' @param plot_beta_dist Default TRUE. Can be time consuming.
#' @param beta_nCpG Number of CpGs rto use for estimating beta value distribution. Default 10000
#' @param prefix If provided, the name of the report and the intermediate files will start with the prefix. 
#' @param n_thr Default 4. Only used if \code{plot_beta_dist} is TRUE
#' @return an interactive html report
#' @examples
#' \dontrun{
#' data('methrix_data')
#' methrix::methrix_report(meth = methrix_data)
#' }
#' @export
methrix_report <- function(meth, output_dir = NULL, recal_stats = FALSE,
    plot_beta_dist = TRUE, beta_nCpG = 10000, prefix=NULL, n_thr = 4) {

    n_covered <- total_CpGs <- n_non_covered <- fract_CpG <- n_CpG <- NULL
    if (!recal_stats) {
        warning("If input methrix is a subsetted version of original methrix object, set recal_stats to TRUE",
            immediate. = TRUE)
        if (is.null(meth@metadata$descriptive_stats)) {
            stop("No previous statistics is available. Set recal_stats to TRUE.")
        }
    }
    if (!is.null(prefix)){
        if (!is.character(prefix))
            stop("The provided prefix is not valid, please provide a string.")
    }
    
    start_proc_time <- proc.time()

    if (is.null(output_dir)) {
        output_dir <- tempdir()
    } else {
        if (!dir.exists(output_dir)) {
            dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
        }
    }

    # Methylation/Coverage per chromosome
    message(paste0("Step 1 of 5: Methylation/Coverage per chromosome"))
    if (!is.null(prefix)){
        of1 <- suppressWarnings(normalizePath(file.path(output_dir, paste0(prefix, "_MC_per_chr.tsv"))))
    } else {
        of1 <- suppressWarnings(normalizePath(file.path(output_dir, "MC_per_chr.tsv")))
    }
    
    if (file.exists(of1)) {
        message("File already present. Skipping step 1..")
    } else {
        if (recal_stats) {
            per_chr_stat <- get_stats(m = meth, per_chr = TRUE)
            gc()
        } else {
            per_chr_stat <- meth@metadata$descriptive_stats$chr_stat
            colnames(per_chr_stat)[which(colnames(per_chr_stat) == "chr")] <- "Chromosome"
        }
        data.table::fwrite(x = per_chr_stat, file = of1, sep = "\t")
    }

    # Global methylation/Coverage
    message(paste0("Step 2 of 5: Global methylation/Coverage per sample\n"))
    if (!is.null(prefix)){
        of2 <- suppressWarnings(normalizePath(file.path(output_dir, paste0(prefix, "_global_MC_per_samp.tsv"))))
    } else {
        of2 <- suppressWarnings(normalizePath(file.path(output_dir, "global_MC_per_samp.tsv")))
    }
    if (file.exists(of2)) {
        message("File already present. Skipping step 2..")
    } else {
        if (recal_stats) {
            genome_stat <- get_stats(m = meth, per_chr = FALSE)
            gc()
        } else {
            genome_stat <- meth@metadata$descriptive_stats$genome_stat
        }
        data.table::fwrite(x = genome_stat, file = of2, sep = "\t")
    }

    # n CpGs covered per chromomse
    message(paste0("Step 3 of 5: Reference CpGs covered per chromosome"))
    if (!is.null(prefix)){
        of3 <- suppressWarnings(normalizePath(file.path(output_dir, paste0(prefix, "_n_covered_per_chr.tsv"))))
    } else {
        of3 <- suppressWarnings(normalizePath(file.path(output_dir, "n_covered_per_chr.tsv")))
    }
    contig_nCpGs <- meth@metadata$ref_CpG
    colnames(contig_nCpGs) <- c("chr", "total_CpGs")
    if (file.exists(of3)) {
        message("File already present. Skipping step 3..")
    } else {
        if (recal_stats) {
            non_cov_tbl <- lapply(seq_len(ncol(meth)), function(i) {
                data.table::as.data.table(as.data.frame(table(rowData(meth)[which(is.na(get_matrix(m = meth,
                  "M")[, i])), "chr"])))
            })
            names(non_cov_tbl) <- colnames(meth)
            gc()
            non_cov_tbl <- data.table::rbindlist(l = non_cov_tbl, use.names = TRUE,
                fill = TRUE, idcol = "Sample_Name")
            colnames(non_cov_tbl) <- c("Sample_Name", "chr", "n_non_covered")
            non_cov_tbl <- merge(non_cov_tbl, contig_nCpGs, by = "chr",
                all.x = TRUE)
            non_cov_tbl[, `:=`(n_covered, total_CpGs - n_non_covered)]
            non_cov_tbl[, `:=`(n_non_covered, NULL)]
        } else {
            non_cov_tbl <- data.table::melt(meth@metadata$descriptive_stats$n_cpgs_covered,
                id.vars = "chr")
            colnames(non_cov_tbl) <- c("chr", "Sample_Name", "n_covered")
            non_cov_tbl <- merge(non_cov_tbl, contig_nCpGs, by = "chr",
                all.x = TRUE)
        }
        data.table::fwrite(x = non_cov_tbl, file = of3, sep = "\t")
    }


    # Common CpGs covered by all samples
    message(paste0("Step 4 of 5: Common reference CpGs covered across all samples"))
    if (!is.null(prefix)){
        of4 <- suppressWarnings(normalizePath(file.path(output_dir, paste0(prefix, "_n_covered_by_all_samples.tsv"))))
    } else {
        of4 <- suppressWarnings(normalizePath(file.path(output_dir, "n_covered_by_all_samples.tsv")))
    }
    if (file.exists(of4)) {
        message("File already present. Skipping step 4..")
    } else {
        # na_vec = apply(get_matrix(meth, type = 'M'), MARGIN = 1, anyNA)
        if (is_h5(meth)) {
            na_vec <- DelayedMatrixStats::rowAnyNAs(get_matrix(meth, type = "M"))
        } else {
            na_vec <- matrixStats::rowAnyNAs(get_matrix(meth, type = "M"))
        }

        mf_chr_summary <- as.data.frame(table(rowData(x = meth)[, "chr"],
            na_vec))

        if (!plot_beta_dist) {
            rm(na_vec)
            gc()
        }

        data.table::setDT(x = mf_chr_summary)
        mf_chr_summary <- mf_chr_summary[na_vec == FALSE]
        mf_chr_summary[, `:=`(na_vec, NULL)]
        colnames(mf_chr_summary) <- c("chr", "n_CpG")
        mf_chr_summary <- merge(meth@metadata$ref_CpG, mf_chr_summary)
        colnames(mf_chr_summary)[2] <- c("total_CpGs")
        mf_chr_summary[, `:=`(fract_CpG, n_CpG/total_CpGs)]
        data.table::fwrite(x = mf_chr_summary, file = of4, sep = "\t")
    }

    # Density plot data
    message(paste0("Step 5 of 5: beta value distribution"))
    if (plot_beta_dist) {
        dens_files <- list.files(path = output_dir, pattern = "*_density\\.tsv\\.gz$")
        if (length(dens_files) == nrow(colData(meth))) {
            message("Files already present. Skipping step 5..")
        } else {
            if (!exists(x = "na_vec")) {
                na_vec <- apply(get_matrix(meth, type = "M"), MARGIN = 1,
                  anyNA)
            }

            row_idx <- sample(which(na_vec == FALSE), size = min(beta_nCpG,
                length(which(na_vec == FALSE))), replace = FALSE)

            lapply(X = seq_len(nrow(colData(meth))), FUN = function(i) {
                
                i_dens <- density(get_matrix(meth, type = "M")[row_idx,
                                                               i], na.rm = TRUE)
                if (!is.null(prefix)){
                    if (!dir.exists(paste0(output_dir, "/", prefix, "/"))) {
                        dir.create(path = paste0(output_dir, "/", prefix, "/"), showWarnings = FALSE, recursive = TRUE)
                    }
                    data.table::fwrite(x = data.table::data.table(x = i_dens$x,
                                                                  y = i_dens$y), file = paste0(output_dir, "/", prefix, "/", rownames(colData(x = meth))[i],
                                                                                               "_density.tsv.gz"), sep = "\t")
                } else {
                    data.table::fwrite(x = data.table::data.table(x = i_dens$x,
                                                                  y = i_dens$y), file = paste0(output_dir, "/", rownames(colData(x = meth))[i],
                                                                                               "_density.tsv.gz"), sep = "\t")
                }
            })
            rm(na_vec)
        }
        gc()
    }

    if (!is.null(prefix)){
        of5 <- suppressWarnings(normalizePath(file.path(output_dir, paste0(prefix, "_contig_lens.tsv"))))
    } else {
        of5 <- suppressWarnings(normalizePath(file.path(output_dir, "contig_lens.tsv")))
    }
    
    data.table::fwrite(x = meth@metadata$chrom_sizes, file = of5, sep = "\t")

    message(paste0("Knitting report"))
    md <- system.file("report", "summarize_methrix.Rmd", package = "methrix")

    if (!is.null(prefix)){
        output_file <- paste0(prefix, "_methrix_reports.html")
    } else {
        output_file <- "methrix_reports.html"
    }
    
    rmarkdown::render(input = md, output_file = "methrix_reports.html",
        output_dir = output_dir, clean = TRUE, params = list(prefix = prefix, n_covered_tsv = of3,
            n_covered_by_all_samples_tsv = of4, mc_per_chr_stat = of1,
            mc_per_sample_stat = of2, chr_lens = of5))
    
    browseURL(url = paste0(output_dir, "/", output_file))

    message(data.table::timetaken(started.at = start_proc_time))
}


