prepare_plot_data <- function(m, ranges = NULL, n_cpgs = 25000, pheno = NULL){

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }
    if (!is.null(n_cpgs)){
        if (!is.numeric(n_cpgs)){
            stop("n_cpgs must be numeric.")
        }
    }

    if (!is.null(ranges)) {
        meth_sub <- subset_methrix(m = m, regions = ranges)
        if (!is.null(n_cpgs)) {
            message("Randomly selecting ", n_cpgs, " sites")
            ids <- sample(x = seq_along(meth_sub), replace = FALSE, size = min(n_cpgs,
                                                                               nrow(meth_sub)))
            meth_sub <- get_matrix(m = meth_sub[ids, ], type = "M", add_loci = FALSE)
        } else {
            meth_sub <- get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
        }
    } else if (!is.null(n_cpgs)) {
        message("Randomly selecting ", n_cpgs, " sites")

        ids <- sample(x = seq_along(m), replace = FALSE, size = min(n_cpgs,
                                                                    nrow(m)))
        meth_sub <- get_matrix(m = m[ids, ], type = "M", add_loci = FALSE)
    } else {
        meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
    }


    if (!is.null(pheno)) {
        if (pheno %in% colnames(colData(m))) {
            colnames(meth_sub) <- as.character(m@colData[, pheno])
        } else {
            stop("Please provide a valid phenotype annotation column.")
        }
    }

    meth_sub <- as.data.frame(meth_sub)
    data.table::setDT(x = meth_sub)
    plot.data <- suppressWarnings(data.table::melt(meth_sub))
    colnames(plot.data) <- c("variable", "Meth")

    gc(verbose = FALSE)
    return(plot.data)


}

get_palette <- function(n_row, col_palette){

    if (!col_palette %in% rownames(RColorBrewer::brewer.pal.info)){

        stop("Please provide a valid RColorBrewer palettte. Possible values are: ", paste0(rownames(RColorBrewer::brewer.pal.info)), sep=", ")
    }
    color_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[col_palette,
                                                                                         "maxcolors"], col_palette))(n_row)
    return(color_pal)
}

#' Violin Plot for \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param n_cpgs Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param col_palette Name of the RColorBrewer palette to use for plotting.
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' data('methrix_data')
#' plot_violin(m = methrix_data)
plot_violin <- function(m, ranges = NULL, n_cpgs = 25000, pheno = NULL,
    col_palette = "RdYlGn") {
    variable <- Meth <- NULL

    plot.data <- prepare_plot_data(m=m, ranges = ranges, n_cpgs = n_cpgs, pheno = pheno)

    col_palette <- get_palette(ncol(m), col_palette)
    # generate the violin plot
    p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = variable, y = Meth,
        fill = variable)) + ggplot2::geom_violin(alpha = 0.8) + ggplot2::theme_classic(base_size = 14) +
        ggplot2::scale_fill_manual(values = col_palette) +
        ggplot2::xlab(pheno) + ggplot2::ylab(expression(beta * "-Value")) +
        theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
            colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_blank(), legend.title = element_blank())

    p
}

#--------------------------------------------------------------------------------------------------------------------------
#' Density Plot of \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param n_cpgs Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param col_palette Name of the RColorBrewer palette to use for plotting.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data('methrix_data')
#' plot_density(m = methrix_data)
plot_density <- function(m, ranges = NULL, n_cpgs = 25000, pheno = NULL,
     col_palette = "RdYlGn") {

    variable <- Meth <- NULL

    plot.data <- prepare_plot_data(m=m, ranges = ranges, n_cpgs = n_cpgs, pheno = pheno)
    col_palette <- get_palette(ncol(m), col_palette)

    # generate the density plot
    p <- ggplot2::ggplot(plot.data, ggplot2::aes(Meth, color = variable)) +
        geom_density(lwd = 1, position = "identity") + ggplot2::theme_classic() +
        ggplot2::xlab("Methylation") + ggplot2::theme_classic(base_size = 14) +
        ggplot2::scale_fill_manual(values = col_palette) +
        ggplot2::xlab(expression(beta * "-Value")) + theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12,
            colour = "black"), axis.title.y = element_blank(), legend.title = element_blank())

    gc(verbose = FALSE)

    p
}

#--------------------------------------------------------------------------------------------------------------------------
#' Principal Component Analysis
#'
#' @param m Input \code{\link{methrix}} object
#' @param top_var Number of variable CpGs to use. Default 1000 Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno Column name of colData(m). Default NULL. Will be used as a factor to color different groups
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param do_plot Should a plot be generated?
#' @param n_pc Default 2.
#' @return PCA results
#' @examples
#' data('methrix_data')
#' methrix_pca(methrix_data, do_plot = FALSE)
#' @export
#'
methrix_pca <- function(m, var = "top", top_var = 1000, ranges = NULL,
    pheno = NULL, do_plot = TRUE, n_pc = 2) {
    var_select <- match.arg(var, c("top", "rand"))

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }

    if (!(is.numeric(top_var) & is.numeric(n_pc))){
        stop("Either top_var or n_pc variables are not numeric.")
    }

    ## subset based on the input ranges
    if (!is.null(ranges)) {
        message("GenomicRanges will be used for the PCA")
        meth_sub <- subset_methrix(m = m, regions = ranges)
        meth_sub <- methrix::get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
    }

    if (is.null(top_var)) {
        message("All CpGs in the dataset will be used for the PCA")
        if (is.null(ranges)) {
            meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
        }
    } else {
        if (!is.numeric(top_var)){
            stop("top_var must be numeric.")
        }
        top_var <- as.integer(as.character(top_var))
        if (var_select == "rand") {
            if (!is.null(ranges)) {
                message("Random CpGs within provided GRanges will be used for the PCA")
                ids <- sample(x = seq_along(meth_sub), replace = FALSE,
                  size = min(top_var, nrow(meth_sub)))
            } else {
                message("Random CpGs will be used for the PCA")
                ids <- sample(x = seq_along(m), replace = FALSE, size = as.integer(as.character(min(top_var,
                  nrow(m)))))
            }
            meth_sub <- get_matrix(m = m[ids, ], type = "M", add_loci = FALSE)
        } else {
            if (!is.null(ranges)) {
                if (is_h5(m)) {
                  sds <- DelayedMatrixStats::rowSds(meth_sub, na.rm = TRUE)
                } else {
                  sds <- matrixStats::rowSds(meth_sub, na.rm = TRUE)
                }
                meth_sub <- meth_sub[order(sds, decreasing = TRUE)[seq_len(min(top_var,
                  nrow(meth_sub)))], ]
            } else {
                meth_sub <- methrix::get_matrix(m = order_by_sd(m)[seq_len(min(top_var,
                  nrow(m)))], type = "M", add_loci = FALSE)
            }
        }
    }

    # Remove NA
    meth_sub <- meth_sub[complete.cases(meth_sub), , drop = FALSE]
    if (nrow(meth_sub) == 0) {
        stop("Zero loci available post NA removal :(")
    }

    meth_pca <- prcomp(x = t(meth_sub), retx = TRUE)

    n_pc <- ncol(meth_pca$x)

    # Variance explained by PC's
    pc_vars <- meth_pca$sdev^2/sum(meth_pca$sdev^2)
    names(pc_vars) <- colnames(meth_pca$x)
    pc_vars <- round(pc_vars, digits = 2)

    #-----------------------------------------------------------------------------------------------------------------------
    # Draw cumulative variance explained by PCs
    
    if (do_plot) {
        par(bty = "n", mgp = c(2.5, 0.5, 0), mar = c(3, 4, 2, 2) + 0.1,
            tcl = -0.25, las = 1)
        plot(pc_vars, type = "h", col = "red", xlab = "", ylab = "variance Explained",
            ylim = c(0, 1), yaxs = "i")
        mtext(side = 1, "Principal component", line = 2)
        cum_var <- cumsum(meth_pca$sdev^2)/sum(meth_pca$sdev^2) * meth_pca$sdev[1]^2/sum(meth_pca$sdev^2)
        lines(cumsum(cum_var), type = "s")
        axis(side = 4, at = pretty(c(0, 1)), labels = pretty(c(0, 1)))
        legend("topright", col = c("red", "black"), lty = 1, c("Per PC", "Cumulative"), bty = "n")
        #lines(x = c(length(meth_pca$sdev), n_pc, n_pc), y = c(cum_var[n_pc], cum_var[n_pc], 0), lty = 3)
        title(main = paste0("Variance explained by ", n_pc, " PC: ", round(sum(c(meth_pca$sdev^2/sum(meth_pca$sdev^2))[seq_len(n_pc)]),
            digits = 2)), adj = 0)
    }
    #-----------------------------------------------------------------------------------------------------------------------

    results <- list(PC_matrix = meth_pca$x, var_explained = pc_vars)

    if (do_plot) {
        plot_pca(pca_res = results, m = m, col_anno = pheno)
    }

    gc(verbose = FALSE)

    return(results)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot PCA results
#'
#' @param pca_res Results from \code{\link{methrix_pca}}
#' @param m optinal methrix object. Default NULL
#' @param col_anno Column name of colData(m). Default NULL. Will be used as a factor to color different groups. Required \code{methrix} object
#' @param shape_anno Column name of colData(m). Default NULL. Will be used as a factor to shape different groups. Required \code{methrix} object
#' @param pc_x Default 'PC1'
#' @param pc_y Default 'PC2'
#' @param show_labels Default FLASE
#' @return ggplot2 object
#' @examples
#' data('methrix_data')
#' mpc = methrix_pca(methrix_data, do_plot = FALSE)
#' plot_pca(mpc)
#' @export
plot_pca <- function(pca_res, m = NULL, col_anno = NULL, shape_anno = NULL,
    pc_x = "PC1", pc_y = "PC2", show_labels = FALSE) {


    X <- Y <- color_me <- shape_me <- row_names <- NULL


    pc_vars <- pca_res$var_explained
    pca_res <- as.data.frame(pca_res$PC_matrix)
    pca_res$row_names <- rownames(pca_res)

    x_lab <- paste0(pc_x, " [", pc_vars[pc_x], " %]")
    y_lab <- paste0(pc_x, " [", pc_vars[pc_y], " %]")

    if (!is.null(col_anno) || !is.null(shape_anno)) {
        if (!is(object = m, class2 = "methrix")) {
            stop("Please provde methrix object while using col_anno or shape_anno")
        }
        pd <- as.data.frame(colData(m))
        pd <- pd[rownames(pca_res), , drop = FALSE]
        pca_res <- cbind(pca_res, pd)
    }

    if (!is.null(col_anno)) {
        col_anno_idx <- which(colnames(pca_res) == col_anno)
        if (length(col_anno_idx) == 0) {
            stop(paste0(col_anno, " not found in provided methrix object"))
        } else {
            colnames(pca_res)[col_anno_idx] <- "color_me"
        }
    }

    if (!is.null(shape_anno)) {
        shape_anno_idx <- which(colnames(pca_res) == shape_anno)
        if (length(shape_anno_idx) == 0) {
            stop(paste0(shape_anno, " not found in provided methrix object"))
        } else {
            colnames(pca_res)[shape_anno_idx] <- "shape_me"
        }
    }

    pc_x_idx <- which(colnames(pca_res) == pc_x)
    pc_y_idx <- which(colnames(pca_res) == pc_y)
    colnames(pca_res)[c(pc_x_idx, pc_y_idx)] <- c("X", "Y")

    if (all(c("color_me", "shape_me") %in% colnames(pca_res))) {
        pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
            shape = shape_me, label = row_names)) + geom_point(size = 3) +
            xlab(pc_x) + ylab(pc_y) + labs(color = col_anno, shape = shape_anno) +
            scale_color_brewer(palette = "Dark2")
    } else if ("color_me" %in% colnames(pca_res)) {
        pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
            label = row_names)) + geom_point(size = 3) + xlab(pc_x) + ylab(pc_y) +
            labs(color = col_anno) + scale_color_brewer(palette = "Dark2")
    } else if ("shape_me" %in% colnames(pca_res)) {
        pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, shape = shape_me,
            label = row_names)) + geom_point(size = 3) + xlab(pc_x) + ylab(pc_y) +
            labs(shape = shape_anno)
    } else {
        pca_gg <- ggplot(data = as.data.frame(pca_res), aes(x = X, y = Y,
            label = row_names)) + geom_point(size = 3, fill = "black",
            color = "gray70") + xlab(pc_x) + ylab(pc_y)
    }

    pca_gg <- pca_gg + xlab(label = x_lab) + ylab(label = y_lab) + theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(colour = "black", size = 12),
            axis.text.y = element_text(colour = "black", size = 12))

    if (show_labels) {
        pca_gg <- pca_gg + geom_label(size = 4)
    }

    pca_gg
}


#--------------------------------------------------------------------------------------------------------------------------
#' Coverage QC Plots
#'
#' @param m Input \code{\link{methrix}} object
#' @param type Choose between 'hist' (histogram) or 'dens' (density plot).
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the plot.
#' @param perGroup Color the plots in a sample-wise manner?
#' @param lim Maximum coverage value to be plotted.
#' @param size.lim The maximum number of observarions (sites*samples) to use. If the dataset is larger that this,
#' random sites will be selected from the genome.
#' @param col_palette Name of the RColorBrewer palette to use for plotting.
#' @return ggplot2 object
#' @examples
#' data('methrix_data')
#' plot_coverage(m = methrix_data)
#' @export

plot_coverage <- function(m, type = c("hist", "dens"), pheno = NULL, perGroup = FALSE,
    lim = 100, size.lim = 1e+06, col_palette = "RdYlGn") {

    value <- variable <- NULL

    if (!is(m, "methrix")){
        stop("A valid methrix object needs to be supplied.")
    }
    colors_palette <- get_palette(ncol(m), col_palette)


    type <- match.arg(arg = type, choices = c("hist", "dens"), several.ok = FALSE)

    if (nrow(m) > size.lim) {
        message("The dataset is bigger than the size limit. A random subset of the object will be used that contains ~",
            size.lim, " observations.")
        n_rows <- trunc(size.lim/nrow(m@colData))
        sel_rows <- sample(seq_len(nrow(m@elementMetadata)), size = n_rows,
            replace = FALSE)

        meth_sub <- methrix::get_matrix(m = m[sel_rows, ], type = "C",
            add_loci = FALSE)

    } else {
        meth_sub <- methrix::get_matrix(m = m, type = "C", add_loci = FALSE)
    }

    if (perGroup) {
        if (is.null(pheno)) {
            stop("For group based plotting, provide group information using the pheno argument.")
        }
        if (pheno %in% colnames(colData(m)) == 0) {
            stop("Phenotype annotation cannot be found in colData(m).")
        }
        colnames(meth_sub) <- m@colData[, pheno]
    }

    meth_sub <- as.data.frame(meth_sub)
    data.table::setDT(x = meth_sub)
    plot.data <- suppressWarnings(expr = data.table::melt(meth_sub))

    plot.data <- plot.data[value <= lim, ]

    # generate the plots
    if (!perGroup) {
        if (type == "dens") {
            p <- ggplot2::ggplot(plot.data, aes(value, color = variable)) +
                ggplot2::geom_density(alpha = 0.5, adjust = 1.5, lwd = 1,
                  position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
                ggplot2::scale_fill_manual(values = colors_palette)

        } else if (type == "hist") {
            p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = variable)) + 
                ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black") + ggplot2::theme_classic() +
                ggplot2::xlab("Coverage")+
                ggplot2::scale_fill_manual(values = colors_palette)
            # print(p)
        }
    } else {
        if (type == "dens") {
            p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, color = variable)) +
                ggplot2::geom_density(alpha = 0.6, adjust = 1.5, lwd = 1,
                  position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
                ggplot2::labs(fill = "Groups") +
                ggplot2::scale_fill_manual(values = colors_palette)
            # print(p)
        } else if (type == "hist") {
            p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = variable)) +
                ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black") + 
                ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
                ggplot2::labs(fill = "Groups") +
                ggplot2::scale_fill_manual(values = colors_palette)
            # print(p)
        }
    }

    gc(verbose = FALSE)

    p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
        colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_blank(), legend.title = element_blank())

}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot descriptive statistics
#' @details plot descriptive statistics results from \code{\link{get_stats}}
#' @param plot_dat results from \code{\link{get_stats}}
#' @param what Can be \code{M} or \code{C}. Default \code{M}
#' @param stat Can be \code{mean} or \code{median}. Default \code{mean}
#' @param ignore_chr Chromsomes to ignore. Default \code{NULL}
#' @param samples Use only these samples. Default \code{NULL}
#' @param n_col number of columns. Passed to `facet_wrap`
#' @param n_row number of rows. Passed to `facet_wrap`
#' @return ggplot2 object
#' @seealso \code{\link{get_stats}}
#' @examples
#' data('methrix_data')
#' gs = get_stats(methrix_data)
#' plot_stats(gs)
#' @export
#'
plot_stats <- function(plot_dat, what = "M", stat = "mean", ignore_chr = NULL,
    samples = NULL, n_col = NULL, n_row = NULL) {

    Chromosome <- . <- Sample_Name <- mean_meth <- sd_meth <- median_meth <- mean_cov <- sd_cov <- NULL
    median_cov <- measurement <- sd_low <- sd_high <- NULL
    what <- match.arg(arg = what, choices = c("M", "C"))
    stat <- match.arg(arg = stat, choices = c("mean", "median"))

    if ("Chromosome" %in% colnames(plot_dat)) {
        if (what == "M") {
            if (stat == "mean") {
                plot_dat <- plot_dat[, .(Chromosome, Sample_Name, mean_meth,
                  sd_meth)]
            } else {
                plot_dat <- plot_dat[, .(Chromosome, Sample_Name, median_meth,
                  sd_meth)]
            }
        } else {
            if (stat == "mean") {
                plot_dat <- plot_dat[, .(Chromosome, Sample_Name, mean_cov,
                  sd_cov)]
            } else {
                plot_dat <- plot_dat[, .(Chromosome, Sample_Name, median_cov,
                  sd_cov)]
            }
        }

        if (!is.null(ignore_chr)) {
            plot_dat <- plot_dat[!Chromosome %in% ignore_chr]
        }

        if (!is.null(samples)) {
            plot_dat <- plot_dat[Sample_Name %in% samples]
        }

        colnames(plot_dat) <- c("Chromosome", "Sample_Name", "measurement",
            "sd")
        plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
        plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
        plot_dat[, `:=`(sd_low, measurement - sd)]
        plot_dat[, `:=`(sd_high, measurement + sd)]
        plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
            no = plot_dat$sd_low)

        plot_dat_gg <- ggplot(data = plot_dat, aes(x = Chromosome, y = measurement)) +
            geom_errorbar(aes(ymin = sd_low, ymax = sd_high), col = "gray70") +
            geom_point(col = "maroon") + facet_wrap(~Sample_Name, nrow = n_row,
            ncol = n_col) + theme_minimal(base_size = 12) + theme(axis.title.x = element_blank(),
            axis.text.x = element_text(hjust = 1, size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_blank())

    } else {
        if (what == "M") {
            if (stat == "mean") {
                plot_dat <- plot_dat[, .(Sample_Name, mean_meth, sd_meth)]
                plot_title <- "Mean methylation"
            } else {
                plot_dat <- plot_dat[, .(Sample_Name, median_meth, sd_meth)]
                plot_title <- "Median methylation"
            }
        } else {
            if (stat == "mean") {
                plot_dat <- plot_dat[, .(Sample_Name, mean_cov, sd_cov)]
                plot_title <- "Mean coverage"
            } else {
                plot_dat <- plot_dat[, .(Sample_Name, median_cov, sd_cov)]
                plot_title <- "Median coverage"
            }
        }

        colnames(plot_dat) <- c("Sample_Name", "measurement", "sd")
        plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
        plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
        plot_dat[, `:=`(sd_low, measurement - sd)]
        plot_dat[, `:=`(sd_high, measurement + sd)]
        plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
            no = plot_dat$sd_low)

        plot_dat_gg <- ggplot(data = plot_dat, aes(x = Sample_Name, y = measurement)) +
            geom_point(col = "maroon", size = 2) + geom_errorbar(aes(ymin = sd_low,
            ymax = sd_high), col = "gray70") + geom_point(col = "maroon") +
            theme_minimal(base_size = 12) + theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                colour = "black"), axis.text.y = element_text(size = 12,
                colour = "black"), axis.title.y = element_blank()) + ggtitle(label = plot_title)
    }

    plot_dat_gg
}
