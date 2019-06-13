#' Violin Plot for \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param n_cpgs Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' data("methrix_data")
#' methrix_violin(m = methrix_data)
methrix_violin <- function(m, ranges = NULL, n_cpgs = 25000, pheno = NULL){

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno %in% colnames(colData(m)) == 0){
      stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table::data.table("id" = rownames(colData(m)),
                             "data" = as.factor(colData(m)[, pheno]))
  }else{
    pheno.plot <- data.table::data.table("id" = rownames(colData(m)),
                             "data" = rownames(colData(m)))
  }

  ## subset based on the input ranges
  if(!is.null(ranges)){
      meth_sub <-get_matrix(m = subset_methrix(m = m, regions = ranges), type = "M", add_loci = FALSE)
  }else if(is.null(n_cpgs)){
      meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
  }else{
    n_cpgs = as.integer(as.character(n_cpgs))
    if(nrow(m) < n_cpgs){
      n_cpgs = nrow(m)
    }
    set.seed(seed = 1024)
    ids = sample(x = 1:nrow(m), replace = FALSE, size = n_cpgs)
    meth_sub <- get_matrix(m = m[ids, ], type = "M", add_loci = FALSE)
  }

  ## melt the object to a long format
  meth.melt <- data.table::melt(meth_sub)
  data.table::setDT(x = meth.melt)

  # merge the pheno column to the others
  plot.data <- merge(x = meth.melt, y = pheno.plot, by.x = "Var2", by.y = "id", all.x = TRUE, all.y = TRUE)

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(x = data, y = value, fill = data))+
    ggplot2::geom_violin(alpha = .8)+
    ggplot2::theme_classic(base_size = 14)+scale_fill_viridis_d()+
    ggplot2::xlab(pheno)+
    ggplot2::ylab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")

  gc(verbose = FALSE)

  p
}

#--------------------------------------------------------------------------------------------------------------------------
#' Density Plot of \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param n_cpgs Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param bw.adjust Multiplicate bandwide adjustment. See \code{\link{geom_density}} for more information
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data("methrix_data")
#' methrix_density(m = methrix_data)
methrix_density <- function(m, ranges = NULL, n_cpgs = 25000, pheno = NULL, bw.adjust = 2){

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno %in% colnames(colData(m)) == 0){
      stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = as.factor(colData(m)[, pheno]))
  }else{
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = rownames(colData(m)))
  }

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub = get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
  }else if(is.null(n_cpgs)){
    meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
  }else{
    n_cpgs = as.integer(as.character(n_cpgs))
    if(nrow(m) < n_cpgs){
      n_cpgs = nrow(m)
    }
    set.seed(seed = 1024)
    ids = sample(x = 1:nrow(m), replace = FALSE, size = n_cpgs)
    meth_sub <- get_matrix(m = m[ids, ], type = "M", add_loci = FALSE)
  }

  ## melt the object to a long format
  meth.melt <- data.table::melt(meth_sub)
  data.table::setDT(x = meth.melt)

  # merge the pheno column to the others
  plot.data <- merge(x = meth.melt, y = pheno.plot, by.x = "Var2", by.y = "id", all.x = TRUE, all.y = TRUE)

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = value, fill = data))+
    ggplot2::geom_density(alpha = .5, adjust = bw.adjust)+
    ggplot2::theme_classic(base_size = 14)+scale_fill_viridis_d()+
    #ggplot2::xlab(pheno)+
    ggplot2::xlab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")

  gc(verbose = FALSE)

  p
}

#--------------------------------------------------------------------------------------------------------------------------
#' Principal Component Analysis
#'
#' @param m Input \code{\link{methrix}} object
#' @param top_var Number of variable CpGs to use. Default 5000. Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param var Choose between random CpG sites ("rand") or most variable CpGs ("top").
#' @param do_plot Should a plot be generated?
#' @param n_pc Number of principal components to return. Default 5.
#' @param do_fast Use the \code{\link{prcomp_irlba}} function for a quick PCA. This might be useful when computing with large datasets.
#' @return PCA results
#' @importFrom  irlba prcomp_irlba
#' @examples
#' data("methrix_data")
#' methrix_pca(methrix_data)
#' @export
#'
methrix_pca <- function(m, var="top",top_var = 1000, ranges = NULL, pheno = NULL, do_plot = TRUE, n_pc = 5,do_fast=FALSE){
  var_select <- match.arg(var,c("top","rand"))
  ## subset based on the input ranges
  if(!is.null(ranges)){
    print("GenomicRanges will be used for the PCA.")
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub = methrix::get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
  }else if(is.null(top_var)){
    print("All CpGs in the dataset will be used for the PCA.")
    meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
  }else{
    print("Top variable CpGs will be used for the PCA.")
    top_var = as.integer(as.character(top_var))
    if(nrow(m) < top_var){
      top_var = nrow(m)
    }
    if(var_select=="rand"){set.seed(seed = 1024)
    ids = sample(x = 1:nrow(m), replace = FALSE, size = as.integer(as.character(top_var)))
      meth_sub <- get_matrix(m = m[ids, ], type = "M", add_loci = FALSE)
      }else{
      meth_sub <-  methrix::get_matrix(m = m, type = "M", add_loci = FALSE)
      mv <- apply(meth_sub,1,sd)
      mv_ord <- order(mv,decreasing = T)
      meth_sub <-  methrix::get_matrix(m = m[mv_ord[1:5000],], type = "M", add_loci = FALSE)
    }
  }

  #Remove NA
  meth_sub = meth_sub[complete.cases(meth_sub),, drop = FALSE]
  if(nrow(meth_sub) == 0){
    stop("Zero loci available post NA removal :(")
  }

  set.seed(seed = 1024)
  if(do_fast==TRUE){
  meth_pca = irlba::prcomp_irlba(x = t(meth_sub), retx = TRUE,n = n_pc)
  }else{
  meth_pca = prcomp(x = t(meth_sub), retx = TRUE)
  }

  if(n_pc > ncol(meth_pca$x)){
    n_pc = ncol(meth_pca$x)
  }

  #-----------------------------------------------------------------------------------------------------------------------
  # Draw cumulative variance explained by PCs
  if(do_plot==TRUE){
    pc_vars = meth_pca$sdev ^ 2 / sum(meth_pca$sdev ^ 2)
  par(bty = "n", mgp = c(2.5,.5,0), mar = c(3,4,2,2)+.1, tcl = -.25, las = 1)
  plot(
    pc_vars,
    type = "h",
    col = "red",
    xlab = "",
    ylab = "variance Explained" ,
    ylim = c(0, 1),
    yaxs = "i"
  )
  mtext(side = 1, "Principal component", line = 2)
  cum_var =
    cumsum(meth_pca$sdev ^ 2) / sum(meth_pca$sdev ^ 2) * meth_pca$sdev[1] ^
    2 / sum(meth_pca$sdev ^ 2)
  lines(cum_var, type = "s")
  axis(
    4,
    at = pretty(c(0, 1)) * meth_pca$sdev[1] ^ 2 / sum(meth_pca$sdev ^ 2),
    labels = pretty(c(0, 1))
  )
  legend("topright", col = c("red", "black"), lty = 1, c("Per PC","Cumulative"), bty = "n")
  lines(x = c(length(meth_pca$sdev), n_pc, n_pc), y = c(cum_var[n_pc], cum_var[n_pc], 0), lty = 3)
  title(main = paste0("Variance explained by ",  n_pc , " PC: ", round(sum(c(meth_pca$sdev^2/sum(meth_pca$sdev^2))[1:n_pc]), digits = 2)), adj = 0)}
  #-----------------------------------------------------------------------------------------------------------------------

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno %in% colnames(colData(m)) == 0){
      stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = as.factor(colData(m)[, pheno]))
  }else{
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = rownames(colData(m)))
  }

  # build the data frame for plotting
  #return(list(pheno.plot, meth_pca))
  plot.data <- data.frame(
    "PC1" = meth_pca$x[,"PC1"],
    "PC2" = meth_pca$x[,"PC2"],
    "ind" = colnames(meth_pca$rotation),
    "labs" = pheno.plot$data, stringsAsFactors = FALSE
  )

  #generate the plot
  if(do_plot){
    p <- ggplot2::ggplot(plot.data,ggplot2::aes(x = PC1, y = PC2, color = labs, label = ind))+
      ggplot2::geom_point(alpha = .8, size = 3)+scale_color_viridis_d()+
      ggplot2::theme_minimal(base_size = 14)+
      ggplot2::xlab(paste0("PC1 (",round(pc_vars[1]*100, digits = 2),"% Variability)"))+
      ggplot2::ylab(paste0("PC2 (",round(pc_vars[2]*100, digits = 2),"% Variability)"))+
      ggplot2::ggtitle(pheno)+
      ggplot2::labs(color = "Annotation")
    print(p)
  }

  gc(verbose = FALSE)

  return(meth_pca$x)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Covergae QC Plots
#'
#' @param m Input \code{\link{methrix}} object
#' @param type Choose between "hist" (histogram) or "dens" (density plot).
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param perSample Color the plots in a sample-wise manner?
#' @param lim Maximum coverage value to be plotted.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data("methrix_data")
#' methrix_coverage(m = methrix_data)
methrix_coverage <- function(m, type = c("hist","dens"), pheno = NULL, perSample = FALSE, lim = 100){

  type = match.arg(arg = type, choices = c("hist", "dens"), several.ok = FALSE)
  #On an average a matrix of 28e6 rows x 10 columns, sizes around 2.4 GB. Copying, and melting would double the memory consumption.
  #We should think of something else.
  meth_sub <- methrix::get_matrix(m = m, type = "C", add_loci = FALSE)

  ## melt the object to a long format
  meth.melt <- data.table::melt(meth_sub)
  data.table::setDT(x = meth.melt)

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno %in% colnames(colData(m)) == 0){
      stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = as.factor(colData(m)[, pheno]))
  }else{
    pheno.plot <- data.table("id" = rownames(colData(m)),
                             "data" = rownames(colData(m)))
  }
  # merge the pheno column to the others
  plot.data <- merge(x = meth.melt, y = pheno.plot, by.x = "Var2", by.y = "id", all.x = TRUE, all.y = TRUE)
  plot.data <- plot.data[value <= lim,]

  #generate the violin plot
  if(!perSample) {
    if(type == "dens"){
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value)) +
        ggplot2::geom_density(alpha = .5, adjust = 1.5) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage")#+
      #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      #print(p)
    } else if(type == "hist") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value)) +
        ggplot2::geom_histogram(alpha = .5, binwidth = 1) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage")
      #print(p)
    }
  } else{
    if(type == "dens") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = data)) +
        ggplot2::geom_density(alpha = .5, adjust = 1.5) +
        ggplot2::theme_classic() +
        #ggplot2::xlab(pheno)+
        ggplot2::xlab("Coverage") +
        #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
        ggplot2::labs(fill = "Samples")
      #print(p)
    } else if (type == "hist") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = data)) +
        ggplot2::geom_histogram(alpha = .5, binwidth = 1) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage") +
        ggplot2::labs(fill = "Samples")
      #print(p)
    }
  }

  gc(verbose = FALSE)

  p

}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot descriptive statistics
#' @details plot descriptive statistics results from get_stat
#' @param plot_dat results from \code{\link{get_stat}}
#' @param what Can be \code{M} or \code{C}. Default \code{M}
#' @param stat Can be \code{mean} or \code{median}. Default \code{mean}
#' @param ignore_chr Chromsomes to ignore. Default \code{NULL}
#' @param samples Use only these samples. Default \code{NULL}
#' @seealso \code{\link{get_stats}}
#' @examples
#' data("methrix_data")
#' gs = get_stats(methrix_data)
#' plot_stats(gs)
#' @export
#'
plot_stats = function(plot_dat, what = "M", stat = "mean", ignore_chr = NULL, samples = NULL){

  what = match.arg(arg = what, choices = c("M", 'C'))
  stat = match.arg(arg = stat, choices = c("mean", 'median'))

  if("Chromosome" %in% colnames(plot_dat)){
    if(what == "M"){
      if(stat == "mean"){
        plot_dat = plot_dat[,.(Chromosome, Sample_Name, mean_meth, sd_meth)]
      }else{
        plot_dat = plot_dat[,.(Chromosome, Sample_Name, median_meth, sd_meth)]
      }
    }else{
      if(stat == "mean"){
        plot_dat = plot_dat[,.(Chromosome, Sample_Name, mean_cov, sd_cov)]
      }else{
        plot_dat = plot_dat[,.(Chromosome, Sample_Name, median_cov, sd_cov)]
      }
    }

    if(!is.null(ignore_chr)){
      plot_dat = plot_dat[!Chromosome %in% ignore_chr]
    }

    if(!is.null(ignore_chr)){
      plot_dat = plot_dat[!Sample_Name %in% samples]
    }

    colnames(plot_dat) = c("Chromosome", "Sample_Name", "measurement", "sd")
    plot_dat[, sd_low := measurement-sd]
    plot_dat[, sd_high := measurement+sd]
    plot_dat$sd_low = ifelse(test = plot_dat$sd_low < 0, yes = 0, no = plot_dat$sd_low)

    samps = plot_dat[,.N,Sample_Name][,Sample_Name]
    chrs = unique(plot_dat[,Chromosome])
    plot_dat = split(plot_dat, as.factor(as.character(plot_dat$Sample_Name)))

    lo = layout(mat = matrix(data = c(1:(length(samps)+1)), byrow = TRUE, ncol = 1), heights = c(rep(4, length(samps)), 2))

    lapply(plot_dat, function(samp){
      y_dat = get_y_lims(c(samp$sd_low, samp$sd_high))
      y_lims = y_dat$y_lims
      y_at = y_dat$y_at

      par(mar = c(0, 3, 2, 0))
      plot(x = NA, NA, xlim = c(1, length(chrs)), frame.plot = FALSE, axes = FALSE, ylim = y_lims, xlab = NA, ylab = NA)
      abline(h = y_at, v = 1:length(samps), lty = 2, col = grDevices::adjustcolor(col = "gray", alpha.f = 0.4))
      segments(x0 = 1:length(chrs), y0 = samp$sd_low, x1 = 1:length(chrs), y1 = samp$sd_high, col = "gray70")
      points(x = 1:length(chrs), y = samp$measurement, pch = 19, col = "maroon")
      axis(side = 2, at = seq(0, 1, 0.25), las = 2, col = grDevices::adjustcolor("gray", alpha.f = 0.9))
      title(main = samp[,Sample_Name][1], adj = 0, font.main= 3, col.main = "royalblue")
    })

    par(mar = c(0, 3, 0, 0))
    plot(NA, NA, xlim = c(1, length(chrs)), ylim = c(0, 0.1), axes = FALSE, xlab = NA, ylab = NA)
    mtext(text = chrs, side = 1, outer = FALSE, at = 1:length(chrs), las = 2, line = -1.5, cex = 0.8)

  }else{
    if(what == "M"){
      if(stat == "mean"){
        plot_dat = plot_dat[,.(Sample_Name, mean_meth, sd_meth)]
        plot_title = "Mean methylation"
      }else{
        plot_dat = plot_dat[,.(Sample_Name, median_meth, sd_meth)]
        plot_title = "Median methylation"
      }
    }else{
      if(stat == "mean"){
        plot_dat = plot_dat[,.(Sample_Name, mean_cov, sd_cov)]
        plot_title = "Mean coverage"
      }else{
        plot_dat = plot_dat[,.(Sample_Name, median_cov, sd_cov)]
        plot_title = "Median coverage"
      }
    }

    colnames(plot_dat) = c("Sample_Name", "measurement", "sd")
    plot_dat[, sd_low := measurement-sd]
    plot_dat[, sd_high := measurement+sd]
    plot_dat$sd_low = ifelse(test = plot_dat$sd_low < 0, yes = 0, no = plot_dat$sd_low)

    samps = plot_dat[,.N,Sample_Name][,Sample_Name]
    y_dat = get_y_lims(c(plot_dat$sd_low, plot_dat$sd_high))
    y_lims = y_dat$y_lims
    y_at = y_dat$y_at

    par(mar = c(5, 3, 2, 0))
    plot(NA, NA, xlim = c(1, length(samps)), frame.plot = FALSE, axes = FALSE, ylim = y_lims, xlab = NA, ylab = NA)
    abline(h = y_at, v = 1:length(samps), lty = 2, col = grDevices::adjustcolor(col = "gray", alpha.f = 0.4))
    segments(x0 = 1:length(samps), y0 = plot_dat$sd_low, x1 = 1:length(samps), y1 = plot_dat$sd_high, col = "gray70")
    points(x = 1:length(samps), y = plot_dat$measurement, pch = 19, col = "maroon")
    axis(side = 2, at = y_at, las = 2, col = "gray70")
    mtext(text = samps, side = 1, outer = FALSE, at = 1:length(samps), las = 2, line = -0.5, cex = 0.8)
    title(main = plot_title, adj = 0, font.main= 3, col.main = "blue")
  }
}

#Tiny script to get axis and limits
get_y_lims = function(vec){

  y_lims = range(vec)
  y_at = pretty(y_lims)

  if(y_at[1] > min(vec, na.rm = TRUE)){
    y_at[1] = min(vec, na.rm = TRUE)
  }
  if(y_at[length(y_at)] < max(vec, na.rm = TRUE)){
    y_at[length(y_at)] = max(vec, na.rm = TRUE)
  }
  y_lims = range(y_at, na.rm = TRUE)

  list(y_lims = y_lims, y_at = y_at)
}

