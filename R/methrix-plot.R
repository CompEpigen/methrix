
#' Violin Plot for \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' data("mm9_bsmap")
#' methrix_violin(m = mm9_bsmap)
methrix_violin <- function(m, ranges = NULL, pheno = NULL, return_ggplot = FALSE){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub = get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
    }else{
      meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
    }

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

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(x = data, y = value, fill = data))+
    ggplot2::geom_violin(alpha = .5)+
    ggplot2::theme_classic()+
    ggplot2::xlab(pheno)+
    ggplot2::ylab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")

  p
}


#' Density Plot of \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param bw.adjust Multiplicate bandwide adjustment. See \code{\link{geom_density}} for more information
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data("mm9_bsmap")
#' methrix_density(m = mm9_bsmap)
methrix_density <- function(m, ranges = NULL, pheno = NULL, return_ggplot = FALSE, bw.adjust = 2){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub = get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
  }else{
    meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
  }

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

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = value, fill = data))+
    ggplot2::geom_density(alpha = .5, adjust = bw.adjust)+
    ggplot2::theme_classic()+
    #ggplot2::xlab(pheno)+
    ggplot2::xlab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")

  p
}


#' Principal Component Analysis
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GRanges}} object
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param do.lot Should a plot be generated?
#' @param n_pc Number of principal components to return.
#' @import FactoMineR
#' @return PCA results
#' @export
#'
methrix_pca <- function(m, ranges = NULL, pheno = NULL, do.plot = TRUE, n_pc = 5){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub = methrix::get_matrix(m = meth_sub, type = "M", add_loci = FALSE)
  }else{
    meth_sub <- methrix::get_matrix(m = m, type = "M", add_loci = FALSE)
  }

  ## melt the object to a long format
  ##To-do: Use irlba for PCA (https://cran.r-project.org/web/packages/irlba/index.html)
  meth.pca <- FactoMineR::PCA(t(meth_sub),graph=F,ncp = n_pc)

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
  plot.data <- data.frame(
    "PC1" = meth.pca$ind$coord[,1],
    "PC2" = meth.pca$ind$coord[,2],
    "ind" = rownames(meth.pca$ind$coord),
    "labs" = pheno.plot$data, stringsAsFactors = FALSE
  )

  #generate the plot
  if(do.plot){
    p <- ggplot2::ggplot(plot.data,ggplot2::aes(x = PC1, y = PC2, color = labs))+
      ggplot2::geom_point(alpha = .5)+
      ggplot2::theme_classic()+
      ggplot2::xlab(paste0("PC1 (",round(meth.pca$eig[1,2],2),"% Variability)"))+
      ggplot2::ylab(paste0("PC2 (",round(meth.pca$eig[2,2],2),"% Variability)"))+
      ggplot2::ggtitle(pheno)+
      ggplot2::labs(color = "Annotation")
    p
  }

  return(meth.pca)
}



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
#' data("mm9_bsmap")
#' methrix_coverage(m = mm9_bsmap)
methrix_coverage <- function(m, type = c("hist","dens"), pheno = NULL, perSample = FALSE, return_ggplot = TRUE, lim = 100){

  type = match.arg(arg = type, choices = c("hist", "dens"), several.ok = FALSE)
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
  plot.data <- plot.data[value <= 100,]

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

  p

}

