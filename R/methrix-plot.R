
#' Violin Plot for \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges \code{\link{GRanges}} object with GenomicRegions which should be used for plotting.
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param return_ggplot Return a \code{\link{ggplot2}} object.
#'
#' @return
#' @export
#'
#' @examples
methrix_violin <- function(m,ranges=NULL,pheno=NULL,return_ggplot=F){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m,ranges)
    meth_sub = methrix::get_matrix(m = meth_sub, type = "M", add_loci = F)
    }else{
      meth_sub <- methrix::get_matrix(m = m, type = "M", add_loci = F)
    }

  ## melt the object to a long format
  meth.melt <- reshape2::melt(meth_sub)

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno%in%colnames(colData(m))==0){stop("Phenotype annotation cannot be found in colData(m).")
      }
      pheno.plot <- data.table("id"=rownames(colData(m)),
                                 "data"=as.factor(colData(m)[,pheno]))
  }else{
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=rownames(colData(m)))
  }

  # merge the pheno column to the others
  plot.data <- merge(meth.melt,pheno.plot,by.x="Var2",by.y="id",all.x=T,all.y=T)

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(data,value,fill=data))+
    ggplot2::geom_violin(alpha=.5)+
    ggplot2::theme_classic()+
    ggplot2::xlab(pheno)+
    ggplot2::ylab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")
  print(p)
  if(return_ggplot==T){
    return(p)
  }
}


#' Density Plot of \eqn{\beta}-Values
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges \code{\link{GRanges}} object with GenomicRegions which should be used for plotting.
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param return_ggplot Return a \code{\link{ggplot2}} object.
#' @param bw.adjust Multiplicate bandwide adjustment. See \code{\link{geom_density}} for more information
#'
#' @return
#' @export
#'
#' @examples
methrix_density <- function(m,ranges=NULL,pheno=NULL,return_ggplot=F,bw.adjust=2){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m,ranges)
    meth_sub = methrix::get_matrix(m = meth_sub, type = "M", add_loci = F)
  }else{
    meth_sub <- methrix::get_matrix(m = m, type = "M", add_loci = F)
  }

  ## melt the object to a long format
  meth.melt <- reshape2::melt(meth_sub)

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno%in%colnames(colData(m))==0){stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=as.factor(colData(m)[,pheno]))
  }else{
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=rownames(colData(m)))
  }

  # merge the pheno column to the others
  plot.data <- merge(meth.melt,pheno.plot,by.x="Var2",by.y="id",all.x=T,all.y=T)

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(value,fill=data))+
    ggplot2::geom_density(alpha=.5,adjust=bw.adjust)+
    ggplot2::theme_classic()+
    #ggplot2::xlab(pheno)+
    ggplot2::xlab(expression(beta*"-Value"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
    ggplot2::labs(fill = "Annotation")
  print(p)
  if(return_ggplot==T){
    return(p)
  }
}


#' Principal Component Analysis
#'
#' @param m Input \code{\link{methrix}} object
#' @param ranges \code{\link{GRanges}} object with GenomicRegions which should be used for subsetting.
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param do.plot Should a plot be generated?
#' @param n_pc Number of principal components to return.
#'
#' @return
#' @export
#'
#' @examples
methrix_pca <- function(m,ranges=NULL,pheno=NULL,do.plot=T,n_pc=5){

  ## subset based on the input ranges
  if(!is.null(ranges)){
    meth_sub <- subset_methrix(m,ranges)
    meth_sub = methrix::get_matrix(m = meth_sub, type = "M", add_loci = F)
  }else{
    meth_sub <- methrix::get_matrix(m = m, type = "M", add_loci = F)
  }

  ## melt the object to a long format
  meth.pca <- FactoMineR::PCA(t(meth_sub),graph=F,ncp = n_pc)

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno%in%colnames(colData(m))==0){stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=as.factor(colData(m)[,pheno]))
  }else{
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=rownames(colData(m)))
  }

  # build the data frame for plotting
  plot.data <- data.frame(
    "PC1" = meth.pca$ind$coord[,1],
    "PC2" = meth.pca$ind$coord[,2],
    "ind" = rownames(meth.pca$ind$coord),
    "labs"=pheno.plot$data
  )

  #generate the violin plot
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(PC1,PC2,color=labs))+
    ggplot2::geom_point(alpha=.5)+
    ggplot2::theme_classic()+
    ggplot2::xlab(paste0("PC1 (",round(meth.pca$eig[1,2],2),"% Variability)"))+
    ggplot2::ylab(paste0("PC2 (",round(meth.pca$eig[2,2],2),"% Variability)"))+
    ggplot2::ggtitle(pheno)+
    ggplot2::labs(color = "Annotation")

  if(do.plot==T){
    print(p)}
    return(meth.pca)
}



#' Covergae QC Plots
#'
#' @param m Input \code{\link{methrix}} object
#' @param type Choose between "hist" (histogram) or "dens" (density plot).
#' @param pheno Column name of colData(m). Will be used as a factor to color different groups in the violin plot.
#' @param perSample Color the plots in a sample-wise manner?
#' @param return_ggplot Return a \code{\link{ggplot2}} object.
#' @param lim Maximum coverage value to be plotted.
#'
#' @return
#' @export
#'
#' @examples
methrix_coverage <- function(m, type=c("hist","dens"),pheno=NULL,perSample=F,return_ggplot=T,lim=100){
  meth_sub <- methrix::get_matrix(m = m, type = "C", add_loci = F)

  ## melt the object to a long format
  meth.melt <- reshape2::melt(meth_sub)

  # add the pheno column of choice
  if(!is.null(pheno)){
    if(pheno%in%colnames(colData(m))==0){stop("Phenotype annotation cannot be found in colData(m).")
    }
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=as.factor(colData(m)[,pheno]))
  }else{
    pheno.plot <- data.table("id"=rownames(colData(m)),
                             "data"=rownames(colData(m)))
  }

  # merge the pheno column to the others
  plot.data <- merge(meth.melt,pheno.plot,by.x="Var2",by.y="id",all.x=T,all.y=T)
  plot.data <- plot.data[plot.data$value<=100,]

  #generate the violin plot
  if(perSample==F){
    if(type=="dens"){
  p <- ggplot2::ggplot(plot.data,ggplot2::aes(value))+
    ggplot2::geom_density(alpha=.5,adjust=1.5)+
    ggplot2::theme_classic()+
    ggplot2::xlab("Coverage")#+
    #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  print(p)
  }else if(type=="hist"){
    p <- ggplot2::ggplot(plot.data,ggplot2::aes(value))+
      ggplot2::geom_histogram(alpha=.5,binwidth = 1)+
      ggplot2::theme_classic()+
      ggplot2::xlab("Coverage")
    print(p)
  }}else{
    if(type=="dens"){
      p <- ggplot2::ggplot(plot.data,ggplot2::aes(value,fill=data))+
        ggplot2::geom_density(alpha=.5,adjust=1.5)+
        ggplot2::theme_classic()+
        #ggplot2::xlab(pheno)+
        ggplot2::xlab("Coverage")+
        #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
        ggplot2::labs(fill = "Samples")
      print(p)}else if(type=="hist"){
        p <- ggplot2::ggplot(plot.data,ggplot2::aes(value,fill=data))+
          ggplot2::geom_histogram(alpha=.5,binwidth = 1)+
          ggplot2::theme_classic()+
          ggplot2::xlab("Coverage")+
          ggplot2::labs(fill = "Samples")
        print(p)
      }
  }
  if(return_ggplot==T){
    return(p)
  }

}

