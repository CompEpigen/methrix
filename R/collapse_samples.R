#' Combines samples based on a column in the colData
#' @details Takes \code{\link{methrix}} object and a column name and combines samples with the same name in the annotataion.
#' @param m \code{\link{methrix}} object
#' @param col a column in the colData, with new sample names. Samples with the same names will be merged and those with NA will be removed
#' @return a \code{\link{methrix}} object
#' @examples
#' data("methrix_data")
#' collapse_samples(m = methrix_data, col="Condition")
#' @export
collapse_samples <- function(m, col=NULL){

  if (!(col %in% colnames(m@colData))){
    stop("The provided column name is not in the dataset. Please provide a name that is in the colData. ")
  }
  m <- m[!is.na(m@colData[,col])]
  groupby <- factor(m@colData[,col])
  sp <- split(seq(along = groupby), groupby)

  if (is_h5(m)){
    stop("Not implemented yet.")
    # new_cov <- sapply(sp, function(i) rowSums(m@assays[[2]][,i, drop = FALSE], na.rm=T))
    # new_beta <- m@assays[[1]]*m@assays[[2]]
    # new_beta <- sapply(sp, function(i) rowSums(new_beta[,i, drop = FALSE], na.rm=T))
    # new_beta <- new_beta/new_cov
  } else {
    new_cov <- sapply(sp, function(i) matrixStats::rowSums2(m@assays[[2]][,i, drop = FALSE], na.rm=TRUE))
    new_beta <- m@assays[[1]]*m@assays[[2]]
    new_beta <- sapply(sp, function(i) matrixStats::rowSums2(new_beta[,i, drop = FALSE], na.rm=TRUE))
    new_beta <- new_beta/new_cov
  }

  new_beta[new_cov==0] <- NA
  new_cov[new_cov==0] <- NA
  combined <- m[,unlist(lapply(sp, '[', 1))]
  dimnames(combined)=dimnames(new_cov)

  combined@assays[[1]] <- new_beta
  combined@assays[[2]] <- new_cov



  combined


}
