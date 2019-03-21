
#' subset methrix objects
#' @name methrix.subset
#' @rdname methrix.subset
#' @param obj An object of class methrix
#' @param based_on The subsetting should be based on either high standard deviation ("sd"), random selection ("random"), regions ("region") or sites with no zeros ("no zeros").
#' @param count The number of sites to include, if the subsetting is based on sd or random
#' @param regions The regions to include as a GRanges object
#' @param include_calc which samples should be included into the sd-based subsetting
#' @return methrix object
#' @exportMethod methrix.subset
#' @import HDF5Array
#' @import DelayedMatrixStats




# methrix_pca <- function(obj, color,  ...){
#   data=as.data.frame(pData(obj))
#   p <- p <- autoplotly(prcomp(t(as.matrix(obj@assays[["M"]]))), data = data,
#                   colour = color,  ...)
#
#   p
# }



methrix.subset <- function(obj=NULL, based_on=c("sd", "random",  "region", "no zeros"), count=NULL, regions=NULL, include_calc=NULL){

  if (class(obj)!="BSseq")
    stop("Input must be a bsseq object.", call. = F)
  if (!(based_on %in% c("sd", "random",  "region", "no zeros")))
    stop("Filtering should be based on sd, random or region", call. = F)
  if (based_on %in% c("sd", "random") & is.null(count))
    stop("You have to specified the number of sites to include", call. = F)
  if (based_on %in% c( "region") & !is.null(count))
    warning(paste0("You are filtering based on: ", based_on, ". The number of sites to included won't be used."))
  if (based_on=="region" & is.null(regions))
    stop("You have to define the regions.")
  if (!is.null(regions) && class(regions)!="GRanges")
    stop("The regions must be in GRanges format.")
  ##add more
#browser()
  if (based_on=="sd"){
    #browser()
    if (is.null(include_calc)){
    sds <- DelayedMatrixStats::rowSds(obj@assays[["M"]])
    obj <- obj[which(sds >= sds[order(sds, decreasing = T, na.last =  NA)][count]),]
    return(obj)} else {
      if (!all(include_calc %in% obj@colData@rownames)){
        stop("The provided samples are not in the dataset")
      } else {
        sds <- DelayedMatrixStats::rowSds(obj@assays[["M"]][,which(obj@colData@rownames %in% include_calc)])
      }
      obj <- obj[which(sds >= sds[order(sds, decreasing = T, na.last =  NA)][count]),]
      return(obj)
    }
  }

  if (based_on=="random"){
    rand <- sample(1:length(obj), size = count, replace = F)
    obj <- obj[rand,]
    return(obj)
  }

  if (based_on=="region"){
    #it has to be with fOverlaps with methrix
    obj <- subsetByOverlaps(obj, ranges = regions)
    return(obj)
  }
  if (based_on=="no zeros"){

    obj <- obj[which(!DelayedMatrixStats::rowAnyNAs(obj@assays[["M"]])),]
    return(obj)
  }

}


#' filter methrix objects
#' @name methrix.filter
#' @rdname methrix.filter
#' @param obj An object of class methrix
#' @param based_on The filtering should be based on either a vector of integers ("vector"), sites that are not covered in any of the samples ("not covered") or region to exclude.
#' @param vector Vector of integers to remove
#' @param regions The regions to remove as a GRanges object
#' @return methrix object
#' @exportMethod methrix.subset
#' @import HDF5Array
#' @import DelayedMatrixStats


#vector means a vector of integers


methrix.filter <- function(obj=NULL, based_on=c("vector", "not covered", "region"),
                          vector= NULL, regions=NULL){

  if (class(obj)!="methrix")
    stop("Input must be a methrix object.", call. = F)
  if (!(based_on %in% c("vector", "not covered", "region")))
    stop("Filtering should be based on sd, random, not covered or region", call. = F)
  if (based_on=="region" & is.null(regions))
    stop("You have to define the regions.")
  if (!is.null(regions) && class(regions)!="GRanges")
    stop("The regions must be in GRanges format.")

  if (based_on=="vector" && !is.null(vector) && !isTRUE(all(vector == floor(vector)))) stop("The supplied vector  must only contain integer values")


  if (based_on=="vector"){
    obj <- obj[-vector,]
    return(obj)
  }
  if (based_on=="not covered"){
    not_zero <- which(DelayedMatrixStats::rowSums2(obj@assays[["Cov"]])!=0)
    obj <- obj[not_zero,]
    return(obj)
  }

  if (based_on=="region"){
    obj <- obj[-findOverlaps(obj, regions)@from,]
    return(obj)
  }

}

#' filters methrix objects on individual site level
#' @name methrix.filter_sites
#' @rdname methrix.filter_sites
#' @param obj An object of class methrix
#' @param based_on The filtering should be based on coverage that is above a limit ("highly covered"), or below ("lowly covered")
#' @param cov.limit Coverage limit
#' @return methrix object
#' @exportMethod methrix.filter_sites
#' @import HDF5Array



methrix.filter_sites <- function(obj, based_on=c("highly covered", "lowly covered"), cov.limit=NULL){

  if (based_on=="highly covered"){
    obj@assays[["M"]][obj@assays[["Cov"]]>=cov.limit]<- NA
    obj@assays[["Cov"]][obj@assays[["Cov"]]>=cov.limit]<- 0
    return(obj)
  }
  if (based_on=="lowly covered"){
    obj@assays[["M"]][obj@assays[["Cov"]] <=cov.limit]<- NA
    obj@assays[["Cov"]][obj@assays[["Cov"]]<=cov.limit]<- 0
    return(obj)
  }

}

##impute function

methrix.impute <- function(obj, methiod=c("average", "knn")){

#obj@assays[["M"]] <- apply(obj@assays[["M"]], function(x), )


}

