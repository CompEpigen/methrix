

library(HDF5Array)
library(plotly)



methrix_pca <- function(obj, include=vector(), color,  ...){
  data=as.data.frame(pData(obj))
  p <- autoplotly(prcomp(t(as.matrix(obj@assays[["M"]][include,]))), data = data,
                  colour = color,  ...)

  p
}



methrix.subset <- function(obj=NULL, based_on=c("sd", "random",  "region", "no zeros"), count=NULL, regions=NULL){

  if (class(obj)!="BSseq")
    stop("Input must be a bsseq object.", call. = F)
  if (!(based_on %in% c("sd", "random",  "region", "no zeros")))
    stop("Filtering should be based on sd, random or region", call. = F)
  if (based_on %in% c("sd", "random") & is.null(count))
    stop("You have to specified the number of sites to include", call. = F)
  if (based_on %in% c( "region") & !is.null(count))
    warning(paste0("You are filtering based on: ", based_on, ". The number of sites to included won't be used."))
  if (based_on=="region" & !is.null(regions))
    stop("You have to define the regions.")
  if (!is.null(regions) && class(regions)!=GRanges)
    stop("The regions must be in GRanges format.")
  ##add more
#browser()
  if (based_on=="sd"){
    sds <- DelayedMatrixStats::rowSds(obj@assays[["M"]])
    obj <- obj[which(sds >= sds[order(sds, decreasing = T, na.last =  NA)][count]),]
    return(obj)
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



########filtering methods
#vector means a vector of integers


metrix.filter <- function(obj=NULL, based_on=c("vector", "not covered", "region"), vector= NULL, regions=NULL){

  if (class(obj)!="BSseq")
    stop("Input must be a bsseq object.", call. = F)
  if (!(based_on %in% c("vector", "not covered", "region")))
    stop("Filtering should be based on sd, random, not covered or region", call. = F)
  if (!is.null(regions) && class(regions)!=GRanges)
    stop("The regions must be in GRanges format.")
  if (!isTRUE(all(vector == floor(vector)))) stop("The supplied vector  must only contain integer values")


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
