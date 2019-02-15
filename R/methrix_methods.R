#' extract methylation matrix
#' @name getMethylationMatrix
#' @rdname getMethylationMatrix
#' @param x An object of class methrix
#' @return methylation methrix
#' @exportMethod getMethylationMatrix
setGeneric(name = "getMethylationMatrix", function(x) standardGeneric("getMethylationMatrix"))

## Accessor methods
#' @rdname getMethylationMatrix
#' @aliases getMethylationMatrix
setMethod(f = "getMethylationMatrix", signature = "methrix", function(x){
  assay(x@SE)
})
