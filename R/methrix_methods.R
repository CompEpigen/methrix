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
  SummarizedExperiment::assay(x = x@SE, i = 1)
})

#' extract methylation matrix
#' @name getCoverageMatrix
#' @rdname getCoverageMatrix
#' @param x An object of class methrix
#' @return coverage methrix
#' @exportMethod getCoverageMatrix
setGeneric(name = "getCoverageMatrix", function(x) standardGeneric("getCoverageMatrix"))

## Accessor methods
#' @rdname getCoverageMatrix
#' @aliases getCoverageMatrix
setMethod(f = "getCoverageMatrix", signature = "methrix", function(x){
  SummarizedExperiment::assay(x = x@SE, i = 2)
})

#' extract CpG Loci
#' @name getCpgLoci
#' @rdname getCpgLoci
#' @param x An object of class methrix
#' @return An object of class DataFrame
#' @exportMethod getCpgLoci
setGeneric(name = "getCpgLoci", function(x) standardGeneric("getCpgLoci"))

## Accessor methods
#' @rdname getCpgLoci
#' @aliases getCpgLoci
setMethod(f = "getCpgLoci", signature = "methrix", function(x){
  SummarizedExperiment::rowData(x = x@SE)
})
