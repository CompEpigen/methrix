#' Class methrix
#' @description S4 class Methrix
#' @slot SE summarized experiment object
#' @slot CpGs number of CpGs
#' @slot samples number of samples
#' @slot h5 are matrices are HDF5Array
#' @exportClass methrix

methrix <- setClass(Class = "methrix",
                       slots=c(SE = "SummarizedExperiment", CpGs = "numeric", samples = "numeric", h5 = "logical"))

setMethod(f = 'show', signature(object = "methrix"), definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(data.table::data.table(ID = c("n_CpGs", "n_samples", "is_HDF5"),
                               summary = c(object@CpGs, object@samples, as.logical(object@h5))))
})
