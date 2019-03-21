#' Class methrix
#' @description S4 class Methrix
#' @slot beta beta matrix
#' @slot cov coverage matrix
#' @slot metadata genome: the name of the BSgenome that was used to extract CpGs, isHDF5: is it stored in HDF5 Array format
#' @slot sample_annotation sample annotation data frame
#' @slot position chromosome, start and strand information in a data.table format
#' @exportClass methrix
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'

# we can include validity checks as well
# we might need to keep track of the processing steps

.methrix <- setClass("methrix", contains="SummarizedExperiment")

methrix <- function(beta, cov, position, is.HDF5=NULL, genome, sample_annotation) {

  se <- SummarizedExperiment(assays=list(beta=beta, cov=cov), metadata=list(genome=genome, is.HDF5=is.HDF5),
                             colData= sample_annotation, rowData=position)

  .methrix(se)
}



# #' Class methrix
# #' @description S4 class Methrix
# #' @slot SE summarized experiment object
# #' @slot CpGs number of CpGs
# #' @slot samples number of samples
# #' @slot h5 are matrices are HDF5Array
# #' @exportClass methrix
#
# methrix <- setClass(Class = "methrix",
#                        slots=c(SE = "SummarizedExperiment", CpGs = "numeric", samples = "numeric", h5 = "logical", genome = "character"))
#
# setMethod(f = 'show', signature(object = "methrix"), definition = function(object){
#   cat(paste('An object of class ', class(object), "\n"))
#   print(data.table::data.table(ID = c("n_CpGs", "n_samples", "is_HDF5"),
#                                summary = c(object@CpGs, object@samples, as.logical(object@h5))))
# })

# #' Class methrixDT
# #' @description S4 class Methrix
# #' @slot betas data.table of beta values
# #' @slot covs data.table of coverage values
# #' @slot chr chr loci
# #' @slot coldata pheno dtata
# #' @exportClass methrixDT
#
# methrix <- setClass(Class = "methrixDT",
#                     slots=c(betas = "data.table", covs = "data.table", chr = "data.table", coldata = "data.table"))
#
# setMethod(f = 'show', signature(object = "methrix"), definition = function(object){
#   cat(paste('An object of class ', class(object), "\n"))
#   print(data.table::data.table(ID = c("n_CpGs", "n_samples", "is_HDF5"),
#                                summary = c(object@CpGs, object@samples, as.logical(object@h5))))
# })
