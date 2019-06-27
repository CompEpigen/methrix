#' Class methrix
#' @description S4 class Methrix
#' @slot beta beta matrix
#' @slot cov coverage matrix
#' @slot metadata genome: the name of the BSgenome that was used to extract CpGs, isHDF5: is it stored in HDF5 Array format
#' @slot sample_annotation sample annotation data frame
#' @slot position chromosome, start and strand information in a data.table format
#' @exportClass methrix
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'

methrix <- setClass(Class = 'methrix', contains = "SummarizedExperiment")
# we can include validity checks as well
# we might need to keep track of the processing steps

setMethod(f = 'show', signature = 'methrix', definition = function(object){
  cat(paste0('An object of class ', class(object), "\n"))
  cat(paste0("   n_CpGs: ", format(nrow(object), big.mark = ","), "\n"))
  cat(paste0("n_samples: ", ncol(object), "\n"))
  cat(paste0("    is_h5: ", is_h5(object), "\n"))
  cat(paste0("Reference: ", object@metadata$genome, "\n"))
})

#Create methrix obj
create_methrix = function(beta_mat = NULL, cov_mat = NULL, cpg_loci = NULL, is_hdf5 = FALSE,
                          genome_name = "hg19", col_data = NULL, h5_dir = NULL, ref_cpg_dt = NULL, chrom_sizes = NULL, desc = NULL){

  if(is_hdf5){
    se = SummarizedExperiment::SummarizedExperiment(assays = list(beta = as(beta_mat, "HDF5Array"), cov = as(cov_mat, "HDF5Array")),
                                                    metadata = list(genome = genome_name, is_h5 = is_hdf5, ref_CpG = ref_cpg_dt,
                                                                    chrom_sizes = chrom_sizes, descriptive_stats = desc),
                                                    colData = col_data, rowData = cpg_loci)
    if(!is.null(h5_dir)){
      tryCatch(HDF5Array::saveHDF5SummarizedExperiment(x = se, dir = h5_dir, replace = TRUE),
               error = function(e) message("The dataset is not saved. Please save manually, using the HDF5Array::saveSummarizedExperiment command. "))
    }
  }else{
    se = SummarizedExperiment::SummarizedExperiment(assays = list(beta = as.matrix(beta_mat), cov = as.matrix(cov_mat)),
                                                    metadata = list(genome = genome_name, is_h5 = is_hdf5, ref_CpG = ref_cpg_dt,
                                                                    chrom_sizes = chrom_sizes, descriptive_stats = desc),
                                                    colData = col_data, rowData = cpg_loci)
  }

  return(methrix(se))
}
