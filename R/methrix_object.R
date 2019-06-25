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
  print(object@metadata$summary)
})

#Create methrix obj
create_methrix = function(beta_mat = NULL, cov_mat = NULL, cpg_loci = NULL, is_hdf5 = FALSE,
                          genome_name = "hg19", col_data = NULL, h5_dir = NULL, ref_cpg_dt = NULL, chrom_sizes = NULL, desc = NULL){

  if(is_hdf5){
    n_non_covered = length(which(DelayedMatrixStats::rowSums2(x = cov_mat) == 0))
    n_non_covered = paste0(n_non_covered, " [", round(n_non_covered/nrow(cov_mat) * 100, digits = 2), "%]")
  }else{
    beta_mat = data.table:::as.matrix.data.table(beta_mat)
    cov_mat = data.table:::as.matrix.data.table(cov_mat)
    n_non_covered = length(which(matrixStats::rowSums2(x = cov_mat) == 0))
    n_non_covered = paste0(n_non_covered, " [", round(n_non_covered/nrow(cov_mat) * 100, digits = 2), "%]")
  }


  se_summary = data.table::data.table(ID = c("n_samples", "n_CpGs", "n_uncovered", "n_chromosomes", "Reference_Build", "is_H5"),
                                      Summary = c(ncol(beta_mat), format(nrow(beta_mat), big.mark = ","),
                                                  n_non_covered, nrow(cpg_loci[,.N,chr]), genome_name, is_hdf5))

  chr_summary = cpg_loci[,.N,chr]

  if(is_hdf5){
    se = SummarizedExperiment::SummarizedExperiment(assays = list(beta = as(beta_mat, "HDF5Array"), cov = as(cov_mat, "HDF5Array")),
                                                    metadata = list(genome = genome_name, is_h5 = is_hdf5, summary = se_summary,
                                                                    chr_summary = chr_summary, ref_CpG = ref_cpg_dt,
                                                                    chrom_sizes = chrom_sizes, descriptive_stats = desc),
                                                    colData = col_data, rowData = cpg_loci)
    if(!is.null(h5_dir)){
      tryCatch(HDF5Array::saveHDF5SummarizedExperiment(x = se, dir = h5_dir, replace = TRUE),
               error = function(e) message("The dataset is not saved. Please save manually, using the HDF5Array::saveSummarizedExperiment command. "))
    }
  }else{
    se = SummarizedExperiment::SummarizedExperiment(assays = list(beta = beta_mat, cov = cov_mat),
                                                    metadata = list(genome = genome_name, is_h5 = is_hdf5, summary = se_summary,
                                                                    chr_summary = chr_summary, ref_CpG = ref_cpg_dt,
                                                                    chrom_sizes = chrom_sizes, descriptive_stats = desc),
                                                    colData = col_data, rowData = cpg_loci)
  }

  return(methrix(se))
}
