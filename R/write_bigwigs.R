#' Exports methrix object as bigWigs
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved. Default getwd()
#' @param samp_names sample names to export
#' @examples
#' data('methrix_data')
#' write_bigwigs(m = methrix_data, output_dir = './temp')
#' @return NULL
#' @importFrom rtracklayer export
#' @export

write_bigwigs = function(m, output_dir = getwd(), samp_names = NULL){
  
  if (!dir.exists(output_dir)) {
    dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  mat_gr <- methrix::get_matrix(m = m, type = "M", add_loci = TRUE, in_granges = TRUE)
  
  seql = m@metadata$chrom_sizes$length
  names(seql) = m@metadata$chrom_sizes$contig
  
  all_samps = names(mcols(mat_gr))  
  
  if(is.null(samp_names)){
    samp_names = all_samps
  }else{
    samp_names = intersect(samp_names, all_samps)
    if(length(samp_names) == 0){
      stop("Incorrect sample names!")
    }
  }
  
  message("----------------------")
  for(samp in samp_names){
    op_bw = paste0(output_dir, "/", samp, ".bw")
    message("*Writing ", op_bw)
    samp_gr = mat_gr[,samp]
    names(mcols(samp_gr)) = "score"
    samp_gr = samp_gr[!is.na(samp_gr$score)]
    seqlengths(samp_gr) = seql[names(seqlengths(samp_gr))]
    rtracklayer::export(samp_gr, con = paste0(output_dir, "/", samp, ".bw"), format="bigWig")
  }
  message("----------------------")
}