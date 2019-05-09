#' Creates a detailed interative html summary report from Methrix object
#' @description Creates a detailed interative html summary report from Methrix object.
#' If the directory contains required files (from previous run), it directly proceeds to generate html report.
#' @param meth \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved. If \code{NULL} creats a \code{tempdir}
#' @param plot_beta_dist Default FALSE. This is a time consuming and writes huge density files required for plotting.
#' @param n_thr Default 4. Only used if \code{plot_beta_dist} is TRUE
#' @return an interactive html report
#' @examples
#' data("mm9_bsmap")
#' methrix::methrix_report(meth = mm9_bsmap)
#' @export
methrix_report = function(meth, output_dir = NULL, plot_beta_dist = FALSE, n_thr = 4){

  start_proc_time = proc.time()

  if(is.null(output_dir)){
    output_dir = tempdir()
  }else{
    if(!dir.exists(output_dir)){
      dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
    }
  }

  m_chr_summary = methrix::get_matrix(m = meth, type = "M", add_loci = TRUE)
  c_chr_summary = methrix::get_matrix(m = meth, type = "C", add_loci = TRUE)

  #summarize matrix per chr
  ##Methylation
  cat(paste0("Step 1 of 7\n"))
  of1 = normalizePath(paste0(output_dir, "/m_per_chr.tsv"))
  if(file.exists(of1)){
    message("File already present. Skipping step 1..")
  }else{
    mean_meth_chr = m_chr_summary[,c(1, 4:ncol(m_chr_summary)), with = FALSE][,lapply(.SD, mean, na.rm = TRUE), chr]
    median_meth_chr = m_chr_summary[,c(1, 4:ncol(m_chr_summary)), with = FALSE][,lapply(.SD, median, na.rm = TRUE), chr]
    meth_per_chr_dt = data.table::rbindlist(l = list(mean = mean_meth_chr, median = median_meth_chr), idcol = "stat", use.names = TRUE)
    data.table::fwrite(x = meth_per_chr_dt, file = of1, sep = "\t")
  }

  ##Coverage
  cat(paste0("Step 2 of 7\n"))
  of2 = normalizePath(paste0(output_dir, "/c_per_chr.tsv"))
  if(file.exists(of2)){
    message("File already present. Skipping step 2..")
  }else{
    mean_cov_chr = c_chr_summary[,c(1, 4:ncol(c_chr_summary)), with = FALSE][,lapply(.SD, function(x){mean(x[!x == 0])}), chr]
    median_cov_chr = c_chr_summary[,c(1, 4:ncol(c_chr_summary)), with = FALSE][,lapply(.SD, function(x){median(x[!x == 0])}), chr]
    cov_per_chr_dt = data.table::rbindlist(l = list(mean = mean_cov_chr, median = median_cov_chr), idcol = "stat", use.names = TRUE)
    data.table::fwrite(x = cov_per_chr_dt, file = of2, sep = "\t")
  }


  #n CpGs covered per chromomse
  cat(paste0("Step 3 of 7\n"))
  of3 = normalizePath(paste0(output_dir, "/n_covered_per_chr.tsv"))
  contig_nCpGs = meth@metadata$ref_CpG
  colnames(contig_nCpGs) = c("chr", "total_CpGs")
  if(file.exists(of3)){
    message("File already present. Skipping step 3..")
  }else{
    n_covered_chr = m_chr_summary[,c(1, 4:ncol(m_chr_summary)), with = FALSE][,lapply(.SD, function(x){length(x[!is.na(x)])}), chr]
    #contig_nCpGs = meth@metadata$ref_CpG
    #colnames(contig_nCpGs) = c("chr", "total_CpGs")
    n_covered_chr = merge(n_covered_chr, contig_nCpGs, by = 'chr', all.x = TRUE)
    data.table::fwrite(x = n_covered_chr, file = of3, sep = "\t")
  }


  #Common CpGs covered by all samples
  cat(paste0("Step 4 of 7\n"))
  of4 = normalizePath(paste0(output_dir, "/n_covered_by_all_samples.tsv"))
  if(file.exists(of4)){
    message("File already present. Skipping step 4..")
  }else{
    mf = methrix::coverage_filter(m = meth, cov_thr = 1, min_samples = nrow(colData(meth)), n_threads = n_thr)
    mf_chr_summary = methrix::getChrSummary(x = mf)
    colnames(mf_chr_summary) = c("chr", "n_CpG")
    rm(mf)
    gc(verbose = FALSE)
    mf_chr_summary = merge(mf_chr_summary, contig_nCpGs, by = 'chr', all.x = TRUE)
    mf_chr_summary[, fract_CpG := n_CpG/total_CpGs]
    data.table::fwrite(x = mf_chr_summary, file = of4, sep = "\t")
  }

  #Global methylation
  cat(paste0("Step 5 of 7\n"))
  of5 = normalizePath(paste0(output_dir, "/global_meth_per_samp.tsv"))
  if(file.exists(of5)){
    message("File already present. Skipping step 5..")
  }else{
    mean_meth = apply(m_chr_summary[,4:ncol(m_chr_summary)], 2, mean, na.rm = TRUE)
    median_meth = apply(m_chr_summary[,4:ncol(m_chr_summary)], 2, median, na.rm = TRUE)
    global_meth = data.table::data.table(Sample_Name = names(mean_meth), mean = mean_meth, median = median_meth)
    data.table::fwrite(x = global_meth, file = of5, sep = "\t")
  }


  #Global methylation
  cat(paste0("Step 6 of 7\n"))
  of6 = normalizePath(path = paste0(output_dir, "/global_cov_per_samp.tsv"))
  if(file.exists(of6)){
    message("File already present. Skipping step 6..")
  }else{
    mean_cov = apply(c_chr_summary[,4:ncol(c_chr_summary)], 2, function(x){mean(x[!x == 0])})
    median_cov = apply(c_chr_summary[,4:ncol(c_chr_summary)], 2, function(x){median(x[!x == 0])})
    global_cov = data.table::data.table(Sample_Name = names(mean_cov), mean = mean_cov, median = median_cov)
    data.table::fwrite(x = global_cov, file = of6, sep = "\t")
  }

  #Density plot data
  cat(paste0("Step 7 of 7\n"))
  if(plot_beta_dist){
    dens_files = list.files(path = output_dir, pattern = "*_density\\.tsv$")
    if(length(dens_files) == nrow(colData(meth))){
      message("Files already present. Skipping step 7..")
    }else{
      parallel::mclapply(X = seq_len(nrow(colData(meth))), FUN = function(i){
        col_idx = i+3
        i_dens = density(unlist(m_chr_summary[,col_idx, with = FALSE])[1:10000], na.rm = TRUE)
        data.table::fwrite(x = data.table::data.table(x = i_dens$x, y = i_dens$y),
                           file = paste0(output_dir, "/", rownames(colData(x = meth))[i], "_density.tsv"), sep = "\t")

      }, mc.cores = n_thr)
    }
  }

  cat(paste0("Knitting report\n"))
  md = system.file('extdata', 'summarize_methrix.Rmd', package = 'methrix')

  rmarkdown::render(input = md, output_file = "methrix_reports.html",
                    output_dir = output_dir, clean = TRUE, params = list(n_covered_tsv = of3,
                                                                         n_covered_by_all_samples_tsv = of4,
                                                                         m_stat_tsv = of1,
                                                                         c_stat_tsv = of2,
                                                                         global_meth_per_samp_tsv = of5,
                                                                         global_cov_per_samp_tsv = of6))
  browseURL(url = paste0(output_dir, "/methrix_reports.html"))

  cat(data.table:::timetaken(started.at = start_proc_time), sep = "\n")
}
