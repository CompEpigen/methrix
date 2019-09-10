## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tinytex.verbose = TRUE)

## ---- message=FALSE, warning=FALSE---------------------------------------
#Load library
library(methrix)
#Genome of your preference to work with
library(BSgenome.Hsapiens.UCSC.hg19)

## ------------------------------------------------------------------------
#Example bedgraph files
bdg_files = list.files(
  path = system.file('extdata', package = 'methrix'),
  pattern = "*bdg\\.gz$",
  full.names = TRUE
)

print(bdg_files)

#Generate some sample annotation table
sample_anno = data.frame(
  row.names = gsub(
    pattern = "\\.bdg\\.gz$",
    replacement = "",
    x = basename(bdg_files)),
  Condition = c("Cancer", 'Cancer', "Normal", "Normal"),
  stringsAsFactors = FALSE)

print(sample_anno)

## ---- warning=FALSE------------------------------------------------------
#First extract genome wide CpGs from the desired genome
hg19_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

## ------------------------------------------------------------------------
#Read the files
meth = methrix::read_bedgraphs(
  files = bdg_files,
  ref_cpgs = hg19_cpgs,
  chr_idx = 1,
  start_idx = 2,
  M_idx = 3,
  U_idx = 4,
  stranded = TRUE,
  collapse_strands = TRUE,
  coldata = sample_anno
)

## ------------------------------------------------------------------------
#Typing meth shows basic summary.
meth

## ---- eval=FALSE---------------------------------------------------------
#  methrix::methrix_report(meth = meth, output_dir = tempdir())

## ------------------------------------------------------------------------
meth = methrix::remove_uncovered(m = meth)
meth

## ------------------------------------------------------------------------
#Coverage matrix
coverage_mat = methrix::get_matrix(m = meth, type = "C")
head(coverage_mat)

## ------------------------------------------------------------------------
#Methylation matrix
meth_mat = methrix::get_matrix(m = meth, type = "M")
head(meth_mat)

## ------------------------------------------------------------------------
#If you prefer you can attach loci info to the matrix
meth_mat_with_loci = methrix::get_matrix(m = meth, type = "M", add_loci = TRUE)
meth_mat_with_loci

## ------------------------------------------------------------------------
#e.g; Retain all loci which are covered at-least in two sample by 3 or more reads
methrix::coverage_filter(m = meth, cov_thr = 3, min_samples = 2)

## ------------------------------------------------------------------------
#Retain sites only from chromosme chr21
methrix::subset_methrix(m = meth, contigs = "chr21")

## ------------------------------------------------------------------------
#e.g; Retain sites only in TP53 loci
target_loci = data.table::data.table(chr = "chr21", start = 27867971, end =  27868103)
print(target_loci)

methrix::subset_methrix(m = meth, regions = target_loci)

## ------------------------------------------------------------------------
methrix::subset_methrix(m = meth, samples = "C1")

#Or you could use [] operator to subset by index
meth[,1]

## ------------------------------------------------------------------------
meth_stats = get_stats(m = meth)
print(meth_stats)

## ------------------------------------------------------------------------
#Draw mean coverage per sample
plot_stats(plot_dat = meth_stats, what = "C", stat = "mean")
#Draw mean methylation per sample
plot_stats(plot_dat = meth_stats, what = "M", stat = "mean")

## ------------------------------------------------------------------------
mpca = methrix_pca(m = meth, do_plot = FALSE)

#Plot PCA results
plot_pca(pca_res = mpca, show_labels = TRUE)

#Color code by an annotation
plot_pca(pca_res = mpca, m = meth, col_anno = "Condition")

## ------------------------------------------------------------------------
#Violin plots
methrix::plot_violin(m = meth)

## ------------------------------------------------------------------------
methrix::plot_coverage(m = meth,type = "dens")

## ------------------------------------------------------------------------
library(bsseq)

## ------------------------------------------------------------------------
bs_seq = methrix::methrix2bsseq(m = meth)

bs_seq

## ------------------------------------------------------------------------
sessionInfo()

