## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
#Load library
library(methrix)
library(BSgenome.Hsapiens.1000genomes.hs37d5) #Genome of your preference to work with

## ------------------------------------------------------------------------
#Get the list of files to process
bed_files = list.files(path = "~/C010-Datasets/Internal/2019-03-25_JMMLC_PBAT/JMML_Methylation/01_raw_data/bedgraphs/CG/", 
                       pattern = "*.bed.gz", full.names = TRUE)
#Normalize the path
bed_files = normalizePath(path = bed_files)
#For vignette purpose use three files
bed_files = bed_files[1:3]

print(bed_files)

## ---- warning=FALSE------------------------------------------------------
#First extract genome wide CpGs from the desired genome
hs37d5_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5", bored = FALSE)

## ------------------------------------------------------------------------
#Read the files
meth = methrix::read_bedgraphs(files = bed_files, 
                               pipeline = "MethylcTools",
                               collapse_starnds = TRUE, 
                               ref_build = "Hs37d5", 
                               vect_batch_size = 3, 
                               ref_cpgs = hs37d5_cpgs)

## ------------------------------------------------------------------------
#Typing meth shows basic summary.
meth

## ------------------------------------------------------------------------
#Get chromosome wise number of CpGs
getChrSummary(x = meth)

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
meth_mat_with_loci = methrix::get_matrix(m = meth, type = "M", add_lcoi = TRUE)
meth_mat_with_loci

## ------------------------------------------------------------------------
meth_clean = methrix::remove_uncovered(m = meth)
meth_clean

## ------------------------------------------------------------------------
#e.g; Retain all loci which are covered at-least in two sample by 3 or more reads
meth_cov_filtered = methrix::coverage_filter(m = meth, cov_thr = 3, min_samples = 2, n_threads = 6)
meth_cov_filtered

## ------------------------------------------------------------------------
#Retain sites only from chromosme 1 and "X"
meth_1x = methrix::subset_methrix(m = meth, contigs = c(1, "X"))
meth_1x

## ------------------------------------------------------------------------
#e.g; Retain sites only in TP53 loci 
tp53_loci = data.table::data.table(chr = "17", start = "7571720", end = "7590868")
meth_tp53 = methrix::subset_methrix(m = meth, regions = tp53_loci)
meth_tp53

## ------------------------------------------------------------------------
meth_d117 = methrix::subset_methrix(m = meth, samples = "tumor00_JMMLC_D117")
meth_d117

## ------------------------------------------------------------------------
meth_sd = methrix::order_by_sd(m = meth)
head(methrix::get_matrix(m = meth_sd, type = "M"))

## ------------------------------------------------------------------------
library(bsseq)

## ------------------------------------------------------------------------
bs_seq = methrix::methrix2bsseq(m = meth)

bs_seq

## ------------------------------------------------------------------------
sessionInfo()

