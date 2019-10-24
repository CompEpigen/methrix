#   Source  and generation of example data
#   Source:
#   WGBS for colon cancer, chr21 and chr22
#
# This is a subset of original `bsseqData` converted to `methrix` containing Whole-genome bisulfite sequencing data (WGBS)
# for colon cancer on chromosome 21 and 22.
#
# references: Hansen, K. D. et al. (2011) Increased methylation variation in epigenetic domains across cancer types.
# Nature Genetics 43, 768-775.
#

library(methrix)

if(!requireNamespace("rtracklayer")) {
  BiocManager::install("rtracklayer")
}
library(rtracklayer)

if(!requireNamespace("DSS")) {
  BiocManager::install("DSS")
}
library(DSS)

if(!requireNamespace("bsseqData")) {
  BiocManager::install("bsseqData")
}
library(bsseqData)

if(!requireNamespace("bsseq")) {
  BiocManager::install("bsseq")
}
library(bsseq)

if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)

#get bsseq data
data(BS.cancer.ex)
BS.cancer.ex

BS.cancer.ex <- BS.cancer.ex[, c(1,2,4,5)]

#extract data
meth <- bsseq::getMeth(BS.cancer.ex, type="raw")
cov <- bsseq::getCoverage(BS.cancer.ex)
gr <- granges(BS.cancer.ex)
gr <- as.data.frame(gr)[,1:2]

#calculate M values
count_meth <- meth*cov
count_unmethylated <- cov - count_meth

df <- list()
location <- "../methrix_data_generation/"
dir.create(location, showWarnings = FALSE, recursive = TRUE)

#create bedgraphs and save them
for (i in colnames(BS.cancer.ex)){
  df[[i]] <- data.frame(chromosome=gr$seqnames, position=gr$start, count_methylated=count_meth[,i],
                        count_unmethylated=count_unmethylated[,i] )
  df[[i]] <- df[[i]][!is.na(df[[i]]$count_methylated),]
  write.table(x = df[[i]], file = paste0(location, i, ".bedGraph"), quote = FALSE, row.names = FALSE, sep = "\t")
  }


#extract pheno
pheno <- pData(BS.cancer.ex)

#create methrix object
hg19_cpg <-  extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

files <- dir(path = location, full.names = TRUE,
             pattern = ".bedGraph$")
colnames(pData(BS.cancer.ex)) <- c("Condition", "Pair")

methrix_obj <- methrix::read_bedgraphs(files = files,coldata = pData(BS.cancer.ex),
                        zero_based = FALSE,
                        stranded = FALSE,
                        ref_build = "hg19",
                        ref_cpgs = hg19_cpg,
                        chr_idx = 1,
                        start_idx = 2,
                        M_idx = 3,
                        U_idx = 4)


#run DSS to define DMRs
dmlTest = DSS::DMLtest(BS.cancer.ex, group1=c("C1", "C2"), group2=c("N1", "N2"))
#dmls2 = callDML(dmlTest, delta=0.1, p.threshold=0.001)
dmrs = DSS::callDMR(dmlTest, delta = 0.1, p.threshold = 0.001)

#extend dmrs
dmrs_gr <- bsseq::data.frame2GRanges(dmrs)
start(dmrs_gr) <- start(dmrs_gr) - 2000
end(dmrs_gr) <- end(dmrs_gr) + 2000

methrix_dta <- subset_methrix(methrix_obj, regions = dmrs_gr)
methrix_data <- remove_uncovered(methrix_dta)
usethis::use_data(methrix_data, overwrite = TRUE)

meth <- get_matrix(methrix_data, type="M")
cov <- get_matrix(methrix_data, type="C")
gr <- methrix_data@elementMetadata


#create M and U values
count_meth_sub <- cov*meth
count_unmethylated_sub <- cov - count_meth_sub
df <- list()
#save bedgraphs
for (i in colnames(count_meth_sub)){
  df[[i]] <- data.frame(chromosome=gr$chr, position=gr$start, count_methylated=count_meth_sub[,i],
                        count_unmethylated=count_unmethylated_sub[,i] )
  df[[i]] <- df[[i]][!is.na(df[[i]]$count_methylated),]
  data.table::fwrite(df[[i]], file = paste0("./inst/extdata/",i,".bedGraph.gz"),
                     compress = "auto", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

