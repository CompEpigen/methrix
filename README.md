## Methrix
An R S4 class for the efficient storage and manipulation of DNA methylation data 

### Summary:
* Faster summarization of Bismark output with `data.table` backend
* Pre-defined loci
* Supports serialized arrays with `HDF5Array` and `saveHDF5SummarizedExperiment`
* Vectorized code (faster, memory expensive) and non-vectorized code (slower, minimal memory)

### Usage
#### Extract CpG loci from the reference genome
```r
> hs37d5_cpgs = methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5")
# Extracting CpGs..
# Here is a Chuck Norris joke while you wait..
# Chuck Norris doesn't chew gum. Chuck Norris chews tin foil.
# Done.
# Extracted 28787054 from 86 contigs.
```
#### Read in bedgraphs and generate a methrix object (which inherits summarizedExperiments)
vect_batch_size = 3 -> Prcess 3 files in a batch

Test data is 3 stranded WGBS bedgraphs from MethylcTools

```r
> x0 = methrix::read_bedgraphs(files = bdg_files[1:3], pipeline = "MethylcTools", collapse_starnds = TRUE, vect_batch_size = 3, ref_build = "Hs37d5", ref_cpgs = hs37d5_cpgs, verbose = FALSE)
# Processing batch 1 of 1
# |--------------------------------------------------|
# |==================================================|
# Missing 0 from: tumor00_JMMLC_D117.CG.bed.gz
# |--------------------------------------------------|
# |==================================================|
# Missing 0 from: tumor00_JMMLC_D123.CG.bed.gz
# |--------------------------------------------------|
# |==================================================|
# Missing 0 from: tumor00_JMMLC_D124.CG.bed.gz
# 00:01:13 elapsed (00:02:58 cpu)

> x0
# An object of class  methrix 
#                 ID  Summary
# 1:       n_samples        3
# 2:          n_CpGs 28162972
# 3: Reference_Build   Hs37d5
# 4:           is_H5    FALSE
```

#### Filter matrices based on coverage statistics (e.g: at-least two samples should have a loci covered by 4 reads)
```r
> x5 = methrix::filter_methrix(m = x0, cov_thr = 5, min_samples = 2, n_threads = 4)
# Retained 8916847 of 28162972 sites

#New filtered methrix object
>  x5
# An object of class  methrix 
#                 ID Summary
# 1:       n_samples       3
# 2:          n_CpGs 8916847
# 3: Reference_Build  Hs37d5
# 4:           is_H5   FALSE 
```

### Quick (and dirty) benchmark

Data: 5 single cell WGBS coverage files form bismark

```r
bench_mark = microbenchmark::microbenchmark(methrix_vectorized = methrix::read_bismark(files = bm_cov[1:5], vectorize = TRUE), methrix = methrix::read_bismark(files = bm_cov[1:5], vectorize = FALSE), bsseq = bsseq::read.bismark(files = bm_cov[1:5]), times = 10)

> bench_mark
Unit: seconds
               expr      min       lq     mean   median        uq       max neval
 methrix_vectorized 8.542675 9.314229 9.822698 9.498361 10.269085 11.821858    10
            methrix 7.283596 7.597634 7.990360 7.672935  8.128339 10.047221    10
              bsseq 2.620453 2.947310 3.808507 3.729243  4.742076  5.786401    10
```

### To do

- [X] Write methods
- [] Full implemetation of data.table backend
- [] Avoid matrices
- [] new class independent of summarized experiment
- [] benchmark memory consumption

