# methrix
An R S4 class for the efficient storage and manipulation of DNA methylation data 

### Summary:
* Faster summarization of Bismark output with `data.table` backend
* Pre-defined loci
* Supports serialized arrays with `HDF5Array` and `saveHDF5SummarizedExperiment`
* Vectorized code (faster, memory expensive) and non-vectorized code (slower, minimal memory)

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

- [] Write methods
- [] Full implemetation of data.table backend
- [] Avoid matrices
- [] new class independent of summarized experiment
- [] benchmark memory consumption

