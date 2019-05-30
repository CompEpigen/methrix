# methrix 0.9.90

* Added plotting functions
* Added DMR calling function
* Support for HDF5 backend

# methrix 0.9.86

* `methrixreport` works on Windows
* Error handling when bedgraph file contains non-reference CpGs
* Subsetting operation uses SE functions
* Code-cleanup, improve log messages
* Added example files for vignette creation
* Update description file
* Added an example dataset
* Added examples for all functions
* Cleared BiocChecks

* Bug fix: Fix NA's caused by `starnd_collapse` from bedgraph files with no M and U info

# methrix 0.9.85

Significant user level improvements:
* `stranded` and `collapse_strands` now works with bedgraphs without strand info (such as from MethylDackel)
* Added `methrix_report` function: one click function to generate comprehensive interactive html report from methrix object
* Added `write_bedgraphs` function: [Still needs some fixes from data.table]

* Bug fix: Replaced NaN's caused by 0/0 to NA's
* Bug fix: `methrix2bsseq` error fix due to NA's
* Add `Tutorial` vignette
* Bug fix: Fixed M values passed to BSseq object creation

# methrix 0.9.8 [2019-03-28]

* Updated `Methrix` object
* Added `getChrSummary` method
* Added `subset_methrix` function for various subset operations
* Added `get_matrix` function for extracting methylation or coverage matrices
* Added `methrix2bsseq` function for easy converting methrix to bsseq object

# methrix 0.9.7 [2019-03-26]

* Added `MethylcTools` to supposrted pipelines
* Added `collapse_starnds` option to summarize CpG sites by strand
* Added `create_methrix` function for object creation
* Updated `methrix` object
* Added `order_by_sd` and `filter_methrix` functions
* Code clean-up.

# methrix 0.9.5 [2019-03-14]

* Added annoying Chuck-Norris jokes while extracting CpGs
* Added `extract_CPGs` function to extract from BSGenome object
* HDF5 matrix support
* Updated code for non-vectorized code
