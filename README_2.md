# methrix DSS Support

DSS-equivalent DML/DMR testing directly on methrix objects, with no BSseq conversion. A single API (`methrix_DMLtest`) auto-dispatches based on whether the object is in-memory or on-disk.

## Files

```
R/
├── methrix_DMLtest.R              # Public API + dispatcher
├── methrix_DMLtest_inmemory.R     # In-memory path (uses DSS internals)
├── methrix_DMLtest_ondisk.R       # On-disk path (full DelayedArray reimplementation)
└── methrix_callDML_callDMR.R      # callDML/callDMR (shared, no DSS dependency)
```

## How It Works

```
methrix_DMLtest(m, group1, group2, ...)
       │
       ├── is_h5(m) == FALSE  →  .DMLtest_inmemory()
       │     • Assays are already plain matrices
       │     • Extracts x, n, chr, pos directly (no BSseq)
       │     • Calls DSS:::compute.mean.*, est.dispersion.BSseq, waldTest.DML
       │     • Uses DSS's C code + mclapply (fast, exact DSS statistics)
       │     • DSS is a dependency
       │
       └── is_h5(m) == TRUE   →  .DMLtest_ondisk()
             • Assays are HDF5-backed DelayedMatrices
             • Column subset, NA→0, M-count computation all stay DELAYED
             • Mean estimation: DelayedArray::rowSums (auto-blocked)
             • Smoothing: per-chromosome chunk realization + R sliding window
             • Dispersion: block-wise MOM (500K rows), global Bayesian shrinkage
             • Wald test: block-wise variance computation
             • No DSS dependency, no full materialization
             • Peak memory: ~50-100 MB for any genome size
```

## Usage

```r
library(methrix)
# source all four files, or add to package R/ directory

# Works identically on both types:
result <- methrix_DMLtest(meth,
  group1 = c("C1", "C2", "C3"),
  group2 = c("N1", "N2", "N3"),
  smoothing = TRUE)

dmls <- methrix_callDML(result, p.threshold = 0.001)
dmrs <- methrix_callDMR(result, p.threshold = 0.001)
```

## Dependencies

| Path | Required packages |
|---|---|
| In-memory | `methrix`, `DSS`, `SummarizedExperiment`, `GenomicRanges` |
| On-disk | `methrix`, `DelayedArray`, `DelayedMatrixStats`, `SummarizedExperiment`, `GenomicRanges`, `BiocParallel` |
| callDML/callDMR | None (base R only) |

## Memory Comparison (28M CpGs × 20 samples, comparing 3 vs 3)

| Pathway | Peak RAM | BSseq created | Materialization |
|---|---|---|---|
| Standard (methrix2bsseq → DSS) | ~18 GB | Yes | Full genome, all 20 samples, twice |
| In-memory path (this) | ~3 GB | No | Only group columns (6 of 20), once |
| On-disk path (this) | ~50-100 MB | No | 500K-row blocks, one chromosome at a time |

## Integration

Add all four `.R` files to the methrix package `R/` directory. Add to NAMESPACE:
```
export(methrix_DMLtest)
export(methrix_callDML)
export(methrix_callDMR)
```

Add to DESCRIPTION Suggests: `DSS` (needed only for in-memory path).
