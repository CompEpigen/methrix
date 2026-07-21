## Correctness tests for methrix_DMLtest() (dispatcher over in-memory / on-disk).
##
## Requires DSS + bsseq; skipped automatically if unavailable so the rest of the
## methrix suite still runs. Assumes both R/methrix_dss.R and
## R/methrix_DMLtest_inmemory.R are loaded.

skip_if_no_dss <- function() {
    testthat::skip_if_not_installed("DSS")
    testthat::skip_if_not_installed("bsseq")
}

example_groups <- function(m) {
    sn <- rownames(SummarizedExperiment::colData(m))
    n  <- length(sn)
    list(g1 = sn[seq_len(floor(n / 2))],
         g2 = sn[(floor(n / 2) + 1L):n])
}

## align two DMLtest frames on shared (chr,pos) and return paired finite stats
paired_stats <- function(a, b) {
    ka <- paste(a$chr, a$pos); kb <- paste(b$chr, b$pos)
    common <- intersect(ka, kb)
    sa <- a$stat[match(common, ka)]; sb <- b$stat[match(common, kb)]
    ok <- is.finite(sa) & is.finite(sb)
    list(n = length(common), a = sa[ok], b = sb[ok])
}

test_that("in-memory path reproduces stock DSS::DMLtest (via bsseq)", {
    skip_if_no_dss()
    data("methrix_data", package = "methrix")
    m <- methrix_data
    g <- example_groups(m)

    ref  <- DSS::DMLtest(methrix2bsseq(m), group1 = g$g1, group2 = g$g2,
                         equal.disp = TRUE, smoothing = FALSE)
    ours <- methrix_DMLtest(m, g$g1, g$g2, equal.disp = TRUE,
                            smoothing = FALSE, verbose = FALSE)  # -> .DMLtest_inmemory

    ps <- paired_stats(ref, ours)
    expect_gt(ps$n, 100)
    ## genome-wide dispersion in both -> effectively identical
    expect_gt(stats::cor(ps$a, ps$b), 0.9999)
    expect_lt(stats::median(abs(ps$a - ps$b)), 1e-3)
})

test_that("on-disk path matches the in-memory path within chunking tolerance", {
    skip_if_no_dss()
    data("methrix_data", package = "methrix")
    m_mem <- methrix_data
    m_h5  <- convert_methrix(m_mem)
    g     <- example_groups(m_mem)

    r_mem <- methrix_DMLtest(m_mem, g$g1, g$g2, equal.disp = TRUE,
                             smoothing = FALSE, verbose = FALSE)
    r_h5  <- methrix_DMLtest(m_h5,  g$g1, g$g2, equal.disp = TRUE,
                             smoothing = FALSE, verbose = FALSE)

    ## same sites are tested (coverage filter is identical)
    expect_identical(paste(r_mem$chr, r_mem$pos), paste(r_h5$chr, r_h5$pos))

    ## per-chromosome dispersion prior on-disk -> near-exact, not bit-identical.
    ## Threshold is looser than the other cor() checks in this file: this is the
    ## only comparison between a single genome-wide dispersion prior (in-memory)
    ## and independent per-chromosome priors (on-disk), which on this bundled
    ## low-coverage 2-chromosome dataset measures ~0.9984 (observed on chr22,
    ## the smaller/lower-coverage chromosome).
    ps <- paired_stats(r_mem, r_h5)
    expect_gt(stats::cor(ps$a, ps$b), 0.995)
    expect_lt(stats::median(abs(ps$a - ps$b)), 0.05)
})

test_that("on-disk result is a valid DMLtest object consumable by callDMR", {
    skip_if_no_dss()
    data("methrix_data", package = "methrix")
    m <- convert_methrix(methrix_data)
    g <- example_groups(methrix_data)

    res <- methrix_DMLtest(m, g$g1, g$g2, equal.disp = TRUE,
                           smoothing = TRUE, verbose = FALSE)

    expect_true("DMLtest" %in% class(res))
    expect_identical(colnames(res),
                     c("chr", "pos", "mu1", "mu2", "diff", "diff.se",
                       "stat", "phi1", "phi2", "pval", "fdr"))
    expect_false(is.unsorted(res$pos[res$chr == res$chr[1]]))  # sorted within chr
    expect_silent(DSS::callDMR(res, p.threshold = 0.05))
})

test_that("sub-chromosome blocking with halo matches whole-chromosome blocking", {
    skip_if_no_dss()
    data("methrix_data", package = "methrix")
    m <- convert_methrix(methrix_data)
    g <- example_groups(methrix_data)

    whole <- methrix_DMLtest(m, g$g1, g$g2, equal.disp = TRUE, smoothing = TRUE,
                             block_size = NULL,  verbose = FALSE)
    split <- methrix_DMLtest(m, g$g1, g$g2, equal.disp = TRUE, smoothing = TRUE,
                             block_size = 2000L, verbose = FALSE)

    ps <- paired_stats(whole, split)
    expect_gt(stats::cor(ps$a, ps$b), 0.999)   # halo keeps edge smoothing correct
})
