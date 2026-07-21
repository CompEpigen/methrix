################################################################################
##
##  methrix_dss.R
##
##  DSS differential-methylation testing driven directly off a methrix object,
##  with NO conversion to bsseq.
##
##  This file provides:
##    * methrix_DMLtest()   - public dispatcher on is_h5(m)
##    * .DMLtest_ondisk()   - on-disk (HDF5) backend: streams one chromosome at
##                            a time so the full genome-wide matrix is never
##                            realised in memory.
##
##  The in-memory backend .DMLtest_inmemory() lives in
##  methrix_DMLtest_inmemory.R and is used unchanged for in-memory objects
##  (single genome-wide dispersion estimate == exact stock-DSS behaviour).
##
##  DSS is a required dependency. The heavy, validated statistics are reused
##  from DSS's non-exported functions (via asNamespace("DSS")) rather than
##  reimplemented: pin the DSS version in DESCRIPTION and rely on
##  tests/testthat/test-methrix_dss.R. Functions used: hasCoverage,
##  compute.mean.noSmooth, compute.mean.Smooth, est.dispersion.BSseq,
##  compute.waldStat.noSmooth, compute.waldStat.Smooth.
##
################################################################################

#' Resolve group labels (names or indices) to sample names of a methrix object
#' @keywords internal
#' @noRd
.methrix_resolve_groups <- function(m, group1, group2) {
    sn <- colnames(m)
    to_names <- function(g, id) {
        if (is.numeric(g)) {
            if (min(g) < 1L || max(g) > length(sn))
                stop("Invalid sample indices in group", id, ".")
            g <- sn[g]
        }
        if (!all(g %in% sn))
            stop("Some sample names in group", id, " are not in the methrix object.")
        g
    }
    g1 <- to_names(group1, 1)
    g2 <- to_names(group2, 2)
    if (length(intersect(g1, g2)) > 0)
        stop("group1 and group2 share samples.")
    list(g1 = g1, g2 = g2)
}

#' Build a per-chromosome (optionally sub-chromosome + halo) block grid over rows
#' @keywords internal
#' @noRd
.methrix_dss_blocks <- function(chr, pos, smoothing, smoothing.span, block_size) {
    ## methrix rows are ordered by chromosome then position; take contiguous runs
    ends   <- cumsum(S4Vectors::runLength(S4Vectors::Rle(chr)))
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    ## halo must cover both the smoothing window and the 200 bp hasCoverage window
    halo   <- if (smoothing) max(as.integer(smoothing.span), 200L) else 200L

    blocks <- list()
    for (k in seq_along(starts)) {
        rng <- starts[k]:ends[k]
        if (is.null(block_size) || length(rng) <= block_size) {
            blocks[[length(blocks) + 1L]] <- list(idx = rng, keep = rep(TRUE, length(rng)))
        } else {
            for (c0 in seq(1L, length(rng), by = block_size)) {
                core <- rng[c0:min(c0 + block_size - 1L, length(rng))]
                ext  <- rng[pos[rng] >= pos[core[1L]] - halo &
                            pos[rng] <= pos[core[length(core)]] + halo]
                blocks[[length(blocks) + 1L]] <- list(idx = ext, keep = ext %in% core)
            }
        }
    }
    blocks
}

#' On-disk (HDF5) backend for methrix_DMLtest: stream one chromosome per block
#' @keywords internal
#' @noRd
#' @importFrom stats p.adjust
.DMLtest_ondisk <- function(m, group1, group2, equal.disp, smoothing,
                            smoothing.span, n_cores, block_size, verbose) {

    if (is.null(n_cores)) n_cores <- 1L
    grp <- .methrix_resolve_groups(m, group1, group2)
    g1  <- match(grp$g1, colnames(m))
    g2  <- match(grp$g2, colnames(m))

    rd  <- SummarizedExperiment::rowData(m)
    chr <- as.character(rd$chr)
    pos <- as.integer(rd$start)
    chr_levels <- unique(chr)

    blocks <- .methrix_dss_blocks(chr, pos, smoothing, smoothing.span, block_size)

    ## DSS internals (fetched once so mclapply children inherit them)
    DSS_ns    <- asNamespace("DSS")
    hasCov    <- get("hasCoverage",             envir = DSS_ns)
    mean_ns   <- get("compute.mean.noSmooth",   envir = DSS_ns)
    mean_sm   <- get("compute.mean.Smooth",     envir = DSS_ns)
    est_disp  <- get("est.dispersion.BSseq",    envir = DSS_ns)
    wald_ns   <- get("compute.waldStat.noSmooth", envir = DSS_ns)
    wald_sm   <- get("compute.waldStat.Smooth",   envir = DSS_ns)
    ## keep DSS's own dispersion parallelism serial; we parallelise over blocks
    serial    <- 1L

    worker <- function(b) {
        mb  <- m[b$idx, ]                                        # lazy subset for H5
        cov <- as.matrix(methrix::get_matrix(mb, type = "C"))    # one contiguous read
        cov[is.na(cov)] <- 0L
        beta <- as.matrix(methrix::get_matrix(mb, type = "M"))
        beta[is.na(beta)] <- 0
        M    <- round(beta * cov)                                # integer methylated counts
        rm(beta)

        ac <- chr[b$idx]; ap <- pos[b$idx]
        x1 <- M[,   g1, drop = FALSE]; n1 <- cov[, g1, drop = FALSE]
        x2 <- M[,   g2, drop = FALSE]; n2 <- cov[, g2, drop = FALSE]
        rm(M, cov)
        storage.mode(x1) <- "integer"; storage.mode(n1) <- "integer"
        storage.mode(x2) <- "integer"; storage.mode(n2) <- "integer"

        ## DSS low-coverage prefilter (200 bp window), per group, as in DMLtest
        ix <- hasCov(n1, ap) & hasCov(n2, ap)
        if (!any(ix)) return(NULL)
        x1 <- x1[ix, , drop = FALSE]; n1 <- n1[ix, , drop = FALSE]
        x2 <- x2[ix, , drop = FALSE]; n2 <- n2[ix, , drop = FALSE]
        ac <- ac[ix]; ap <- ap[ix]; keep <- b$keep[ix]

        nreps1 <- ncol(x1); nreps2 <- ncol(x2)

        ## ---- DSS engine, matrices in, no BSseq ----
        if (!smoothing) {
            ep1 <- mean_ns(x1, n1); ep2 <- mean_ns(x2, n2)
        } else {
            ep1 <- mean_sm(x1, n1, ac, ap, smoothing.span)
            ep2 <- mean_sm(x2, n2, ac, ap, smoothing.span)
        }

        if (equal.disp || nreps1 == 1L || nreps2 == 1L) {
            phi1 <- phi2 <- est_disp(cbind(x1, x2), cbind(n1, n2), cbind(ep1, ep2), serial)
        } else {
            phi1 <- est_disp(x1, n1, ep1, serial)
            phi2 <- est_disp(x2, n2, ep2, serial)
        }

        wt1 <- 1 / (1 + (n1 - 1) * phi1); wt1 <- wt1 / mean(wt1)
        wt2 <- 1 / (1 + (n2 - 1) * phi2); wt2 <- wt2 / mean(wt2)
        x1w <- x1 * wt1; n1w <- n1 * wt1
        x2w <- x2 * wt2; n2w <- n2 * wt2

        if (!smoothing) {
            ep1 <- mean_ns(x1w, n1w); ep2 <- mean_ns(x2w, n2w)
            w   <- wald_ns(ep1[, 1], ep2[, 1], n1w, n2w, phi1, phi2)
        } else {
            ep1 <- mean_sm(x1w, n1w, ac, ap, smoothing.span)
            ep2 <- mean_sm(x2w, n2w, ac, ap, smoothing.span)
            w   <- wald_sm(ep1[, 1], ep2[, 1], n1w, n2w, phi1, phi2,
                           smoothing.span, ac, ap)
        }

        res <- data.frame(chr = ac, pos = ap, w, stringsAsFactors = FALSE)
        res <- res[keep & !is.na(res$stat), , drop = FALSE]      # drop halo + NA sites
        res$fdr <- NULL                                          # FDR is global; set later
        res
    }

    if (verbose)
        message("-Streaming ", length(blocks), " block(s) over ",
                length(chr_levels), " chromosome(s); on-disk (HDF5) methrix, n_cores = ",
                n_cores)

    parts <- if (n_cores > 1L)
        parallel::mclapply(blocks, worker, mc.cores = n_cores)
    else
        lapply(blocks, worker)

    err <- vapply(parts, function(x) inherits(x, "try-error"), logical(1L))
    if (any(err))
        stop("Block processing failed: ",
             attr(parts[[which(err)[1]]], "condition")$message)

    parts <- parts[!vapply(parts, is.null, logical(1L))]
    if (length(parts) == 0L)
        stop("No CpG sites passed the coverage filter in any block.")

    out <- do.call(rbind, parts)
    out <- out[order(match(out$chr, chr_levels), out$pos), , drop = FALSE]
    out$fdr <- p.adjust(out$pval, method = "fdr")                # single genome-wide BH
    rownames(out) <- NULL
    class(out) <- c("data.frame", "DMLtest")
    if (verbose) message("-Done. Tested ", nrow(out), " CpG sites.")
    out
}

#' Differential methylation test (DSS) on an in-memory or on-disk \code{methrix}
#'
#' @details
#' Runs the DSS single-CpG Wald test (\code{DSS::DMLtest}) without converting the
#' \code{\link{methrix}} object to a \code{bsseq} object, dispatching on storage:
#' \itemize{
#'   \item \strong{in-memory} objects go to \code{.DMLtest_inmemory()}, which
#'     performs a single genome-wide dispersion estimate and reproduces stock DSS.
#'   \item \strong{on-disk} (HDF5-backed) objects go to \code{.DMLtest_ondisk()},
#'     which streams one chromosome per block so the full matrix is never
#'     realised. The dispersion-shrinkage prior is then estimated per block; with
#'     whole-chromosome blocks this is a near-exact approximation of the
#'     genome-wide result. Set \code{block_size} to sub-split large chromosomes; a
#'     halo of \code{max(smoothing.span, 200)} bp is added and trimmed so
#'     smoothing and coverage filtering stay correct at block edges. FDR is
#'     computed once over the full retained set of sites.
#' }
#' The result is column- and class-identical to \code{DSS::DMLtest} in both
#' cases and can be passed unchanged to \code{DSS::callDML} / \code{DSS::callDMR}.
#'
#' @param m \code{\link{methrix}} object (in-memory or HDF5-backed).
#' @param group1,group2 Sample names or column indices defining the two groups.
#' @param equal.disp Assume equal dispersion in both groups. Default \code{FALSE}.
#' @param smoothing Smooth methylation levels before testing. Default \code{FALSE}.
#' @param smoothing.span Smoothing window in bp. Default 500.
#' @param n_cores Cores for the parallel step (blocks for on-disk, dispersion for
#'   in-memory). \code{NULL} lets the in-memory backend auto-detect. Default 1.
#' @param block_size (on-disk only) max CpGs realised per block; \code{NULL}
#'   (default) = one whole chromosome per block. Ignored for in-memory objects.
#' @param verbose Report progress. Default \code{TRUE}.
#'
#' @return A \code{data.frame} of class \code{c("data.frame", "DMLtest")} with
#'   columns \code{chr, pos, mu1, mu2, diff, diff.se, stat, phi1, phi2, pval,
#'   fdr}, sorted by chromosome then position.
#'
#' @examples
#' \dontrun{
#' data('methrix_data')
#' sn <- rownames(SummarizedExperiment::colData(methrix_data))
#' ## in-memory (dispatches to .DMLtest_inmemory)
#' r1 <- methrix_DMLtest(methrix_data, sn[1:3], sn[4:6], smoothing = TRUE)
#' ## on-disk (dispatches to .DMLtest_ondisk, streamed)
#' m_h5 <- convert_methrix(methrix_data)
#' r2 <- methrix_DMLtest(m_h5, sn[1:3], sn[4:6], smoothing = TRUE, n_cores = 2)
#' dmrs <- DSS::callDMR(r2, p.threshold = 0.05)
#' }
#' @seealso \code{\link{methrix2bsseq}}, \code{DSS::DMLtest}, \code{DSS::callDMR}
#' @export
methrix_DMLtest <- function(m, group1, group2,
                            equal.disp = FALSE, smoothing = FALSE,
                            smoothing.span = 500,
                            n_cores = 1L, block_size = NULL,
                            verbose = TRUE) {

    if (!is(m, "methrix"))
        stop("A valid methrix object needs to be supplied.")
    if (!requireNamespace("DSS", quietly = TRUE))
        stop("Package 'DSS' is required for methrix_DMLtest().")

    ## normalise groups to names so both backends receive the same form
    grp    <- .methrix_resolve_groups(m, group1, group2)

    if (!methrix:::is_h5(m)) {
        if (verbose) message("-in-memory methrix: using .DMLtest_inmemory()")
        return(.DMLtest_inmemory(m, group1 = grp$g1, group2 = grp$g2,
                                 equal.disp = equal.disp, smoothing = smoothing,
                                 smoothing.span = smoothing.span, ncores = n_cores))
    }

    .DMLtest_ondisk(m, group1 = grp$g1, group2 = grp$g2,
                    equal.disp = equal.disp, smoothing = smoothing,
                    smoothing.span = smoothing.span,
                    n_cores = n_cores, block_size = block_size, verbose = verbose)
}


################################################################################
##
##  methrix_DMLtest_inmemory.R
##
##  In-memory methrix → DSS internals (no BSseq).
##
##  For in-memory objects the assays are already plain matrices.
##  We extract them directly, compute M-counts, and call DSS's
##  non-exported functions which are fast (C code, mclapply).
##
##  DSS is a required dependency for this path.
##
################################################################################

#' @keywords internal
#' @importFrom GenomicRanges seqnames
.DMLtest_inmemory <- function(m, group1, group2, equal.disp, smoothing,
                              smoothing.span, ncores) {

  if (is.null(ncores))
    ncores <- if (.Platform$OS.type == "windows") 1L
  else max(parallel::detectCores() - 3L, 1L)

  all_s <- colnames(m)
  idx1 <- match(group1, all_s)
  idx2 <- match(group2, all_s)

  ## ---- Extract plain matrices (already in memory, instant) ----
  rr = as.data.frame(rowData(x = m))
  rr$seqnames <- as.character(rr$chr)
  rr$chr <- NULL

  if(is.null(rr$end)){
    rr$end <- rr$start + 1
  }

  rr <- GenomicRanges::makeGRangesFromDataFrame(rr, keep.extra.columns = FALSE)

  allchr <- as.character(GenomicRanges::seqnames(rr))
  allpos <- as.integer(GenomicRanges::start(rr))

  beta <- SummarizedExperiment::assay(m, "beta")
  cov  <- SummarizedExperiment::assay(m, "cov")

  cov[is.na(cov)]   <- 0L
  beta[is.na(beta)]  <- 0

  mcounts <- round(beta * cov)

  x1 <- mcounts[, idx1, drop = FALSE]
  n1 <- cov[, idx1, drop = FALSE]
  x2 <- mcounts[, idx2, drop = FALSE]
  n2 <- cov[, idx2, drop = FALSE]

  storage.mode(x1) <- "integer"; storage.mode(n1) <- "integer"
  storage.mode(x2) <- "integer"; storage.mode(n2) <- "integer"

  rm(beta, cov, mcounts)

  ## ---- Coverage filter via DSS internal ----
  DSS_ns <- asNamespace("DSS")
  hasCoverage_fn <- get("hasCoverage", envir = DSS_ns)

  keep <- hasCoverage_fn(n1, allpos) & hasCoverage_fn(n2, allpos)
  x1 <- x1[keep, , drop = FALSE]; n1 <- n1[keep, , drop = FALSE]
  x2 <- x2[keep, , drop = FALSE]; n2 <- n2[keep, , drop = FALSE]
  allchr <- allchr[keep]; allpos <- allpos[keep]

  nreps1 <- ncol(x1); nreps2 <- ncol(x2)

  if ((nreps1 == 1 || nreps2 == 1) && !equal.disp && !smoothing)
    stop("No replicates in at least one group. ",
         "Set smoothing=TRUE or equal.disp=TRUE.")

  ## ---- Call DSS internal functions ----
  compute.mean.noSmooth <- get("compute.mean.noSmooth", envir = DSS_ns)
  compute.mean.Smooth   <- get("compute.mean.Smooth",   envir = DSS_ns)
  est.dispersion.BSseq  <- get("est.dispersion.BSseq",  envir = DSS_ns)
  waldTest.DML          <- get("waldTest.DML",           envir = DSS_ns)

  if (!smoothing) {
    ep1 <- compute.mean.noSmooth(x1, n1)
    ep2 <- compute.mean.noSmooth(x2, n2)

    message("Estimating dispersion for each CpG site ...")
    if (equal.disp || nreps1 == 1 || nreps2 == 1) {
      phi1 <- phi2 <- est.dispersion.BSseq(
        cbind(x1, x2), cbind(n1, n2), cbind(ep1, ep2), ncores)
    } else {
      phi1 <- est.dispersion.BSseq(x1, n1, ep1, ncores)
      phi2 <- est.dispersion.BSseq(x2, n2, ep2, ncores)
    }

    wt1 <- 1 / (1 + (n1 - 1) * phi1); wt1 <- wt1 / mean(wt1)
    wt2 <- 1 / (1 + (n2 - 1) * phi2); wt2 <- wt2 / mean(wt2)
    x1w <- x1 * wt1; n1w <- n1 * wt1
    x2w <- x2 * wt2; n2w <- n2 * wt2

    ep1 <- compute.mean.noSmooth(x1w, n1w)
    ep2 <- compute.mean.noSmooth(x2w, n2w)

    result <- waldTest.DML(x1w, n1w, ep1, phi1,
                           x2w, n2w, ep2, phi2,
                           smoothing = FALSE,
                           allchr = allchr, allpos = allpos)
  } else {
    message("Smoothing ...")
    ep1 <- compute.mean.Smooth(x1, n1, allchr, allpos, smoothing.span)
    ep2 <- compute.mean.Smooth(x2, n2, allchr, allpos, smoothing.span)

    message("Estimating dispersion for each CpG site ...")
    if (equal.disp) {
      phi1 <- phi2 <- est.dispersion.BSseq(
        cbind(x1, x2), cbind(n1, n2), cbind(ep1, ep2), ncores)
    } else {
      phi1 <- est.dispersion.BSseq(x1, n1, ep1, ncores)
      phi2 <- est.dispersion.BSseq(x2, n2, ep2, ncores)
    }

    wt1 <- 1 / (1 + (n1 - 1) * phi1); wt1 <- wt1 / mean(wt1)
    wt2 <- 1 / (1 + (n2 - 1) * phi2); wt2 <- wt2 / mean(wt2)
    x1w <- x1 * wt1; n1w <- n1 * wt1
    x2w <- x2 * wt2; n2w <- n2 * wt2

    ep1 <- compute.mean.Smooth(x1w, n1w, allchr, allpos, smoothing.span)
    ep2 <- compute.mean.Smooth(x2w, n2w, allchr, allpos, smoothing.span)

    message("Computing test statistics ...")
    result <- waldTest.DML(x1w, n1w, ep1, phi1,
                           x2w, n2w, ep2, phi2,
                           smoothing = TRUE,
                           smoothing.span = smoothing.span,
                           allchr = allchr, allpos = allpos)
  }

  class(result) <- c("data.frame", "DMLtest")
  message("Done. Tested ", nrow(result), " CpG sites.")
  invisible(result)
}

################################################################################
## Multi-factor extension (outline)
##
## DMLfit.multiFactor / DMLtest.multiFactor follow the same block pattern: the
## fit (DSS:::DMLfit.oneCG) is per-CpG and mean smoothing is per-chromosome, so
## the identical chromosome-streaming driver applies. The accumulator differs:
## concatenate each block's fit matrix (fit$fit) and GRanges (fit$gr), keep a
## single shared design/formula/X, then run the coef/contrast test on the
## assembled object. Its p-value adjustment is also global -> apply once.
################################################################################


