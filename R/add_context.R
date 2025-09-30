#' Add Reference Base and Sequence Context to Methrix Object
#'
#' @description
#' Adds reference base and sequence context information to an existing methrix object.
#' Supports both FASTA files and BSgenome objects as reference sources.
#'
#' @param m A methrix object
#' @param ref_fasta Path to reference FASTA file (.fa or .fasta). Must have a .fai index
#'   (create with `samtools faidx`). Required if `bsgenome` is not provided.
#' @param bsgenome A BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38). Required if
#'   `ref_fasta` is not provided.
#' @param verbose Print progress messages (default: TRUE)
#'
#' @details
#' This function extracts the reference base and sequence context (CG, CHG, CHH)
#' for each CpG site in the methrix object and adds them to the rowData.
#'
#' When using a FASTA file, the function uses efficient C code with 50KB chunk caching
#' for fast extraction.
#'
#' When using a BSgenome object, the function extracts sequences using R/Bioconductor
#' functions.
#'
#' The context is determined from a 3-base window starting at each position, with
#' reverse complement handling for minus strand sites.
#'
#' @return A methrix object with `ref_base` and `context` columns added to rowData
#'
#' @examples
#' \dontrun{
#' # Using FASTA file
#' m <- read_modkit_v2(files = bedfiles, chrom_sizes = "hg38.chrom.sizes")
#' m <- add_context(m, ref_fasta = "hg38.fa")
#'
#' # Using BSgenome
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' m <- add_context(m, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
#' }
#'
#' @export
add_context <- function(m, ref_fasta = NULL, bsgenome = NULL, verbose = TRUE) {

  # Validate inputs
  if (!is(m, "methrix")) {
    stop("'m' must be a methrix object")
  }

  if (is.null(ref_fasta) && is.null(bsgenome)) {
    stop("Either 'ref_fasta' or 'bsgenome' must be provided")
  }

  if (!is.null(ref_fasta) && !is.null(bsgenome)) {
    stop("Provide only one of 'ref_fasta' or 'bsgenome', not both")
  }

  # Get coordinates from methrix object
  gr <- SummarizedExperiment::rowRanges(m)
  chr <- as.character(GenomicRanges::seqnames(gr))
  pos <- GenomicRanges::start(gr)  # Already 0-based in methrix
  strand <- as.character(GenomicRanges::strand(gr))

  if (verbose) {
    message(sprintf("Adding context for %d sites", length(chr)))
  }

  # Extract context using appropriate method
  if (!is.null(ref_fasta)) {
    # Use C implementation for FASTA files (fast)
    if (!file.exists(ref_fasta)) {
      stop("FASTA file not found: ", ref_fasta)
    }

    fai_file <- paste0(ref_fasta, ".fai")
    if (!file.exists(fai_file)) {
      stop("FASTA index not found: ", fai_file, "\n",
           "Create it with: samtools faidx ", ref_fasta)
    }

    result <- .Call("add_context_c",
                    chr,
                    as.integer(pos),
                    strand,
                    ref_fasta,
                    as.logical(verbose),
                    PACKAGE = "methrix")

    ref_base <- result$ref_base
    context <- result$context

  } else {
    # Use BSgenome (R implementation)
    if (verbose) {
      message("Extracting context from BSgenome...")
    }

    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("Package 'BSgenome' is required for using BSgenome objects")
    }

    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("Package 'Biostrings' is required for using BSgenome objects")
    }

    # Process by chromosome for efficiency
    ref_base <- character(length(chr))
    context <- character(length(chr))

    unique_chrs <- unique(chr)
    sites_processed <- 0

    for (cur_chr in unique_chrs) {
      chr_idx <- which(chr == cur_chr)

      if (verbose) {
        message(sprintf("Processing %s...", cur_chr))
      }

      # Extract reference base
      ref_base[chr_idx] <- as.character(BSgenome::getSeq(bsgenome,
                                                          cur_chr,
                                                          pos[chr_idx] + 1,  # BSgenome is 1-based
                                                          pos[chr_idx] + 1))

      # Extract 3-base context
      context_seq <- BSgenome::getSeq(bsgenome,
                                      cur_chr,
                                      pos[chr_idx] + 1,
                                      pos[chr_idx] + 3)

      # Handle strand
      for (i in seq_along(chr_idx)) {
        idx <- chr_idx[i]
        seq <- as.character(context_seq[i])

        # Reverse complement for minus strand
        if (strand[idx] == "-") {
          seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
        }

        # Determine context
        if (nchar(seq) >= 3) {
          if (substr(seq, 1, 2) == "CG") {
            context[idx] <- "CG"
          } else if (substr(seq, 1, 1) == "C" && substr(seq, 3, 3) == "G") {
            context[idx] <- "CHG"
          } else if (substr(seq, 1, 1) == "C") {
            context[idx] <- "CHH"
          } else {
            context[idx] <- "CNN"
          }
        } else {
          context[idx] <- "NNN"
        }
      }

      sites_processed <- sites_processed + length(chr_idx)

      # Periodic interrupt check
      if (sites_processed %% 10000 == 0) {
        if (verbose) {
          message(sprintf("  Processed %d sites", sites_processed))
        }
      }
    }

    if (verbose) {
      message("Context extraction complete")
    }
  }

  # Add to rowData
  SummarizedExperiment::rowData(m)$ref_base <- ref_base
  SummarizedExperiment::rowData(m)$context <- context

  if (verbose) {
    message("Added ref_base and context to rowData")
  }

  return(m)
}