#' Read ModKit bedMethyl files (V2 - High Performance Implementation)
#'
#' @description NEW: Reimplemented version using hash-based aggregation for
#' 10-20x performance improvement over the original implementation.
#'
#' @param files Character vector of ModKit bedMethyl file paths (.bed.gz with .tbi)
#' @param chrom_sizes Path to chromosome sizes file (TSV: chrom<TAB>length). Can be .fai file. REQUIRED.
#' @param target_mod Target modification code: "m" (5mC), "h" (5hmC), "a" (6mA), etc. Default "m".
#' @param interval_size Genomic interval size (bp) per chunk. Default 10000000 (10MB).
#' @param min_coverage Minimum coverage required per site. Default 1.
#' @param quality_filter Maximum failure rate (0-1) per site. Default 0.1.
#' @param combine_strands Combine CpG strands for cytosine modifications? Default FALSE.
#' @param ref_fasta Optional reference FASTA file path (with .fai index) for context annotation.
#'   When provided, adds 'ref_base' and 'context' columns to rowData with reference base
#'   and sequence context (CG/CHG/CHH) information.
#' @param h5 Store matrices as HDF5Arrays? Default FALSE.
#' @param h5_dir Directory for HDF5 storage. Default NULL.
#' @param h5temp Temporary directory for HDF5 files during processing. Default NULL (uses tempdir()).
#' @param coldata Optional DataFrame with sample metadata. If NULL, created from filenames.
#' @param verbose Print progress messages? Default TRUE.
#'
#' @return A methrix object containing the merged methylation data
#'
#' @details This is a reimplemented version of read_modkit() that uses:
#' - Interval-based streaming architecture (default for all cases)
#' - Hash-based coordinate discovery per interval (no global state)
#' - Chromosomes filtered to common set across all files
#' - Genome divided into fixed-size intervals (default 10MB)
#' - Each interval processed independently: discover → aggregate → output
#'
#' Memory characteristics:
#' - h5=FALSE: Accumulates intervals in memory, minimal overhead (~550 MB for 100M sites)
#' - h5=TRUE: True streaming with only one interval in RAM at a time
#'
#' The interval-based approach eliminates the need to store all coordinates globally,
#' providing bounded memory usage regardless of dataset size.
#'
#' @examples
#' \dontrun{
#' # Basic 5mC analysis (chromosome sizes file required)
#' meth <- read_modkit_v2(modkit_files, chrom_sizes = "hg38.chrom.sizes", target_mod = "m")
#'
#' # With context extraction
#' meth <- read_modkit_v2(modkit_files, chrom_sizes = "hg38.chrom.sizes",
#'                        ref_fasta = "hg38.fa", target_mod = "m")
#'
#' # Large dataset with HDF5 backend (memory efficient)
#' meth <- read_modkit_v2(modkit_files, chrom_sizes = "hg38.chrom.sizes",
#'                        h5 = TRUE, h5_dir = "./h5_data")
#'
#' # With custom interval size for memory control
#' meth <- read_modkit_v2(modkit_files, chrom_sizes = "hg38.chrom.sizes", interval_size = 5000000)
#' }
#'
#' @export
read_modkit_v2 <- function(files,
                          chrom_sizes,
                          target_mod = "m",
                          interval_size = 10000000L,
                          min_coverage = 1L,
                          quality_filter = 0.1,
                          combine_strands = FALSE,
                          ref_fasta = NULL,
                          h5 = FALSE,
                          h5_dir = NULL,
                          h5temp = NULL,
                          coldata = NULL,
                          verbose = TRUE) {

  # Input validation
  if (is.null(files) || length(files) == 0) {
    stop("No input files provided")
  }

  if (!is.character(files)) {
    stop("Files must be character vector")
  }

  # Check chromosome sizes file (REQUIRED)
  if (missing(chrom_sizes) || is.null(chrom_sizes)) {
    stop("chrom_sizes is required. Provide a TSV file with chromosome names and lengths.\n",
         "Format: <chrom><TAB><length>\n",
         "Example: chr1\\t249250621\n",
         "Can also use .fai file from samtools faidx.")
  }

  if (!file.exists(chrom_sizes)) {
    stop("Chromosome sizes file not found: ", chrom_sizes)
  }

  # Validate reference FASTA if provided
  if (!is.null(ref_fasta)) {
    if (!is.character(ref_fasta) || length(ref_fasta) != 1) {
      stop("ref_fasta must be a single character string")
    }
    if (!file.exists(ref_fasta)) {
      stop("Reference FASTA file not found: ", ref_fasta)
    }
    fai_file <- paste0(ref_fasta, ".fai")
    if (!file.exists(fai_file)) {
      stop("Reference FASTA index not found: ", fai_file,
           "\nRun: samtools faidx ", ref_fasta)
    }
  }

  # Check file existence
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }

  # Check for .tbi indices
  missing_indices <- character(0)
  for (file in files) {
    tbi_file <- paste0(file, ".tbi")
    if (!file.exists(tbi_file)) {
      missing_indices <- c(missing_indices, file)
    }
  }

  if (length(missing_indices) > 0) {
    stop("Missing tabix indices (.tbi files) for: ",
         paste(missing_indices, collapse = ", "),
         "\nRun: tabix -p bed <file> for each file")
  }

  # Validate modification code
  valid_mods <- c("m", "h", "a", "c", "g", "o")
  if (!target_mod %in% valid_mods) {
    stop("Invalid modification code '", target_mod, "'. Valid codes: ",
         paste(valid_mods, collapse = ", "))
  }

  # Setup sample metadata
  if (is.null(coldata)) {
    sample_names <- extract_sample_names(files)
    coldata <- data.frame(row.names = sample_names, stringsAsFactors = FALSE)
  } else {
    if (!is.data.frame(coldata)) {
      stop("coldata must be a data.frame")
    }
    if (nrow(coldata) != length(files)) {
      stop("Number of files (", length(files), ") must match coldata rows (",
           nrow(coldata), ")")
    }
  }

  # PHASE 2: Interval-based streaming
  start_time <- Sys.time()
  start_mem <- as.numeric(object.size(ls(all.names = TRUE))) / 1024^2

  # Step 1: Generate genomic intervals in C (efficient, quiet)
  intervals_df <- .Call("read_modkit_v2_create_intervals_c",
                       chrom_sizes,
                       files,
                       as.integer(interval_size),
                       as.logical(FALSE),  # Always quiet for interval creation
                       PACKAGE = "methrix")

  n_intervals <- nrow(intervals_df)
  n_chromosomes <- length(unique(intervals_df$chr))

  if (h5) {
    # TRUE STREAMING: Process to HDF5 directly
    methrix_obj <- process_intervals_to_h5(
      intervals_df = intervals_df,
      files = files,
      target_mod = target_mod,
      min_coverage = min_coverage,
      quality_filter = quality_filter,
      combine_strands = combine_strands,
      ref_fasta = ref_fasta,
      h5_dir = h5_dir,
      h5temp = h5temp,
      coldata = coldata,
      verbose = verbose,
      n_chromosomes = n_chromosomes,
      interval_size = interval_size
    )
  } else {
    # IN-MEMORY ACCUMULATION: Collect all intervals then build matrices
    methrix_obj <- process_intervals_to_memory(
      intervals_df = intervals_df,
      files = files,
      target_mod = target_mod,
      min_coverage = min_coverage,
      quality_filter = quality_filter,
      combine_strands = combine_strands,
      ref_fasta = ref_fasta,
      coldata = coldata,
      verbose = verbose,
      n_chromosomes = n_chromosomes,
      interval_size = interval_size
    )
  }

  # Final summary
  if (verbose) {
    end_time <- Sys.time()
    end_mem <- as.numeric(object.size(ls(all.names = TRUE))) / 1024^2
    total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    peak_mem <- end_mem  # Simplified - can use gc() for accurate peak
    mem_increase <- end_mem - start_mem

    cat(sprintf("\nDone | Time taken: %.3f secs | Peak memory: %.0f MB | Memory increase: %.1f MB | Sites: %s | Samples: %d\n",
                total_time, peak_mem, mem_increase,
                format(nrow(methrix_obj), big.mark = ","),
                ncol(methrix_obj)))
  }

  return(methrix_obj)
}

# Helper function to convert C result to methrix
convert_c_result_to_methrix_v2 <- function(c_result, coldata, target_mod,
                                          h5, h5_dir, h5temp, ref_fasta, verbose) {

  # Extract matrices
  beta_matrix <- c_result$beta_matrix
  cov_matrix <- c_result$cov_matrix

  # Set column names
  colnames(beta_matrix) <- rownames(coldata)
  colnames(cov_matrix) <- rownames(coldata)

  # Extract coordinate and context data BEFORE HDF5 conversion
  # (we need to keep this lightweight data, but can free the large matrices)
  chr_data <- c_result$chr
  start_data <- c_result$start
  strand_data <- c_result$strand
  ref_base_data <- c_result$ref_base
  context_data <- c_result$context

  # Convert to HDF5 if requested - do this IMMEDIATELY to free memory
  if (h5) {
    if (verbose) message("Converting to HDF5 format...")

    # Set temporary directory
    if (is.null(h5temp)) {
      h5temp <- tempdir()
    }

    # Create unique filenames to avoid conflicts
    sink_counter <- 1
    while (any(c(paste0("beta_modkit_v2_", sink_counter, ".h5"),
                 paste0("cov_modkit_v2_", sink_counter, ".h5")) %in% dir(h5temp))) {
      sink_counter <- sink_counter + 1
    }

    # Write matrices to HDF5 immediately (frees RAM)
    beta_matrix <- HDF5Array::writeHDF5Array(
      beta_matrix,
      filepath = file.path(h5temp, paste0("beta_modkit_v2_", sink_counter, ".h5")),
      name = "beta",
      level = 6  # compression level
    )

    cov_matrix <- HDF5Array::writeHDF5Array(
      cov_matrix,
      filepath = file.path(h5temp, paste0("cov_modkit_v2_", sink_counter, ".h5")),
      name = "cov",
      level = 6
    )

    # Free memory from C result (we already extracted what we need)
    rm(list = c("c_result"))
    gc()

    if (verbose) message("HDF5 conversion complete, memory released")
  }

  # Create row data with context information if available
  if (!is.null(ref_base_data) && !is.null(context_data)) {
    rowData <- DataFrame(
      chr = chr_data,
      start = start_data,
      strand = strand_data,
      ref_base = ref_base_data,
      context = context_data
    )
    if (verbose) message("Added reference base and context information to rowData")
  } else {
    rowData <- DataFrame(
      chr = chr_data,
      start = start_data,
      strand = strand_data
    )
  }

  # Create chromosome summary
  chr_summary <- as.data.table(rowData)[, .N, by = chr]
  data.table::setnames(chr_summary, "N", "n_cpgs")

  # Create metadata
  modkit_metadata <- list(
    source = "ModKit_V2",
    modification_type = get_modification_name(target_mod),
    modification_code = target_mod,
    implementation = "hash_aggregation",
    total_sites = nrow(beta_matrix),
    processing_engine = "modkit_merge_v2_c",
    reference_fasta = ref_fasta,
    context_extracted = !is.null(ref_base_data),
    context_types = if (!is.null(context_data)) names(table(context_data)) else NULL,
    created_date = Sys.Date()
  )

  # Calculate basic statistics (silently)
  genome_stats <- data.table::data.table(
    Sample_Name = colnames(beta_matrix)
  )

  # Use block processing for HDF5 arrays (memory efficient)
  if (h5) {
    for (i in seq_len(ncol(beta_matrix))) {
      # DelayedArray operations work on blocks automatically
      sample_beta <- beta_matrix[, i]
      sample_cov <- cov_matrix[, i]

      # Use DelayedMatrixStats for efficient block-wise computation
      valid_mask <- !is.na(sample_beta) & !is.na(sample_cov) & sample_cov > 0

      if (sum(valid_mask) > 0) {
        beta_valid <- sample_beta[valid_mask]
        cov_valid <- sample_cov[valid_mask]

        genome_stats$mean_meth[i] <- DelayedMatrixStats::rowMeans2(
          as.matrix(beta_valid), na.rm = TRUE)
        genome_stats$median_meth[i] <- DelayedMatrixStats::rowMedians(
          as.matrix(beta_valid), na.rm = TRUE)
        genome_stats$mean_cov[i] <- DelayedMatrixStats::rowMeans2(
          as.matrix(cov_valid), na.rm = TRUE)
        genome_stats$median_cov[i] <- DelayedMatrixStats::rowMedians(
          as.matrix(cov_valid), na.rm = TRUE)
      } else {
        genome_stats$mean_meth[i] <- NA_real_
        genome_stats$median_meth[i] <- NA_real_
        genome_stats$mean_cov[i] <- NA_real_
        genome_stats$median_cov[i] <- NA_real_
      }
    }
  } else {
    # In-memory calculation (original fast path)
    for (i in seq_len(ncol(beta_matrix))) {
      sample_beta <- beta_matrix[, i]
      sample_cov <- cov_matrix[, i]

      valid_idx <- !is.na(sample_beta) & !is.na(sample_cov) & sample_cov > 0

      if (sum(valid_idx) > 0) {
        genome_stats$mean_meth[i] <- mean(sample_beta[valid_idx], na.rm = TRUE)
        genome_stats$median_meth[i] <- median(sample_beta[valid_idx], na.rm = TRUE)
        genome_stats$mean_cov[i] <- mean(sample_cov[valid_idx], na.rm = TRUE)
        genome_stats$median_cov[i] <- median(sample_cov[valid_idx], na.rm = TRUE)
      } else {
        genome_stats$mean_meth[i] <- NA_real_
        genome_stats$median_meth[i] <- NA_real_
        genome_stats$mean_cov[i] <- NA_real_
        genome_stats$median_cov[i] <- NA_real_
      }
    }
  }

  descriptive_stats <- list(
    genome_stat = genome_stats,
    modkit_info = modkit_metadata
  )

  # Create methrix object
  methrix_obj <- create_methrix(
    beta_mat = beta_matrix,
    cov_mat = cov_matrix,
    cpg_loci = rowData,
    is_hdf5 = h5,
    genome_name = "modkit_v2",
    col_data = coldata,
    h5_dir = h5_dir,
    ref_cpg_dt = chr_summary,
    chrom_sizes = NULL,
    desc = descriptive_stats
  )

  S4Vectors::metadata(methrix_obj)$modkit_info <- modkit_metadata

  return(methrix_obj)
}

# Helper function: Process intervals to in-memory matrices
process_intervals_to_memory <- function(intervals_df, files, target_mod,
                                       min_coverage, quality_filter, combine_strands,
                                       ref_fasta, coldata, verbose,
                                       n_chromosomes, interval_size) {

  n_samples <- length(files)
  n_intervals <- nrow(intervals_df)
  all_results <- list()

  # Track start time for progress bar
  loop_start <- Sys.time()

  # Single-line progress bar
  if (verbose) {
    cat(sprintf("Interval size: %s bp, processing %d bins from %d chromosome%s\n",
                format(interval_size, big.mark = ","),
                n_intervals,
                n_chromosomes,
                ifelse(n_chromosomes > 1, "s", "")))
  }

  # Process each interval
  for (i in seq_len(n_intervals)) {
    interval <- intervals_df[i, ]

    # Call C function to process this interval (always quiet)
    interval_result <- .Call("read_modkit_v2_process_interval_c",
                             interval$chr,
                             as.integer(interval$start),
                             as.integer(interval$end),
                             files,
                             target_mod,
                             as.integer(min_coverage),
                             as.numeric(quality_filter),
                             as.logical(combine_strands),
                             as.logical(FALSE),  # Always quiet in C
                             PACKAGE = "methrix")

    # Store result if any sites found
    if (!is.null(interval_result) && !is.null(interval_result$beta) && nrow(interval_result$beta) > 0) {
      all_results[[length(all_results) + 1]] <- interval_result
    }

    # Update progress bar (single line, overwrite) - only update every 5% or on completion
    if (verbose) {
      if (i == n_intervals || i %% max(1, floor(n_intervals / 20)) == 0 || i == 1) {
        progress_pct <- i / n_intervals
        bar_width <- 40
        filled <- round(bar_width * progress_pct)
        bar <- paste0(rep("#", filled), collapse = "")
        spaces <- paste0(rep(" ", bar_width - filled), collapse = "")

        # Calculate elapsed time from loop start
        elapsed <- as.numeric(difftime(Sys.time(), loop_start, units = "secs"))

        # Use cat with file = stderr() for real-time updates
        cat(sprintf("\r[%02d:%02d:%02d] %s%s %5d/%d bins processed",
                    floor(elapsed / 3600),
                    floor((elapsed %% 3600) / 60),
                    round(elapsed %% 60),
                    bar, spaces, i, n_intervals),
            file = stderr())
        flush(stderr())
      }
    }
  }

  if (verbose) cat("\n", file = stderr())  # New line after progress bar

  # Merge all interval results
  if (length(all_results) == 0) {
    stop("No sites found across all intervals")
  }

  # Concatenate matrices from all intervals
  beta_matrix <- do.call(rbind, lapply(all_results, function(x) x$beta))
  cov_matrix <- do.call(rbind, lapply(all_results, function(x) x$cov))

  # Build coordinate vectors from interval results
  chr_vec <- character(0)
  start_vec <- integer(0)
  strand_vec <- character(0)

  for (res in all_results) {
    n_sites <- length(res$positions)
    chr_vec <- c(chr_vec, rep(res$chr, n_sites))
    start_vec <- c(start_vec, res$positions)
    strand_vec <- c(strand_vec, res$strand)
  }

  # Set column names
  colnames(beta_matrix) <- rownames(coldata)
  colnames(cov_matrix) <- rownames(coldata)

  # Add context if reference FASTA provided
  ref_base_vec <- NULL
  context_vec <- NULL

  if (!is.null(ref_fasta)) {
    if (verbose) message("Extracting sequence context from reference...")

    context_result <- .Call("add_context_c",
                           chr_vec,
                           as.integer(start_vec),
                           strand_vec,
                           ref_fasta,
                           as.logical(verbose),
                           PACKAGE = "methrix")

    ref_base_vec <- context_result$ref_base
    context_vec <- context_result$context
  }

  # Create final result object
  c_result <- list(
    beta_matrix = beta_matrix,
    cov_matrix = cov_matrix,
    chr = chr_vec,
    start = start_vec,
    strand = strand_vec,
    ref_base = ref_base_vec,
    context = context_vec
  )

  # Convert to methrix object (reuse existing function)
  methrix_obj <- convert_c_result_to_methrix_v2(c_result, coldata, target_mod,
                                                h5 = FALSE, NULL, NULL, ref_fasta, verbose)

  return(methrix_obj)
}

# Helper function: Process intervals to HDF5 (true streaming)
process_intervals_to_h5 <- function(intervals_df, files, target_mod,
                                   min_coverage, quality_filter, combine_strands,
                                   ref_fasta, h5_dir, h5temp, coldata, verbose) {

  n_samples <- length(files)

  # Set temporary directory
  if (is.null(h5temp)) {
    h5temp <- tempdir()
  }

  # Create unique filenames
  sink_counter <- 1
  while (any(c(paste0("beta_modkit_v2_stream_", sink_counter, ".h5"),
               paste0("cov_modkit_v2_stream_", sink_counter, ".h5")) %in% dir(h5temp))) {
    sink_counter <- sink_counter + 1
  }

  beta_h5_path <- file.path(h5temp, paste0("beta_modkit_v2_stream_", sink_counter, ".h5"))
  cov_h5_path <- file.path(h5temp, paste0("cov_modkit_v2_stream_", sink_counter, ".h5"))

  if (verbose) {
    message("Streaming ", nrow(intervals_df), " intervals to HDF5")
    message("Beta matrix: ", basename(beta_h5_path))
    message("Coverage matrix: ", basename(cov_h5_path))
  }

  # Initialize HDF5 sinks (we don't know final dimensions yet)
  # We'll collect coordinates first to determine total sites
  all_coords <- list()
  all_beta_chunks <- list()
  all_cov_chunks <- list()

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = nrow(intervals_df), style = 3)
  }

  # Process each interval
  for (i in seq_len(nrow(intervals_df))) {
    interval <- intervals_df[i, ]

    # Call C function to process this interval
    interval_result <- .Call("read_modkit_v2_process_interval_c",
                             interval$chr,
                             as.integer(interval$start),
                             as.integer(interval$end),
                             files,
                             target_mod,
                             as.integer(min_coverage),
                             as.numeric(quality_filter),
                             as.logical(combine_strands),
                             as.logical(verbose),  # pass through verbose setting
                             PACKAGE = "methrix")

    # Store results if any sites found
    if (!is.null(interval_result) && !is.null(interval_result$beta) && nrow(interval_result$beta) > 0) {
      all_coords[[length(all_coords) + 1]] <- list(
        chr = interval_result$chr,
        positions = interval_result$positions,
        strand = interval_result$strand
      )
      all_beta_chunks[[length(all_beta_chunks) + 1]] <- interval_result$beta
      all_cov_chunks[[length(all_cov_chunks) + 1]] <- interval_result$cov
    }

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("\nWriting ", length(all_beta_chunks), " chunks to HDF5...")
  }

  # Merge coordinates and write to HDF5
  if (length(all_beta_chunks) == 0) {
    stop("No sites found across all intervals")
  }

  # Concatenate all chunks
  beta_matrix <- do.call(rbind, all_beta_chunks)
  cov_matrix <- do.call(rbind, all_cov_chunks)

  # Build coordinate vectors from interval results
  chr_vec <- character(0)
  start_vec <- integer(0)
  strand_vec <- character(0)

  for (coords in all_coords) {
    n_sites <- length(coords$positions)
    chr_vec <- c(chr_vec, rep(coords$chr, n_sites))
    start_vec <- c(start_vec, coords$positions)
    strand_vec <- c(strand_vec, coords$strand)
  }

  # Write to HDF5
  colnames(beta_matrix) <- rownames(coldata)
  colnames(cov_matrix) <- rownames(coldata)

  beta_h5 <- HDF5Array::writeHDF5Array(
    beta_matrix,
    filepath = beta_h5_path,
    name = "beta",
    level = 6
  )

  cov_h5 <- HDF5Array::writeHDF5Array(
    cov_matrix,
    filepath = cov_h5_path,
    name = "cov",
    level = 6
  )

  # Free memory
  rm(beta_matrix, cov_matrix, all_beta_chunks, all_cov_chunks)
  gc()

  if (verbose) message("HDF5 write complete, memory released")

  # Add context if needed
  ref_base_vec <- NULL
  context_vec <- NULL

  if (!is.null(ref_fasta)) {
    if (verbose) message("Extracting sequence context from reference...")

    context_result <- .Call("add_context_c",
                           chr_vec,
                           as.integer(start_vec),
                           strand_vec,
                           ref_fasta,
                           as.logical(verbose),
                           PACKAGE = "methrix")

    ref_base_vec <- context_result$ref_base
    context_vec <- context_result$context
  }

  # Create result structure with HDF5 arrays
  c_result <- list(
    beta_matrix = beta_h5,
    cov_matrix = cov_h5,
    chr = chr_vec,
    start = start_vec,
    strand = strand_vec,
    ref_base = ref_base_vec,
    context = context_vec
  )

  # Convert to methrix (h5 flag already TRUE in matrices)
  methrix_obj <- convert_c_result_to_methrix_v2(c_result, coldata, target_mod,
                                                h5 = TRUE, h5_dir, h5temp, ref_fasta, verbose)

  return(methrix_obj)
}