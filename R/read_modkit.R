#' Read ModKit bedMethyl files and create methrix object
#' 
#' @description Reads ModKit bedMethyl files using tabix streaming and creates a methrix object. 
#' 
#' @param files Character vector of ModKit bedMethyl file paths (.bed.gz with .tbi)
#' @param target_mod Target modification code: "m" (5mC), "h" (5hmC), "a" (6mA), "c" (5fC), "g" (5caC). Default "m".
#' @param bin_size Bin size for fixed genomic binning in bases. Default 50000000 (50MB).
#' @param quality_filter Maximum failure rate (0-1) per site. Default 0.1.
#' @param min_coverage Minimum coverage required per site. Default 1.
#' @param combine_strands Combine CpG strands for cytosine modifications? Default FALSE.
#' @param ref_fasta Optional reference FASTA file path (with .fai index) for context annotation.
#'   When provided, adds 'ref_base' and 'context' columns to rowData with reference base
#'   and sequence context (CG/CHG/CHH) information. Required for use_fixed_bins=TRUE.
#' @param use_fixed_bins Use fixed genomic bins instead of data-driven bins? Default TRUE.
#'   When TRUE, requires ref_fasta to determine chromosome lengths.
#' @param n_cores Number of cores for parallel processing. Default 1.
#' @param h5 Store matrices as HDF5Arrays? Default FALSE.
#' @param h5_dir Directory for HDF5 storage. Default NULL.
#' @param coldata Optional DataFrame with sample metadata. If NULL, created from filenames.
#' @param auto_index Automatically create tabix indices for files missing .tbi files? Default TRUE.
#' @param verbose Print progress messages? Default TRUE.
#' 
#' @return A methrix object containing the merged methylation data
#' 
#' @details This function uses an optimized C implementation for merging ModKit
#' bedMethyl files. It automatically discovers chromosomes with data and uses
#' a streaming merge approach for memory efficiency. Requires BGziped bedmethyl files.
#' 
#' @examples
#' \dontrun{
#' # Basic 5mC analysis
#' meth <- read_modkit(modkit_files, target_mod = "m")
#' 
#' # 6mA analysis with parallel processing
#' meth_6mA <- read_modkit(modkit_files, target_mod = "a", n_cores = 4)
#' 
#' # With reference genome context
#' meth_ctx <- read_modkit(modkit_files, target_mod = "m", 
#'                        ref_fasta = "genome.fa")
#' 
#' # High-performance with HDF5 storage
#' meth_h5 <- read_modkit(modkit_files, target_mod = "m", 
#'                       h5 = TRUE, n_cores = 8)
#'                       
#' # Disable automatic indexing if needed
#' meth <- read_modkit(modkit_files, auto_index = FALSE)
#' }
#' 
#' @export
read_modkit <- function(files,
                        target_mod = "m",
                        bin_size = 50000000L,
                        quality_filter = 0.1,
                        min_coverage = 1L,
                        combine_strands = FALSE,
                        ref_fasta = NULL,
                        use_fixed_bins = TRUE,
                        n_cores = 1L,
                        h5 = FALSE,
                        h5_dir = NULL,
                        coldata = NULL,
                        auto_index = TRUE,
                        verbose = TRUE) {
  
  # Input validation
  if (is.null(files) || length(files) == 0) {
    stop("No input files provided")
  }
  
  if (!is.character(files)) {
    stop("Files must be character vector")
  }
  
  # Check file existence
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }
  
  # Check for .tbi indices and create if missing
  missing_indices <- character(0)
  for (file in files) {
    tbi_file <- paste0(file, ".tbi")
    if (!file.exists(tbi_file)) {
      missing_indices <- c(missing_indices, file)
    }
  }
  
  if (length(missing_indices) > 0) {
    if (auto_index) {
      if (verbose) {
        message("Creating tabix indices for ", length(missing_indices), " files...")
      }
      
      for (file in missing_indices) {
        if (verbose) {
          message("  Indexing ", basename(file), "...")
        }
        
        tryCatch({
          # Use Rsamtools to create tabix index
          Rsamtools::indexTabix(file, format = "bed", zeroBased = TRUE)
        }, error = function(e) {
          stop("Failed to create tabix index for ", file, ": ", e$message, 
               "\nEnsure file is bgzip compressed. Run: bgzip ", file)
        })
      }
      
      if (verbose) {
        message("Tabix indexing completed")
      }
      
    } else {
      stop("Missing tabix indices (.tbi files): ", paste(paste0(missing_indices, ".tbi"), collapse = ", "),
           "\nRun: tabix -p bed <file> for each file, or set auto_index = TRUE")
    }
  }
  
  # Validate modification code
  valid_mods <- c("m", "h", "a", "c", "g", "o")
  if (!target_mod %in% valid_mods) {
    stop("Invalid modification code '", target_mod, "'. Valid codes: ", 
         paste(valid_mods, collapse = ", "))
  }
  
  # Validate numeric parameters
  if (!is.numeric(quality_filter) || quality_filter < 0 || quality_filter > 1) {
    stop("quality_filter must be between 0 and 1")
  }
  
  if (!is.integer(min_coverage) || min_coverage < 1) {
    min_coverage <- as.integer(min_coverage)
    if (min_coverage < 1) stop("min_coverage must be >= 1")
  }
  
  if (!is.integer(bin_size) || bin_size < 1000) {
    bin_size <- as.integer(bin_size)
    if (bin_size < 1000) stop("bin_size must be >= 1000")
  }
  
  if (!is.integer(n_cores) || n_cores < 1) {
    n_cores <- as.integer(n_cores)
    if (n_cores < 1) stop("n_cores must be >= 1")
  }
  
  # Validate fixed binning requirements
  if (use_fixed_bins && is.null(ref_fasta)) {
    stop("use_fixed_bins=TRUE requires ref_fasta to determine chromosome lengths")
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
           ". Run: samtools faidx ", ref_fasta)
    }
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
    # Ensure valid R names
    rownames(coldata) <- make.names(rownames(coldata), unique = TRUE)
  }
  
  if (verbose) {
    message("Reading ModKit bedMethyl files...")
    message("Files: ", length(files))
    message("Target modification: ", get_modification_name(target_mod))
    message("Samples: ", nrow(coldata))
    if (!is.null(ref_fasta)) {
      message("Reference FASTA: ", basename(ref_fasta))
      message("Will extract reference base and sequence context")
    }
    if (n_cores > 1) {
      message("Parallel cores: ", n_cores)
    }
  }
  
  # Call C function
  start_time <- Sys.time()
  
  c_result <- .Call("read_modkit_c",
                    files,
                    target_mod,
                    bin_size,
                    n_cores,
                    quality_filter,
                    min_coverage,
                    combine_strands,
                    ref_fasta,
                    use_fixed_bins,
                    verbose, PACKAGE = "methrix")
  
  processing_time <- Sys.time() - start_time
  
  if (verbose) {
    message("C processing completed in: ", format(processing_time))
    message("Sites discovered: ", formatC(nrow(c_result$beta_matrix), format = "d", big.mark = ","))
    if (!is.null(c_result$ref_base)) {
      message("Reference context extracted for all sites")
    }
  }
  
  # Convert C results to methrix object
  methrix_obj <- convert_c_result_to_methrix(c_result, coldata, target_mod, 
                                             h5, h5_dir, ref_fasta, verbose)
  
  if (verbose) {
    total_time <- Sys.time() - start_time
    message("Total processing time: ", format(total_time))
    message("Created methrix object: ", nrow(methrix_obj), " sites × ", 
            ncol(methrix_obj), " samples")
    
    # Report context summary if available
    if (!is.null(c_result$ref_base)) {
      context_summary <- table(c_result$context)
      message("Context distribution:")
      for (ctx in names(context_summary)) {
        message("  ", ctx, ": ", formatC(context_summary[ctx], format = "d", big.mark = ","))
      }
    }
  }
  
  return(methrix_obj)
}


extract_sample_names <- function(files) {
  sample_names <- basename(files)
  # Remove common ModKit suffixes
  sample_names <- gsub("\\.wf_mods.*$", "", sample_names)
  sample_names <- gsub("\\.bed.*$", "", sample_names)
  sample_names <- gsub("\\.bedmethyl.*$", "", sample_names)
  sample_names <- make.names(sample_names, unique = TRUE)
  return(sample_names)
}

get_modification_name <- function(mod_code) {
  mod_names <- c(
    "m" = "5-methylcytosine (5mC)",
    "h" = "5-hydroxymethylcytosine (5hmC)", 
    "a" = "6-methyladenine (6mA)",
    "c" = "5-formylcytosine (5fC)",
    "g" = "5-carboxylcytosine (5caC)",
    "o" = "8-oxoguanine"
  )
  return(mod_names[mod_code])
}

convert_c_result_to_methrix <- function(c_result, coldata, target_mod, h5, h5_dir, ref_fasta, verbose) {
  
  # Extract matrices
  beta_matrix <- c_result$beta_matrix
  cov_matrix <- c_result$cov_matrix
  
  # Set column names
  colnames(beta_matrix) <- rownames(coldata)
  colnames(cov_matrix) <- rownames(coldata)
  
  # Convert to HDF5 if requested
  if (h5) {
    if (verbose) message("Converting to HDF5 format...")
    beta_matrix <- as(beta_matrix, "HDF5Array")
    cov_matrix <- as(cov_matrix, "HDF5Array")
  }
  
  # UPDATED: Create row data with context information if available
  if (!is.null(c_result$ref_base) && !is.null(c_result$context)) {
    rowData <- DataFrame(
      chr = c_result$chr,
      start = c_result$start,
      strand = c_result$strand,
      ref_base = c_result$ref_base,      # NEW
      context = c_result$context         # NEW
    )
    if (verbose) message("Added reference base and context information to rowData")
    
    # Validate context extraction for cytosine modifications
    if (target_mod %in% c("m", "h", "c", "g")) {
      non_c_sites <- sum(c_result$ref_base != "C" & c_result$strand == "+", na.rm = TRUE)
      if (non_c_sites > 0) {
        warning("Found ", non_c_sites, " sites where reference base is not C on + strand. ",
                "This may indicate coordinate system issues or non-cytosine modifications.")
      }
    }
  } else {
    rowData <- DataFrame(
      chr = c_result$chr,
      start = c_result$start,
      strand = c_result$strand
    )
  }
  
  # Create chromosome summary for metadata
  chr_summary <- as.data.table(rowData)[, .N, by = chr]
  data.table::setnames(chr_summary, "N", "n_cpgs")
  
  # Create chromosome sizes (approximate from data range)
  chrom_sizes <- as.data.table(rowData)[, .(
    contig = chr,
    length = max(start) - min(start) + 1L
  ), by = chr][, .(contig, length)]
  
  # Create ModKit-specific metadata
  modkit_metadata <- list(
    source = "ModKit",
    modification_type = get_modification_name(target_mod),
    modification_code = target_mod,
    coordinate_system = "union_observed",
    total_sites = nrow(beta_matrix),
    reference_fasta = ref_fasta,
    processing_engine = "modkit_merge_improved_c",
    created_date = Sys.Date(),
    # NEW: Context extraction info
    context_extracted = !is.null(c_result$ref_base),
    context_types = if (!is.null(c_result$context)) names(table(c_result$context)) else NULL
  )
  
  # Create basic descriptive statistics
  if (verbose) message("Calculating basic statistics...")
  
  # Calculate genome-wide stats per sample
  genome_stats <- data.table::data.table(
    Sample_Name = colnames(beta_matrix)
  )
  
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
  
  descriptive_stats <- list(
    genome_stat = genome_stats,
    modkit_info = modkit_metadata
  )
  
  # Create methrix object using existing infrastructure
  methrix_obj <- create_methrix(
    beta_mat = beta_matrix,
    cov_mat = cov_matrix,
    cpg_loci = rowData,
    is_hdf5 = h5,
    genome_name = ifelse(is.null(ref_fasta), yes = "modkit", no = basename(ref_fasta)),
    col_data = coldata,
    h5_dir = h5_dir,
    ref_cpg_dt = chr_summary,
    chrom_sizes = chrom_sizes,
    desc = descriptive_stats
  )
  
  # Add ModKit-specific metadata to the methrix object
  S4Vectors::metadata(methrix_obj)$modkit_info <- modkit_metadata
  
  return(methrix_obj)
}