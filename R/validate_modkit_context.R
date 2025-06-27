validate_context_extraction <- function(methrix_obj, target_mod, verbose = TRUE) {
  
  if (!"ref_base" %in% colnames(rowData(methrix_obj))) {
    return(list(
      has_context = FALSE,
      validation_passed = TRUE,
      message = "No context information available"
    ))
  }
  
  rd <- rowData(methrix_obj)
  
  # Check for cytosine modifications
  if (target_mod %in% c("m", "h", "c", "g")) {
    
    # Count non-cytosine sites on positive strand
    plus_strand <- rd$strand == "+"
    non_c_plus <- sum(rd$ref_base[plus_strand] != "C", na.rm = TRUE)
    total_plus <- sum(plus_strand, na.rm = TRUE)
    
    # Count non-guanine sites on negative strand (should be G for C on opposite strand)
    minus_strand <- rd$strand == "-"
    non_g_minus <- sum(rd$ref_base[minus_strand] != "G", na.rm = TRUE)
    total_minus <- sum(minus_strand, na.rm = TRUE)
    
    # Context distribution
    context_dist <- table(rd$context)
    
    validation_results <- list(
      has_context = TRUE,
      target_mod = target_mod,
      total_sites = nrow(rd),
      plus_strand_sites = total_plus,
      minus_strand_sites = total_minus,
      non_c_plus_strand = non_c_plus,
      non_g_minus_strand = non_g_minus,
      context_distribution = context_dist,
      validation_passed = (non_c_plus / total_plus < 0.01) && (non_g_minus / total_minus < 0.01)
    )
    
    if (verbose) {
      message("Context validation for ", get_modification_name(target_mod), ":")
      message("  Total sites: ", formatC(validation_results$total_sites, format = "d", big.mark = ","))
      message("  + strand: ", formatC(total_plus, format = "d", big.mark = ","), 
              " (", round(100 * total_plus / validation_results$total_sites, 1), "%)")
      message("  - strand: ", formatC(total_minus, format = "d", big.mark = ","), 
              " (", round(100 * total_minus / validation_results$total_sites, 1), "%)")
      
      if (non_c_plus > 0) {
        message("  WARNING: ", non_c_plus, " non-C sites on + strand (", 
                round(100 * non_c_plus / total_plus, 2), "%)")
      }
      if (non_g_minus > 0) {
        message("  WARNING: ", non_g_minus, " non-G sites on - strand (", 
                round(100 * non_g_minus / total_minus, 2), "%)")
      }
      
      message("  Context distribution:")
      for (ctx in names(context_dist)) {
        message("    ", ctx, ": ", formatC(context_dist[ctx], format = "d", big.mark = ","),
                " (", round(100 * context_dist[ctx] / validation_results$total_sites, 1), "%)")
      }
    }
    
  } else {
    # For non-cytosine modifications (like 6mA)
    validation_results <- list(
      has_context = TRUE,
      target_mod = target_mod,
      total_sites = nrow(rd),
      context_distribution = table(rd$context),
      validation_passed = TRUE,
      message = paste("Context extraction for", get_modification_name(target_mod))
    )
    
    if (verbose) {
      message("Context extraction for ", get_modification_name(target_mod), ":")
      message("  Total sites: ", formatC(validation_results$total_sites, format = "d", big.mark = ","))
      message("  Reference bases: ", paste(unique(rd$ref_base), collapse = ", "))
      message("  Context types: ", paste(names(validation_results$context_distribution), collapse = ", "))
    }
  }
  
  return(validation_results)
}

#' Get context summary for methrix object
#' 
#' @description Summarize sequence context information in a methrix object
#' 
#' @param m A methrix object with context information
#' @param by_chromosome Include per-chromosome breakdown? Default FALSE.
#' 
#' @return data.frame with context summary
#' @export
#' 
#' @examples
#' \dontrun{
#' # Read ModKit data with context
#' meth <- read_modkit(files, ref_fasta = "genome.fa")
#' 
#' # Get context summary
#' context_summary(meth)
#' 
#' # Get per-chromosome breakdown
#' context_summary(meth, by_chromosome = TRUE)
#' }
context_summary <- function(m, by_chromosome = FALSE) {
  
  if (!is(m, "methrix")) {
    stop("Input must be a methrix object")
  }
  
  rd <- rowData(m)
  
  if (!"context" %in% colnames(rd)) {
    stop("No context information found in methrix object. Use ref_fasta parameter in read_modkit()")
  }
  
  if (by_chromosome) {
    # Per-chromosome summary
    summary_dt <- as.data.table(rd)[, .N, by = .(chr, context)]
    summary_wide <- data.table::dcast(summary_dt, chr ~ context, value.var = "N", fill = 0)
    
    # Add totals
    context_cols <- setdiff(colnames(summary_wide), "chr")
    summary_wide$total <- rowSums(summary_wide[, ..context_cols], na.rm = TRUE)
    
    # Calculate percentages
    for (col in context_cols) {
      pct_col <- paste0(col, "_pct")
      summary_wide[[pct_col]] <- round(100 * summary_wide[[col]] / summary_wide$total, 1)
    }
    
    return(as.data.frame(summary_wide))
    
  } else {
    # Genome-wide summary
    context_counts <- table(rd$context)
    context_df <- data.frame(
      context = names(context_counts),
      count = as.numeric(context_counts),
      percentage = round(100 * as.numeric(context_counts) / sum(context_counts), 1),
      stringsAsFactors = FALSE
    )
    
    return(context_df)
  }
}

check_reference_compatibility <- function(ref_fasta, modkit_files, n_sites = 1000) {
  
  # Check file existence
  if (!file.exists(ref_fasta)) {
    stop("Reference FASTA file not found: ", ref_fasta)
  }
  
  fai_file <- paste0(ref_fasta, ".fai")
  if (!file.exists(fai_file)) {
    stop("Reference FASTA index not found: ", fai_file)
  }
  
  # Read a small sample without reference first
  if (length(modkit_files) > 0 && file.exists(modkit_files[1])) {
    
    message("Sampling ", n_sites, " sites from ", basename(modkit_files[1]), " for compatibility check...")
    
    # Read small sample without reference
    sample_data <- read_modkit(modkit_files[1], 
                               ref_fasta = NULL,
                               bin_size = 10000000L,  # Large bins for sampling
                               verbose = FALSE)
    
    # Sample random sites
    if (nrow(sample_data) > n_sites) {
      sample_idx <- sample(nrow(sample_data), n_sites)
      sample_data <- sample_data[sample_idx, ]
    }
    
    # Extract coordinates for testing
    coords <- rowData(sample_data)
    
    # Try to read same data with reference
    tryCatch({
      sample_with_ref <- read_modkit(modkit_files[1],
                                     ref_fasta = ref_fasta,
                                     bin_size = 10000000L,
                                     verbose = FALSE)
      
      message("Reference compatibility check: PASSED")
      return(list(
        compatible = TRUE,
        message = "Reference FASTA is compatible with ModKit data",
        sample_chromosomes = unique(coords$chr),
        total_sampled_sites = nrow(coords)
      ))
      
    }, error = function(e) {
      warning("Reference compatibility check failed: ", e$message)
      return(list(
        compatible = FALSE,
        message = paste("Reference compatibility check failed:", e$message),
        sample_chromosomes = unique(coords$chr),
        total_sampled_sites = nrow(coords)
      ))
    })
    
  } else {
    stop("No valid ModKit files provided for compatibility check")
  }
}

preview_modkit_file <- function(file, n_lines = 1000) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  # Check if tabix indexed
  tbi_file <- paste0(file, ".tbi")
  has_index <- file.exists(tbi_file)
  
  # Read first n_lines
  if (grepl("\\.gz$", file)) {
    conn <- gzfile(file, "r")
  } else {
    conn <- file(file, "r")
  }
  
  lines <- readLines(conn, n = n_lines)
  close(conn)
  
  if (length(lines) == 0) {
    stop("File appears to be empty")
  }
  
  # Parse header and data lines
  header_lines <- grep("^#", lines)
  data_lines <- setdiff(seq_along(lines), header_lines)
  
  if (length(data_lines) == 0) {
    stop("No data lines found in file")
  }
  
  # Parse first data line to check format
  first_data <- lines[data_lines[1]]
  fields <- strsplit(first_data, "\t")[[1]]
  
  # Check if it looks like ModKit bedMethyl format
  is_modkit <- length(fields) >= 11 && 
    grepl("^(chr)?[0-9XYMTxy]+$", fields[1]) &&  # chromosome
    grepl("^[0-9]+$", fields[2]) &&               # start position
    grepl("^[0-9]+$", fields[3])                  # end position
  
  # Extract sample of chromosomes and positions
  chr_sample <- character()
  pos_sample <- integer()
  mod_sample <- character()
  
  for (i in head(data_lines, min(100, length(data_lines)))) {
    line_fields <- strsplit(lines[i], "\t")[[1]]
    if (length(line_fields) >= 4) {
      chr_sample <- c(chr_sample, line_fields[1])
      pos_sample <- c(pos_sample, as.integer(line_fields[2]))
      mod_sample <- c(mod_sample, line_fields[4])
    }
  }
  
  result <- list(
    file = basename(file),
    has_tabix_index = has_index,
    total_lines_read = length(lines),
    header_lines = length(header_lines),
    data_lines = length(data_lines),
    appears_to_be_modkit = is_modkit,
    n_fields = length(fields),
    sample_chromosomes = unique(chr_sample),
    position_range = if (length(pos_sample) > 0) range(pos_sample) else c(NA, NA),
    modification_codes = unique(mod_sample),
    first_data_line = first_data
  )
  
  return(result)
}