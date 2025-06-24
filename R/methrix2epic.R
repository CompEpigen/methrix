# Helper function to detect genome build from methrix metadata
detect_genome_build <- function(m) {
  genome_info <- S4Vectors::metadata(m)$genome
  
  if (is.null(genome_info) || is.na(genome_info)) {
    return(NULL)
  }
  
  genome_info <- as.character(genome_info)
  
  if (grepl("hg19|GRCh37|Hs37|hs37|b37", genome_info, ignore.case = TRUE)) {
    return("hg19")
  } else if (grepl("hg38|GRCh38|Hs38|hs38|b38", genome_info, ignore.case = TRUE)) {
    return("hg38")
  } else {
    return(NULL)
  }
}

# Helper function to get array probe annotations
get_array_annotations <- function(array_type, genome) {
  
  all_annotations <- data.table::data.table()
  
  # Get 450K annotations
  if ("450K" %in% array_type) {
    tryCatch({
      if (genome == "hg19") {
        if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
          stop("Package 'IlluminaHumanMethylation450kanno.ilmn12.hg19' is required. Please install it using:\n",
               "BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19')")
        }
        
        # Load annotation object directly
        utils::data("Locations", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19", envir = environment())
        ann450k <- Locations
        
      } else if (genome == "hg38") {
        if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg38", quietly = TRUE)) {
          stop("Package 'IlluminaHumanMethylation450kanno.ilmn12.hg38' is required. Please install it using:\n",
               "BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg38')")
        }
        
        # Load annotation object directly
        utils::data("Locations", package = "IlluminaHumanMethylation450kanno.ilmn12.hg38", envir = environment())
        ann450k <- Locations
      }
      
      # Convert to data.table
      ann450k_dt <- data.table::as.data.table(as.data.frame(ann450k))
      ann450k_dt[, `:=`(
        probe_id = rownames(ann450k),
        array_type = "450K"
      )]
      
      # Select and rename columns
      ann450k_dt <- ann450k_dt[, .(
        chr = as.character(chr),
        start = as.integer(pos),
        end = as.integer(pos),
        probe_id,
        array_type
      )]
      
      all_annotations <- rbind(all_annotations, ann450k_dt)
      
    }, error = function(e) {
      warning("Could not load 450K annotations: ", e$message)
    })
  }
  
  # Get EPIC annotations
  if ("EPIC" %in% array_type) {
    tryCatch({
      if (genome == "hg19") {
        if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
          stop("Package 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19' is required. Please install it using:\n",
               "BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b4.hg19')")
        }
        
        # Load annotation object directly
        utils::data("Locations", package = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", envir = environment())
        annEPIC <- Locations
        
      } else if (genome == "hg38") {
        if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b5.hg38", quietly = TRUE)) {
          stop("Package 'IlluminaHumanMethylationEPICanno.ilm10b5.hg38' is required. Please install it using:\n",
               "BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b5.hg38')")
        }
        
        # Load annotation object directly
        utils::data("Locations", package = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38", envir = environment())
        annEPIC <- Locations
      }
      
      # Convert to data.table
      annEPIC_dt <- data.table::as.data.table(as.data.frame(annEPIC))
      annEPIC_dt[, `:=`(
        probe_id = rownames(annEPIC),
        array_type = "EPIC"
      )]
      
      # Select and rename columns
      annEPIC_dt <- annEPIC_dt[, .(
        chr = as.character(chr),
        start = as.integer(pos),
        end = as.integer(pos),
        probe_id,
        array_type
      )]
      
      all_annotations <- rbind(all_annotations, annEPIC_dt)
      
    }, error = function(e) {
      warning("Could not load EPIC annotations: ", e$message)
    })
  }
  
  if (nrow(all_annotations) == 0) {
    stop("No array annotations could be loaded\nPlease install IlluminaHumanMethylationEPICanno.ilm10b4.hg19 or IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
  }
  
  # Remove duplicates (probes present in both arrays)
  all_annotations <- unique(all_annotations, by = c("chr", "start"))
  
  # Handle chromosome naming consistency
  if (!any(grepl("^chr", all_annotations$chr))) {
    all_annotations[, chr := paste0("chr", chr)]
  }
  
  # Set key for efficient joining
  data.table::setkey(all_annotations, chr, start)
  
  return(all_annotations)
}

# Helper function for exact coordinate matching
match_exact_coordinates <- function(methrix_coords, array_annotations) {
  
  # Ensure consistent chromosome naming
  methrix_coords <- data.table::copy(methrix_coords)
  if (!any(grepl("^chr", methrix_coords$chr))) {
    methrix_coords[, chr := paste0("chr", chr)]
  }
  
  # Set keys for efficient joining
  data.table::setkey(methrix_coords, chr, start)
  data.table::setkey(array_annotations, chr, start)
  
  # Find exact matches
  matches <- array_annotations[methrix_coords, which = TRUE, nomatch = 0]
  
  # Return indices of matching methrix coordinates
  matched_indices <- which(!is.na(matches))
  
  return(list(
    methrix_indices = matched_indices,
    array_indices = matches[matched_indices],
    match_info = array_annotations[matches[matched_indices]]
  ))
}

# Helper function for tolerance-based coordinate matching
match_tolerance_coordinates <- function(methrix_coords, array_annotations, tolerance) {
  
  # Convert to GRanges for overlap-based matching
  methrix_coords_copy <- data.table::copy(methrix_coords)
  if (!any(grepl("^chr", methrix_coords_copy$chr))) {
    methrix_coords_copy[, chr := paste0("chr", chr)]
  }
  
  # Create GRanges objects
  methrix_gr <- GenomicRanges::makeGRangesFromDataFrame(
    methrix_coords_copy,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "start"  # Point ranges for CpGs
  )
  
  # Expand array coordinates by tolerance
  array_coords_expanded <- data.table::copy(array_annotations)
  array_coords_expanded[, `:=`(
    start = start - tolerance,
    end = end + tolerance
  )]
  
  array_gr <- GenomicRanges::makeGRangesFromDataFrame(
    array_coords_expanded,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  # Find overlaps
  overlaps <- GenomicRanges::findOverlaps(methrix_gr, array_gr)
  
  # Extract matching indices
  methrix_indices <- S4Vectors::queryHits(overlaps)
  array_indices <- S4Vectors::subjectHits(overlaps)
  
  # Remove duplicates (keep first match for each methrix coordinate)
  unique_methrix <- !duplicated(methrix_indices)
  methrix_indices <- methrix_indices[unique_methrix]
  array_indices <- array_indices[unique_methrix]
  
  return(list(
    methrix_indices = methrix_indices,
    array_indices = array_indices,
    match_info = array_annotations[array_indices]
  ))
}

# Helper function to validate inputs
validate_inputs <- function(m, array_type, genome, match_tolerance, return_eset) {
  
  # Check methrix object
  if (!is(m, "methrix")) {
    stop("Input 'm' must be a methrix object")
  }
  
  # Check array_type
  valid_array_types <- c("450K", "EPIC", "both")
  if (!array_type %in% valid_array_types) {
    stop("array_type must be one of: ", paste(valid_array_types, collapse = ", "))
  }
  
  # Check genome if provided
  if (!is.null(genome)) {
    valid_genomes <- c("hg19", "hg38")
    if (!genome %in% valid_genomes) {
      stop("genome must be one of: ", paste(valid_genomes, collapse = ", "))
    }
  }
  
  # Check match_tolerance
  if (!is.numeric(match_tolerance) || match_tolerance < 0 || match_tolerance != as.integer(match_tolerance)) {
    stop("match_tolerance must be a non-negative integer")
  }
  
  # Check return_eset
  if (!is.logical(return_eset) || length(return_eset) != 1) {
    stop("return_eset must be a single logical value (TRUE or FALSE)")
  }
  
  # Check if methrix has any CpGs
  if (nrow(m) == 0) {
    stop("Input methrix object has no CpG sites")
  }
  
  return(TRUE)
}

# Main function

#' Subset methrix object to Illumina array probes
#'
#' @description Subsets a methrix object to include only CpG sites that 
#' correspond to Illumina 450K or EPIC array probes, enabling direct comparison 
#' with array-based methylation data.
#' 
#' @details This function matches CpG coordinates from a methrix object with 
#' probe coordinates from Illumina methylation arrays. It supports both 450K 
#' and EPIC arrays across hg19 and hg38 genome builds. The function can 
#' auto-detect the genome build from methrix metadata or use a user-specified build.
#' 
#' @param m A methrix object
#' @param array_type Character vector specifying which array probes to include. 
#' Options: "450K", "EPIC", or "both". Default: c("450K", "EPIC", "both")
#' @param genome Character specifying genome build ("hg19" or "hg38"). 
#' If NULL (default), auto-detected from methrix metadata
#' @param match_tolerance Integer. Tolerance in base pairs for coordinate matching. 
#' Default: 0 (exact matching)
#' @param keep_metadata Logical. Whether to add array probe metadata to output rowData. 
#' Default: TRUE
#' @param return_eset Logical. Whether to return an ExpressionSet object instead of methrix. 
#' Default: FALSE
#' @param verbose Logical. Whether to print progress messages. Default: TRUE
#' 
#' @return A subsetted methrix object or ExpressionSet (if return_eset=TRUE) containing 
#' only array probe sites. If keep_metadata=TRUE, rowData/featureData will include 
#' probe_id and array_type columns.
#' 
#' @examples
#' \dontrun{
#' data('methrix_data')
#' 
#' # Subset to EPIC array probes (auto-detect genome)
#' epic_data <- methrix2epic(methrix_data, array_type = "EPIC")
#' 
#' # Subset to 450K probes with specific genome
#' k450_data <- methrix2epic(methrix_data, array_type = "450K", genome = "hg19")
#' 
#' # Subset to both arrays with 2bp tolerance
#' both_data <- methrix2epic(methrix_data, array_type = "both", match_tolerance = 2)
#' 
#' # Return as ExpressionSet instead of methrix
#' epic_eset <- methrix2epic(methrix_data, array_type = "EPIC", return_eset = TRUE)
#' 
#' # Check what was subsetted
#' S4Vectors::metadata(epic_data)$array_subset
#' 
#' # Check ExpressionSet results  
#' summarize_array_subset(epic_eset)
#' }
#' 
#' @export
#' @importFrom data.table data.table setkey foverlaps as.data.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors metadata queryHits subjectHits
#' @importFrom utils data
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame assayDataNew
#' 
methrix2epic <- function(m, array_type = c("450K", "EPIC", "both"), 
                         genome = NULL, match_tolerance = 0, 
                         keep_metadata = TRUE, return_eset = FALSE, verbose = TRUE) {
  
  start_time <- proc.time()
  
  # Match and validate array_type argument
  array_type <- match.arg(array_type, choices = c("450K", "EPIC", "both"))
  if (array_type == "both") {
    array_type <- c("450K", "EPIC")
  }
  
  # Validate all inputs
  validate_inputs(m, array_type, genome, match_tolerance, return_eset)
  
  # Detect or validate genome build
  if (is.null(genome)) {
    genome <- detect_genome_build(m)
    if (is.null(genome)) {
      stop("Cannot determine genome build from methrix metadata. Please specify 'hg19' or 'hg38'")
    }
    if (verbose) message("Auto-detected genome build: ", genome)
  } else {
    genome <- match.arg(genome, choices = c("hg19", "hg38"))
    if (verbose) message("Using specified genome build: ", genome)
  }
  
  # Get array annotations
  if (verbose) message("Loading array annotations for ", paste(array_type, collapse = " and "), " arrays...")
  array_annotations <- get_array_annotations(array_type, genome)
  
  if (verbose) message("Found ", nrow(array_annotations), " unique array probe coordinates")
  
  # Get methrix coordinates
  methrix_coords <- data.table::as.data.table(rowData(m))
  
  # Find matching coordinates
  if (verbose) message("Matching coordinates", 
                       if (match_tolerance > 0) paste(" with", match_tolerance, "bp tolerance") else " (exact matching)", "...")
  
  if (match_tolerance == 0) {
    match_results <- match_exact_coordinates(methrix_coords, array_annotations)
  } else {
    match_results <- match_tolerance_coordinates(methrix_coords, array_annotations, match_tolerance)
  }
  
  matched_indices <- match_results$methrix_indices
  
  # Check if any matches were found
  if (length(matched_indices) == 0) {
    stop("No matching probes found between methrix and array annotations.\n",
         "This might be due to:\n",
         "- Incorrect genome build\n",
         "- Different chromosome naming conventions\n",
         "- No overlap between your data and array probe locations")
  }
  
  # Subset methrix object
  m_subset <- m[matched_indices, ]
  
  # Set rownames to probe IDs for the methrix object
  probe_ids <- match_results$match_info$probe_id
  rownames(m_subset) <- probe_ids
  
  # Add array metadata if requested
  if (keep_metadata) {
    probe_metadata <- match_results$match_info[, .(probe_id, array_type)]
    
    # Combine with existing rowData
    new_rowData <- cbind(
      data.table::as.data.table(rowData(m_subset)),
      probe_metadata
    )
    
    rowData(m_subset) <- new_rowData
  }
  
  # Update metadata
  array_coverage <- round(length(matched_indices) / nrow(array_annotations) * 100, 2)
  methrix_coverage <- round(length(matched_indices) / nrow(m) * 100, 2)
  
  S4Vectors::metadata(m_subset)$array_subset <- list(
    array_type = array_type,
    genome = genome,
    n_probes_found = length(matched_indices),
    n_probes_total_array = nrow(array_annotations),
    n_probes_total_methrix = nrow(m),
    array_coverage_percent = array_coverage,
    methrix_coverage_percent = methrix_coverage,
    match_tolerance = match_tolerance,
    processing_time = proc.time() - start_time
  )
  
  if (verbose) {
    message("Successfully subsetted to ", length(matched_indices), " array probes")
    message("Array coverage: ", array_coverage, "% (", length(matched_indices), "/", nrow(array_annotations), ")")
    message("Methrix coverage: ", methrix_coverage, "% (", length(matched_indices), "/", nrow(m), ")")
    message("Processing completed in ", round((proc.time() - start_time)[3], 2), " seconds")
  }
  
  # Convert to ExpressionSet if requested
  if (return_eset) {
    if (!requireNamespace("Biobase", quietly = TRUE)) {
      stop("Package 'Biobase' is required for ExpressionSet output. Please install it using:\n",
           "BiocManager::install('Biobase')")
    }
    
    if (verbose) message("Converting to ExpressionSet...")
    
    # Extract methylation matrix (beta values) and ensure it's a clean matrix
    exprs_matrix <- get_matrix(m_subset, type = "M", add_loci = FALSE)
    
    # Convert to clean matrix and remove any special attributes
    exprs_matrix <- as.matrix(exprs_matrix)
    class(exprs_matrix) <- "matrix"
    
    # Ensure rownames and colnames are character vectors (probe IDs should already be set)
    rownames(exprs_matrix) <- as.character(rownames(m_subset))
    colnames(exprs_matrix) <- as.character(colnames(m_subset))
    
    # Create phenoData (sample metadata)
    pheno_df <- as.data.frame(colData(m_subset))
    pheno_data <- Biobase::AnnotatedDataFrame(data = pheno_df)
    
    # Create featureData (probe/CpG metadata)
    feature_df <- as.data.frame(rowData(m_subset))
    # Make sure probe IDs are in the featureData rownames too
    rownames(feature_df) <- rownames(m_subset)
    feature_data <- Biobase::AnnotatedDataFrame(data = feature_df)
    
    # Create ExpressionSet using assayData parameter
    assay_data <- Biobase::assayDataNew("environment", exprs = exprs_matrix)
    eset <- Biobase::ExpressionSet(
      assayData = assay_data,
      phenoData = pheno_data,
      featureData = feature_data
    )
    
    # Add array subset metadata to experimentData
    array_subset_info <- S4Vectors::metadata(m_subset)$array_subset
    Biobase::experimentData(eset)@other$array_subset <- array_subset_info
    Biobase::experimentData(eset)@other$original_genome <- S4Vectors::metadata(m_subset)$genome
    
    return(eset)
  }
  
  return(m_subset)
}

# Additional utility function to summarize array subset results
summarize_array_subset <- function(obj) {
  
  # Handle both methrix and ExpressionSet objects
  if (is(obj, "methrix")) {
    array_info <- S4Vectors::metadata(obj)$array_subset
    obj_type <- "methrix"
  } else if (is(obj, "ExpressionSet")) {
    if (!requireNamespace("Biobase", quietly = TRUE)) {
      stop("Package 'Biobase' is required to summarize ExpressionSet objects")
    }
    array_info <- Biobase::experimentData(obj)@other$array_subset
    obj_type <- "ExpressionSet"
  } else {
    stop("Input must be a methrix object or ExpressionSet created using methrix2epic()")
  }
  
  if (is.null(array_info)) {
    stop("This object was not created using methrix2epic()")
  }
  
  cat("Array Subset Summary (", obj_type, " object)\n", sep = "")
  cat("==========================================\n")
  cat("Array type(s):", paste(array_info$array_type, collapse = ", "), "\n")
  cat("Genome build:", array_info$genome, "\n")
  cat("Match tolerance:", array_info$match_tolerance, "bp\n")
  cat("Probes found:", array_info$n_probes_found, "\n")
  cat("Array coverage:", array_info$array_coverage_percent, "%\n")
  cat("Original methrix coverage:", array_info$methrix_coverage_percent, "%\n")
  cat("Processing time:", round(array_info$processing_time[3], 2), "seconds\n")
  
  invisible(array_info)
}