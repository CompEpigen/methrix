#include "modkit_r_interface.h"

// Simple test function to verify registration works
SEXP test_modkit_c() {
  return ScalarString(mkChar("ModKit C interface is working"));
}

// Global result structure for accumulating data
static modkit_r_result_t *g_result = NULL;
static int g_current_site = 0;

// Main R interface function
SEXP read_modkit_c(SEXP files, SEXP target_mod, SEXP bin_size, 
                   SEXP n_cores, SEXP quality_filter, SEXP min_coverage,
                   SEXP combine_strands, SEXP ref_fasta, SEXP verbose) {
  
  // Parse R arguments into config structure
  config_t config = parse_r_config(files, target_mod, bin_size, n_cores,
                                   quality_filter, min_coverage, combine_strands,
                                   ref_fasta, verbose);
  
  // Initialize reference if provided
  if (config.ref_fasta) {
    if (init_reference(&config) != 0) {
      cleanup_config(&config);
      error("Failed to initialize reference FASTA: %s", config.ref_fasta);
    }
  }
  
  // Auto-discover chromosomes (keep existing implementation)
  if (discover_chromosomes(&config) != 0) {
    cleanup_reference();
    cleanup_config(&config);
    error("Failed to discover chromosomes from input files");
  }
  
  // Generate genomic bins (keep existing implementation)  
  if (generate_bins(&config) != 0) {
    cleanup_reference();
    cleanup_config(&config);
    error("Failed to generate genomic bins");
  }
  
  // Estimate total sites for memory allocation
  int estimated_sites = estimate_total_sites(&config);
  
  // Initialize result structure
  g_result = init_modkit_result(estimated_sites, config.n_files);
  if (!g_result) {
    cleanup_reference();
    cleanup_config(&config);
    error("Failed to allocate memory for results");
  }
  g_current_site = 0;
  
  // Process all bins (modified to accumulate in memory)
  int process_result = process_all_bins_to_memory(&config);
  
  if (process_result != 0) {
    cleanup_modkit_result(g_result);
    cleanup_reference();
    cleanup_config(&config);
    error("Failed to process genomic bins");
  }
  
  // Convert result to R objects
  SEXP r_result = convert_modkit_result_to_r(g_result);
  
  // Cleanup
  cleanup_modkit_result(g_result);
  cleanup_reference();
  cleanup_config(&config);
  g_result = NULL;
  
  return r_result;
}

// Parse R arguments into config structure
config_t parse_r_config(SEXP files, SEXP target_mod, SEXP bin_size, 
                        SEXP n_cores, SEXP quality_filter, SEXP min_coverage,
                        SEXP combine_strands, SEXP ref_fasta, SEXP verbose) {
  
  config_t config = {0};
  
  // Parse files
  config.n_files = length(files);
  config.input_files = (char**)R_alloc(config.n_files, sizeof(char*));
  for (int i = 0; i < config.n_files; i++) {
    const char *file_path = CHAR(STRING_ELT(files, i));
    config.input_files[i] = (char*)R_alloc(strlen(file_path) + 1, sizeof(char));
    strcpy(config.input_files[i], file_path);
  }
  
  // Parse other parameters
  config.target_mod = CHAR(STRING_ELT(target_mod, 0))[0];
  config.bin_size = asInteger(bin_size);
  config.n_threads = asInteger(n_cores);
  config.quality_filter = asReal(quality_filter);
  config.min_coverage = asInteger(min_coverage);
  config.combine_strands = asLogical(combine_strands);
  config.verbose = asLogical(verbose);
  
  // Handle optional reference FASTA
  if (!isNull(ref_fasta) && length(ref_fasta) > 0) {
    const char *ref_path = CHAR(STRING_ELT(ref_fasta, 0));
    config.ref_fasta = (char*)R_alloc(strlen(ref_path) + 1, sizeof(char));
    strcpy(config.ref_fasta, ref_path);
    config.include_context = 1;  // Enable context when ref provided
  } else {
    config.ref_fasta = NULL;
    config.include_context = 0;
  }
  
  // Set defaults
  config.region = NULL;
  config.output_file = NULL;
  config.validate_cytosines = 0;
  
  return config;
}

// Initialize result structure (UPDATED with context fields)
modkit_r_result_t* init_modkit_result(int n_sites, int n_samples) {
  modkit_r_result_t *result = (modkit_r_result_t*)R_alloc(1, sizeof(modkit_r_result_t));
  
  // Allocate coordinate vectors
  result->chr_names = (char**)R_alloc(n_sites, sizeof(char*));
  result->positions = (int32_t*)R_alloc(n_sites, sizeof(int32_t));
  result->strands = (char*)R_alloc(n_sites, sizeof(char));
  result->n_sites = n_sites;  // Set to allocated size initially
  
  // NEW: Initialize context fields
  result->include_context = 0;  // Will be set if reference is available
  result->ref_bases = NULL;
  result->contexts = NULL;
  
  // Allocate matrices (column-major for R)
  result->beta_matrix = (double*)R_alloc(n_sites * n_samples, sizeof(double));
  result->cov_matrix = (int32_t*)R_alloc(n_sites * n_samples, sizeof(int32_t));
  result->n_samples = n_samples;
  
  // Initialize matrices with NA/0
  for (int i = 0; i < n_sites * n_samples; i++) {
    result->beta_matrix[i] = R_NaN;
    result->cov_matrix[i] = 0;
  }
  
  // Allocate sample names
  result->sample_names = (char**)R_alloc(n_samples, sizeof(char*));
  
  // Initialize other fields
  result->chr_list = NULL;
  result->chr_counts = NULL;
  result->n_chromosomes = 0;
  result->total_sites_processed = 0;
  result->processing_notes = NULL;
  result->n_notes = 0;
  
  return result;
}

// Estimate total sites for memory allocation
int estimate_total_sites(config_t *config) {
  // Simple fixed size approach - start with reasonable buffer
  // Will grow dynamically if needed
  return 1000000;  // 1M sites initial allocation
}

// UPDATED: Modified function to accumulate results and extract context
void accumulate_coordinate_row(const coord_t *coord, file_stream_t *streams, int n_files,
                               config_t *config) {
  
  if (!g_result || g_current_site >= g_result->n_sites) {
    error("Result buffer overflow or not initialized");
  }
  
  // Store coordinate information
  int coord_len = strlen(coord->chr) + 1;
  g_result->chr_names[g_current_site] = (char*)R_alloc(coord_len, sizeof(char));
  strcpy(g_result->chr_names[g_current_site], coord->chr);
  g_result->positions[g_current_site] = coord->start;
  g_result->strands[g_current_site] = coord->strand;
  
  // NEW: Extract reference context if FASTA provided
  if (config->include_context && ref_fai) {
    // Initialize context arrays if this is the first site with context
    if (!g_result->include_context) {
      g_result->ref_bases = (char*)R_alloc(g_result->n_sites, sizeof(char));
      g_result->contexts = (char**)R_alloc(g_result->n_sites, sizeof(char*));
      g_result->include_context = 1;
      
      // Initialize previously processed sites with 'N' and "CNN"
      for (int i = 0; i < g_current_site; i++) {
        g_result->ref_bases[i] = 'N';
        g_result->contexts[i] = (char*)R_alloc(4, sizeof(char));
        strcpy(g_result->contexts[i], "CNN");
      }
    }
    
    // Get reference base
    char ref_base = get_ref_base_cached(coord->chr, coord->start);
    g_result->ref_bases[g_current_site] = ref_base;
    
    // Get sequence context
    char *context = get_context(coord->chr, coord->start, coord->strand);
    int ctx_len = strlen(context) + 1;
    g_result->contexts[g_current_site] = (char*)R_alloc(ctx_len, sizeof(char));
    strcpy(g_result->contexts[g_current_site], context);
    free(context);  // Free the context string returned by get_context
  } else if (g_result->include_context) {
    // If we have context arrays but no reference for this site
    g_result->ref_bases[g_current_site] = 'N';
    g_result->contexts[g_current_site] = (char*)R_alloc(4, sizeof(char));
    strcpy(g_result->contexts[g_current_site], "CNN");
  }
  
  // Store matrix data (column-major indexing for R)
  for (int i = 0; i < n_files; i++) {
    int matrix_idx = g_current_site + i * g_result->n_sites;  // Column-major
    
    if (streams[i].has_output_data) {
      g_result->beta_matrix[matrix_idx] = streams[i].output_beta;
      g_result->cov_matrix[matrix_idx] = streams[i].output_cov;
    } else {
      g_result->beta_matrix[matrix_idx] = R_NaN;
      g_result->cov_matrix[matrix_idx] = 0;
    }
  }
  
  g_current_site++;
}

// UPDATED: Convert C result to R objects (with context data)
SEXP convert_modkit_result_to_r(modkit_r_result_t *result) {
  
  // Adjust n_sites to actual processed sites
  int actual_sites = g_current_site;
  
  // Determine number of elements in result list
  int n_elements = result->include_context ? 9 : 7;
  
  // Create R list with components
  SEXP r_result = PROTECT(allocVector(VECSXP, n_elements));
  SEXP names = PROTECT(allocVector(STRSXP, n_elements));
  
  // 1. Chromosome vector
  SEXP chr_vec = PROTECT(allocVector(STRSXP, actual_sites));
  for (int i = 0; i < actual_sites; i++) {
    SET_STRING_ELT(chr_vec, i, mkChar(result->chr_names[i]));
  }
  SET_VECTOR_ELT(r_result, 0, chr_vec);
  SET_STRING_ELT(names, 0, mkChar("chr"));
  
  // 2. Position vector
  SEXP pos_vec = PROTECT(allocVector(INTSXP, actual_sites));
  memcpy(INTEGER(pos_vec), result->positions, actual_sites * sizeof(int32_t));
  SET_VECTOR_ELT(r_result, 1, pos_vec);
  SET_STRING_ELT(names, 1, mkChar("start"));
  
  // 3. Strand vector
  SEXP strand_vec = PROTECT(allocVector(STRSXP, actual_sites));
  for (int i = 0; i < actual_sites; i++) {
    char strand_str[2] = {result->strands[i], '\0'};
    SET_STRING_ELT(strand_vec, i, mkChar(strand_str));
  }
  SET_VECTOR_ELT(r_result, 2, strand_vec);
  SET_STRING_ELT(names, 2, mkChar("strand"));
  
  // 4. Beta matrix - only copy actual data
  SEXP beta_mat = PROTECT(allocMatrix(REALSXP, actual_sites, result->n_samples));
  for (int i = 0; i < actual_sites; i++) {
    for (int j = 0; j < result->n_samples; j++) {
      int src_idx = i + j * result->n_sites;  // Source column-major index
      int dst_idx = i + j * actual_sites;     // Destination column-major index
      REAL(beta_mat)[dst_idx] = result->beta_matrix[src_idx];
    }
  }
  SET_VECTOR_ELT(r_result, 3, beta_mat);
  SET_STRING_ELT(names, 3, mkChar("beta_matrix"));
  
  // 5. Coverage matrix - only copy actual data
  SEXP cov_mat = PROTECT(allocMatrix(INTSXP, actual_sites, result->n_samples));
  for (int i = 0; i < actual_sites; i++) {
    for (int j = 0; j < result->n_samples; j++) {
      int src_idx = i + j * result->n_sites;  // Source column-major index
      int dst_idx = i + j * actual_sites;     // Destination column-major index
      INTEGER(cov_mat)[dst_idx] = result->cov_matrix[src_idx];
    }
  }
  SET_VECTOR_ELT(r_result, 4, cov_mat);
  SET_STRING_ELT(names, 4, mkChar("cov_matrix"));
  
  // 6. Sample names
  SEXP sample_names = PROTECT(allocVector(STRSXP, result->n_samples));
  for (int i = 0; i < result->n_samples; i++) {
    if (result->sample_names[i]) {
      SET_STRING_ELT(sample_names, i, mkChar(result->sample_names[i]));
    }
  }
  SET_VECTOR_ELT(r_result, 5, sample_names);
  SET_STRING_ELT(names, 5, mkChar("sample_names"));
  
  // 7. Metadata
  SEXP metadata = PROTECT(allocVector(VECSXP, 2));
  SEXP meta_names = PROTECT(allocVector(STRSXP, 2));
  
  // Total sites processed
  SEXP total_sites = PROTECT(allocVector(INTSXP, 1));
  INTEGER(total_sites)[0] = actual_sites;
  SET_VECTOR_ELT(metadata, 0, total_sites);
  SET_STRING_ELT(meta_names, 0, mkChar("total_sites"));
  
  // Processing info
  SEXP proc_info = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(proc_info, 0, mkChar("ModKit C engine processing completed"));
  SET_VECTOR_ELT(metadata, 1, proc_info);
  SET_STRING_ELT(meta_names, 1, mkChar("processing_info"));
  
  setAttrib(metadata, R_NamesSymbol, meta_names);
  SET_VECTOR_ELT(r_result, 6, metadata);
  SET_STRING_ELT(names, 6, mkChar("metadata"));
  
  // NEW: Add reference base and context vectors if available
  int unprotect_count = 12;
  if (result->include_context) {
    // 8. Reference base vector
    SEXP ref_base_vec = PROTECT(allocVector(STRSXP, actual_sites));
    for (int i = 0; i < actual_sites; i++) {
      char base_str[2] = {result->ref_bases[i], '\0'};
      SET_STRING_ELT(ref_base_vec, i, mkChar(base_str));
    }
    SET_VECTOR_ELT(r_result, 7, ref_base_vec);
    SET_STRING_ELT(names, 7, mkChar("ref_base"));
    
    // 9. Context vector
    SEXP context_vec = PROTECT(allocVector(STRSXP, actual_sites));
    for (int i = 0; i < actual_sites; i++) {
      SET_STRING_ELT(context_vec, i, mkChar(result->contexts[i]));
    }
    SET_VECTOR_ELT(r_result, 8, context_vec);
    SET_STRING_ELT(names, 8, mkChar("context"));
    
    unprotect_count = 14;  // Added 2 more PROTECT calls
  }
  
  // Set names on result list
  setAttrib(r_result, R_NamesSymbol, names);
  
  UNPROTECT(unprotect_count);
  return r_result;
}

// Process all bins to memory (implementation moved from modkit_r_processing.c)
int process_all_bins_to_memory(config_t *config) {
  
  // Extract sample names from file paths
  for (int i = 0; i < config->n_files; i++) {
    const char *basename = strrchr(config->input_files[i], '/');
    basename = basename ? basename + 1 : config->input_files[i];
    
    char *sample_name = (char*)R_alloc(256, sizeof(char));
    strncpy(sample_name, basename, 255);
    sample_name[255] = '\0';
    
    // Remove file extensions
    char *dot = strchr(sample_name, '.');
    if (dot) *dot = '\0';
    
    g_result->sample_names[i] = sample_name;
  }
  
  // Process bins sequentially (parallel processing simplified for R compatibility)
  if (config->verbose) {
    Rprintf("Processing %d bins sequentially...\n", config->n_bins);
    if (config->include_context) {
      Rprintf("Reference FASTA provided - extracting sequence context...\n");
    }
  }
  
  for (int bin_idx = 0; bin_idx < config->n_bins; bin_idx++) {
    int result = process_single_bin_to_memory(&config->bins[bin_idx], 
                                              config->input_files,
                                              config->n_files, config);
    
    if (result != 0) {
      return result;
    }
    
    if (config->verbose && (bin_idx + 1) % 100 == 0) {
      Rprintf("Processed %d/%d bins\n", bin_idx + 1, config->n_bins);
    }
    
    // Check for user interruption
    R_CheckUserInterrupt();
  }
  
  return 0;
}

// Single bin processing to accumulate in memory
int process_single_bin_to_memory(genomic_bin_t *bin, char **input_files, int n_files,
                                 config_t *config) {
  
  file_stream_t *streams = (file_stream_t*)R_alloc(n_files, sizeof(file_stream_t));
  
  // Initialize streams for this bin
  int active_streams = 0;
  for (int i = 0; i < n_files; i++) {
    memset(&streams[i], 0, sizeof(file_stream_t));
    if (init_file_stream_for_bin(&streams[i], input_files[i], i, bin, config) == 0) {
      if (streams[i].has_current_record) {
        active_streams++;
      }
    }
  }
  
  if (active_streams == 0) {
    cleanup_streams(streams, n_files);
    return 0;  // No data in this bin
  }
  
  // Process data in this bin using streaming merge
  while (any_stream_has_data(streams, n_files)) {
    coord_t min_coord = find_minimum_coordinate(streams, n_files);
    collect_coordinate_data(streams, n_files, &min_coord, config);
    
    // Check if we need more space BEFORE accumulating
    if (g_current_site >= g_result->n_sites) {
      reallocate_result_if_needed();
    }
    
    // Modified: accumulate instead of writing to file
    accumulate_coordinate_row(&min_coord, streams, n_files, config);
    
    reset_output_data(streams, n_files);
  }
  
  cleanup_streams(streams, n_files);
  return 0;
}

// UPDATED: Reallocate result structure if we underestimated the size (with context arrays)
void reallocate_result_if_needed() {
  if (!g_result || g_current_site < g_result->n_sites) {
    return;  // No reallocation needed
  }
  
  // Increase size by 50%
  int new_size = (int)(g_result->n_sites * 1.5);
  
  if (g_result->n_samples <= 0) {
    error("Invalid number of samples during reallocation");
  }
  
  // Reallocate coordinate arrays
  char **new_chr_names = (char**)R_alloc(new_size, sizeof(char*));
  int32_t *new_positions = (int32_t*)R_alloc(new_size, sizeof(int32_t));
  char *new_strands = (char*)R_alloc(new_size, sizeof(char));
  
  // NEW: Reallocate context arrays if needed
  char *new_ref_bases = NULL;
  char **new_contexts = NULL;
  if (g_result->include_context) {
    new_ref_bases = (char*)R_alloc(new_size, sizeof(char));
    new_contexts = (char**)R_alloc(new_size, sizeof(char*));
  }
  
  // Reallocate matrices
  double *new_beta_matrix = (double*)R_alloc(new_size * g_result->n_samples, sizeof(double));
  int32_t *new_cov_matrix = (int32_t*)R_alloc(new_size * g_result->n_samples, sizeof(int32_t));
  
  // Copy existing data
  memcpy(new_chr_names, g_result->chr_names, g_result->n_sites * sizeof(char*));
  memcpy(new_positions, g_result->positions, g_result->n_sites * sizeof(int32_t));
  memcpy(new_strands, g_result->strands, g_result->n_sites * sizeof(char));
  
  if (g_result->include_context) {
    memcpy(new_ref_bases, g_result->ref_bases, g_result->n_sites * sizeof(char));
    memcpy(new_contexts, g_result->contexts, g_result->n_sites * sizeof(char*));
  }
  
  memcpy(new_beta_matrix, g_result->beta_matrix, 
         g_result->n_sites * g_result->n_samples * sizeof(double));
  memcpy(new_cov_matrix, g_result->cov_matrix,
         g_result->n_sites * g_result->n_samples * sizeof(int32_t));
  
  // Initialize new space with NA/0
  for (int i = g_result->n_sites * g_result->n_samples; 
       i < new_size * g_result->n_samples; i++) {
    new_beta_matrix[i] = R_NaN;
    new_cov_matrix[i] = 0;
  }
  
  // Update result structure
  g_result->chr_names = new_chr_names;
  g_result->positions = new_positions;
  g_result->strands = new_strands;
  g_result->beta_matrix = new_beta_matrix;
  g_result->cov_matrix = new_cov_matrix;
  g_result->n_sites = new_size;
  
  if (g_result->include_context) {
    g_result->ref_bases = new_ref_bases;
    g_result->contexts = new_contexts;
  }
  
  Rprintf("Reallocated result arrays to %d sites\n", new_size);
}

// Cleanup function
void cleanup_modkit_result(modkit_r_result_t *result) {
  // R handles memory cleanup automatically for R_alloc'd memory
  // This function is mainly for completeness
  if (result) {
    result->n_sites = 0;
    result->n_samples = 0;
    result->include_context = 0;
  }
}