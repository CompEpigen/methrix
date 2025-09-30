#include "modkit_merge_v2.h"

// ============================================================================
// R INTERFACE
// ============================================================================

SEXP read_modkit_v2_c(SEXP files, SEXP chrom_sizes_file, SEXP target_mod,
                      SEXP interval_size, SEXP min_coverage,
                      SEXP quality_filter, SEXP combine_strands, SEXP ref_fasta,
                      SEXP verbose) {

    // Parse R arguments
    int n_files = length(files);
    char **input_files = (char**)R_alloc(n_files, sizeof(char*));
    for (int i = 0; i < n_files; i++) {
        input_files[i] = (char*)CHAR(STRING_ELT(files, i));
    }

    const char *chrom_sizes_path = CHAR(STRING_ELT(chrom_sizes_file, 0));
    const char *ref_fasta_path = CHAR(STRING_ELT(ref_fasta, 0));
    int has_ref_fasta = (strlen(ref_fasta_path) > 0);

    discovery_config_t config;
    config.target_mod = CHAR(STRING_ELT(target_mod, 0))[0];
    config.min_coverage = asInteger(min_coverage);
    config.quality_filter = asReal(quality_filter);
    config.combine_strands = asLogical(combine_strands);
    config.verbose = asLogical(verbose);

    int32_t interval_bp = asInteger(interval_size);

    if (config.verbose) {
        Rprintf("==============================================\n");
        Rprintf("ModKit Merge V2: High-Performance Implementation\n");
        Rprintf("==============================================\n");
        Rprintf("Files: %d\n", n_files);
        Rprintf("Target modification: %c\n", config.target_mod);
        Rprintf("Chromosome sizes: %s\n", chrom_sizes_path);
        if (has_ref_fasta) {
            Rprintf("Reference FASTA: %s\n", ref_fasta_path);
        }
        Rprintf("Interval size: %d bp\n", interval_bp);
        Rprintf("\n");
    }

    // ========================================================================
    // PHASE 0: READ AND FILTER CHROMOSOME SIZES
    // ========================================================================

    chrom_sizes_t *all_sizes = read_chrom_sizes(chrom_sizes_path);
    if (!all_sizes) {
        error("Failed to read chromosome sizes file");
    }

    if (config.verbose) {
        Rprintf("Chromosome sizes file: %d chromosome(s)\n", all_sizes->n_chroms);
    }

    // Filter to only chromosomes present in tabix files
    chrom_sizes_t *filtered_sizes = filter_chrom_sizes_by_tabix(
        all_sizes, input_files, n_files, config.verbose
    );

    free_chrom_sizes(all_sizes);

    if (!filtered_sizes) {
        error("No chromosomes match between sizes file and input files");
    }

    // ========================================================================
    // PHASE 1: COORDINATE DISCOVERY
    // ========================================================================

    coordinate_set_t *coord_set = discover_coordinates(
        input_files, n_files, filtered_sizes, &config
    );
    if (!coord_set) {
        free_chrom_sizes(filtered_sizes);
        error("Failed to discover coordinates");
    }

    if (coord_set->n_sites == 0) {
        free_coordinate_set(coord_set);
        free_chrom_sizes(filtered_sizes);
        error("No sites discovered. Check target modification code and filters.");
    }

    // ========================================================================
    // PHASE 2: CHUNKING
    // ========================================================================

    if (config.verbose) {
        Rprintf("\nPhase 2: Creating Chunks\n");
    }

    chunk_list_t *chunks = create_chunks_by_interval(
        coord_set, filtered_sizes, interval_bp
    );
    if (!chunks) {
        free_coordinate_set(coord_set);
        free_chrom_sizes(filtered_sizes);
        error("Failed to create chunks");
    }

    if (config.verbose) {
        Rprintf("  Created %d chunks\n", chunks->n_chunks);
        Rprintf("  Average sites per chunk: %.0f\n",
                (double)coord_set->n_sites / chunks->n_chunks);
        Rprintf("\n");
    }

    free_chrom_sizes(filtered_sizes);  // No longer needed after chunking

    // ========================================================================
    // PHASE 3: CHUNK PROCESSING
    // ========================================================================

    if (config.verbose) {
        Rprintf("Phase 3: Processing Chunks\n");
    }

    // Allocate final R matrices
    SEXP beta_matrix = PROTECT(allocMatrix(REALSXP, coord_set->n_sites, n_files));
    SEXP cov_matrix = PROTECT(allocMatrix(INTSXP, coord_set->n_sites, n_files));

    double *beta_ptr = REAL(beta_matrix);
    int *cov_ptr = INTEGER(cov_matrix);

    // Initialize with NA/0
    for (size_t i = 0; i < (size_t)coord_set->n_sites * n_files; i++) {
        beta_ptr[i] = R_NaN;
        cov_ptr[i] = 0;
    }

    // Process each chunk
    for (int c = 0; c < chunks->n_chunks; c++) {
        chunk_result_t *result = process_chunk(&chunks->chunks[c], coord_set,
                                              input_files, n_files, &config);

        if (!result) {
            UNPROTECT(2);
            free_chunk_list(chunks);
            free_coordinate_set(coord_set);
            error("Failed to process chunk %d", c);
        }

        // Copy chunk result to final matrices (convert to R column-major layout)
        int32_t chunk_start = chunks->chunks[c].start_site_idx;
        for (int32_t local_idx = 0; local_idx < result->n_sites; local_idx++) {
            int32_t global_idx = chunk_start + local_idx;

            for (int s = 0; s < n_files; s++) {
                // Source: row-major (local_idx * n_files + s)
                size_t src_idx = (size_t)local_idx * n_files + s;

                // Destination: column-major (global_idx + s * n_total_sites)
                size_t dst_idx = global_idx + (size_t)s * coord_set->n_sites;

                // Only copy if we have data
                if (get_bit(result->has_data_bits, src_idx)) {
                    beta_ptr[dst_idx] = result->beta_values[src_idx];
                    cov_ptr[dst_idx] = result->cov_values[src_idx];
                }
            }
        }

        free_chunk_result(result);

        // Progress update
        if (config.verbose && (c + 1) % 10 == 0) {
            Rprintf("  Completed %d/%d chunks (%.1f%%)\n",
                    c + 1, chunks->n_chunks,
                    100.0 * (c + 1) / chunks->n_chunks);
        }
    }

    if (config.verbose) {
        Rprintf("  Completed all %d chunks\n\n", chunks->n_chunks);
    }

    // ========================================================================
    // PHASE 4: CREATE R RESULT
    // ========================================================================

    if (config.verbose) {
        Rprintf("Phase 4: Creating R Result\n");
    }

    // Create coordinate vectors
    SEXP chr_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));
    SEXP pos_vec = PROTECT(allocVector(INTSXP, coord_set->n_sites));
    SEXP strand_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));

    for (size_t i = 0; i < coord_set->n_sites; i++) {
        SET_STRING_ELT(chr_vec, i, mkChar(coord_set->sites[i].chr));
        INTEGER(pos_vec)[i] = coord_set->sites[i].pos;

        char strand_str[2] = {coord_set->sites[i].strand, '\0'};
        SET_STRING_ELT(strand_vec, i, mkChar(strand_str));
    }

    if (config.verbose) {
        Rprintf("  Total sites: %zu\n", coord_set->n_sites);
        Rprintf("  Total samples: %d\n\n", n_files);
    }

    // ========================================================================
    // PHASE 5: CONTEXT EXTRACTION (if ref_fasta provided)
    // ========================================================================

    SEXP ref_base_vec = R_NilValue;
    SEXP context_vec = R_NilValue;

    if (has_ref_fasta) {
        // Initialize reference FASTA
        if (init_reference_v2(ref_fasta_path) != 0) {
            Rprintf("Warning: Failed to initialize reference FASTA\n");
            Rprintf("Continuing without context extraction\n");
        } else {
            // Extract context for all sites
            char *ref_bases = NULL;
            char *contexts = NULL;

            if (extract_context_for_sites(coord_set, &ref_bases, &contexts, config.verbose) == 0) {
                // Create R vectors for context
                ref_base_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));
                context_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));

                for (size_t i = 0; i < coord_set->n_sites; i++) {
                    char ref_str[2] = {ref_bases[i], '\0'};
                    SET_STRING_ELT(ref_base_vec, i, mkChar(ref_str));
                    SET_STRING_ELT(context_vec, i, mkChar(&contexts[i * 4]));
                }

                free(ref_bases);
                free(contexts);
            }

            cleanup_reference_v2();
        }
    }

    // Create result list (with or without context)
    int n_elements = has_ref_fasta && (ref_base_vec != R_NilValue) ? 7 : 5;
    SEXP result_list = PROTECT(allocVector(VECSXP, n_elements));
    SET_VECTOR_ELT(result_list, 0, chr_vec);
    SET_VECTOR_ELT(result_list, 1, pos_vec);
    SET_VECTOR_ELT(result_list, 2, strand_vec);
    SET_VECTOR_ELT(result_list, 3, beta_matrix);
    SET_VECTOR_ELT(result_list, 4, cov_matrix);

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, n_elements));
    SET_STRING_ELT(names, 0, mkChar("chr"));
    SET_STRING_ELT(names, 1, mkChar("start"));
    SET_STRING_ELT(names, 2, mkChar("strand"));
    SET_STRING_ELT(names, 3, mkChar("beta_matrix"));
    SET_STRING_ELT(names, 4, mkChar("cov_matrix"));

    if (n_elements == 7) {
        SET_VECTOR_ELT(result_list, 5, ref_base_vec);
        SET_VECTOR_ELT(result_list, 6, context_vec);
        SET_STRING_ELT(names, 5, mkChar("ref_base"));
        SET_STRING_ELT(names, 6, mkChar("context"));
    }

    setAttrib(result_list, R_NamesSymbol, names);

    // Print final stats
    if (config.verbose) {
        Rprintf("Done!\n");
        Rprintf("==============================================\n");
    }

    // Cleanup
    free_chunk_list(chunks);
    free_coordinate_set(coord_set);

    int n_protect = (n_elements == 7) ? 9 : 7;
    UNPROTECT(n_protect);
    return result_list;
}

// ============================================================================
// STREAMING HDF5 SUPPORT: CHUNK METADATA
// ============================================================================

SEXP read_modkit_v2_init_chunks_c(SEXP files, SEXP chrom_sizes_file, SEXP target_mod,
                                   SEXP interval_size, SEXP min_coverage,
                                   SEXP quality_filter, SEXP combine_strands, SEXP ref_fasta,
                                   SEXP verbose) {

    // Parse R arguments (same as main function)
    int n_files = length(files);
    char **input_files = (char**)R_alloc(n_files, sizeof(char*));
    for (int i = 0; i < n_files; i++) {
        input_files[i] = (char*)CHAR(STRING_ELT(files, i));
    }

    const char *chrom_sizes_path = CHAR(STRING_ELT(chrom_sizes_file, 0));
    const char *ref_fasta_path = CHAR(STRING_ELT(ref_fasta, 0));
    int has_ref_fasta = (strlen(ref_fasta_path) > 0);

    discovery_config_t config;
    config.target_mod = CHAR(STRING_ELT(target_mod, 0))[0];
    config.min_coverage = asInteger(min_coverage);
    config.quality_filter = asReal(quality_filter);
    config.combine_strands = asLogical(combine_strands);
    config.verbose = asLogical(verbose);

    int32_t interval_bp = asInteger(interval_size);

    if (config.verbose) {
        Rprintf("==============================================\n");
        Rprintf("ModKit Merge V2: Streaming Mode Initialization\n");
        Rprintf("==============================================\n");
    }

    // PHASE 0: Read chromosome sizes
    chrom_sizes_t *all_sizes = read_chrom_sizes(chrom_sizes_path);
    if (!all_sizes) {
        error("Failed to read chromosome sizes file");
    }

    chrom_sizes_t *filtered_sizes = filter_chrom_sizes_by_tabix(
        all_sizes, input_files, n_files, config.verbose
    );
    free_chrom_sizes(all_sizes);

    if (!filtered_sizes) {
        error("No chromosomes match between sizes file and input files");
    }

    // PHASE 1: Coordinate discovery
    coordinate_set_t *coord_set = discover_coordinates(
        input_files, n_files, filtered_sizes, &config
    );
    if (!coord_set) {
        free_chrom_sizes(filtered_sizes);
        error("Failed to discover coordinates");
    }

    if (coord_set->n_sites == 0) {
        free_coordinate_set(coord_set);
        free_chrom_sizes(filtered_sizes);
        error("No sites discovered");
    }

    // PHASE 2: Create chunks
    chunk_list_t *chunks = create_chunks_by_interval(
        coord_set, filtered_sizes, interval_bp
    );
    if (!chunks) {
        free_coordinate_set(coord_set);
        free_chrom_sizes(filtered_sizes);
        error("Failed to create chunks");
    }

    free_chrom_sizes(filtered_sizes);

    if (config.verbose) {
        Rprintf("Created %d chunks for streaming processing\n", chunks->n_chunks);
        Rprintf("==============================================\n\n");
    }

    // Store coord_set and chunks as external pointers
    SEXP coord_ptr = PROTECT(R_MakeExternalPtr(coord_set, R_NilValue, R_NilValue));
    SEXP chunks_ptr = PROTECT(R_MakeExternalPtr(chunks, R_NilValue, R_NilValue));

    // Extract context if ref_fasta provided
    SEXP ref_base_vec = R_NilValue;
    SEXP context_vec = R_NilValue;

    if (has_ref_fasta) {
        if (init_reference_v2(ref_fasta_path) != 0) {
            Rprintf("Warning: Failed to initialize reference FASTA\n");
        } else {
            char *ref_bases = NULL;
            char *contexts = NULL;

            if (extract_context_for_sites(coord_set, &ref_bases, &contexts, config.verbose) == 0) {
                ref_base_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));
                context_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));

                for (size_t i = 0; i < coord_set->n_sites; i++) {
                    char ref_str[2] = {ref_bases[i], '\0'};
                    SET_STRING_ELT(ref_base_vec, i, mkChar(ref_str));
                    SET_STRING_ELT(context_vec, i, mkChar(&contexts[i * 4]));
                }

                free(ref_bases);
                free(contexts);
            }
            cleanup_reference_v2();
        }
    }

    // Create coordinate vectors
    SEXP chr_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));
    SEXP pos_vec = PROTECT(allocVector(INTSXP, coord_set->n_sites));
    SEXP strand_vec = PROTECT(allocVector(STRSXP, coord_set->n_sites));

    for (size_t i = 0; i < coord_set->n_sites; i++) {
        SET_STRING_ELT(chr_vec, i, mkChar(coord_set->sites[i].chr));
        INTEGER(pos_vec)[i] = coord_set->sites[i].pos;
        char strand_str[2] = {coord_set->sites[i].strand, '\0'};
        SET_STRING_ELT(strand_vec, i, mkChar(strand_str));
    }

    // Return list with metadata and pointers
    int n_elements = has_ref_fasta && (ref_base_vec != R_NilValue) ? 9 : 7;
    SEXP result_list = PROTECT(allocVector(VECSXP, n_elements));
    SET_VECTOR_ELT(result_list, 0, coord_ptr);      // External pointer to coord_set
    SET_VECTOR_ELT(result_list, 1, chunks_ptr);     // External pointer to chunks
    SET_VECTOR_ELT(result_list, 2, chr_vec);
    SET_VECTOR_ELT(result_list, 3, pos_vec);
    SET_VECTOR_ELT(result_list, 4, strand_vec);
    SET_VECTOR_ELT(result_list, 5, ScalarInteger(chunks->n_chunks));
    SET_VECTOR_ELT(result_list, 6, ScalarInteger(n_files));

    SEXP names = PROTECT(allocVector(STRSXP, n_elements));
    SET_STRING_ELT(names, 0, mkChar("coord_ptr"));
    SET_STRING_ELT(names, 1, mkChar("chunks_ptr"));
    SET_STRING_ELT(names, 2, mkChar("chr"));
    SET_STRING_ELT(names, 3, mkChar("start"));
    SET_STRING_ELT(names, 4, mkChar("strand"));
    SET_STRING_ELT(names, 5, mkChar("n_chunks"));
    SET_STRING_ELT(names, 6, mkChar("n_samples"));

    if (n_elements == 9) {
        SET_VECTOR_ELT(result_list, 7, ref_base_vec);
        SET_VECTOR_ELT(result_list, 8, context_vec);
        SET_STRING_ELT(names, 7, mkChar("ref_base"));
        SET_STRING_ELT(names, 8, mkChar("context"));
    }

    setAttrib(result_list, R_NamesSymbol, names);

    int n_protect = (n_elements == 9) ? 8 : 6;
    UNPROTECT(n_protect);
    return result_list;
}

// ============================================================================
// STREAMING HDF5 SUPPORT: PROCESS SINGLE CHUNK
// ============================================================================

SEXP read_modkit_v2_process_chunk_c(SEXP coord_ptr, SEXP chunks_ptr, SEXP chunk_id,
                                     SEXP files, SEXP target_mod, SEXP min_coverage,
                                     SEXP quality_filter, SEXP combine_strands, SEXP verbose) {

    // Extract pointers
    coordinate_set_t *coord_set = (coordinate_set_t*)R_ExternalPtrAddr(coord_ptr);
    chunk_list_t *chunks = (chunk_list_t*)R_ExternalPtrAddr(chunks_ptr);

    if (!coord_set || !chunks) {
        error("Invalid pointers - chunk data may have been freed");
    }

    int chunk_idx = asInteger(chunk_id);
    if (chunk_idx < 0 || chunk_idx >= chunks->n_chunks) {
        error("Invalid chunk_id: %d (must be 0 to %d)", chunk_idx, chunks->n_chunks - 1);
    }

    // Parse arguments
    int n_files = length(files);
    char **input_files = (char**)R_alloc(n_files, sizeof(char*));
    for (int i = 0; i < n_files; i++) {
        input_files[i] = (char*)CHAR(STRING_ELT(files, i));
    }

    discovery_config_t config;
    config.target_mod = CHAR(STRING_ELT(target_mod, 0))[0];
    config.min_coverage = asInteger(min_coverage);
    config.quality_filter = asReal(quality_filter);
    config.combine_strands = asLogical(combine_strands);
    config.verbose = asLogical(verbose);

    // Process the requested chunk
    chunk_result_t *result = process_chunk(&chunks->chunks[chunk_idx], coord_set,
                                          input_files, n_files, &config);

    if (!result) {
        error("Failed to process chunk %d", chunk_idx);
    }

    // Convert chunk result to R matrices
    SEXP beta_matrix = PROTECT(allocMatrix(REALSXP, result->n_sites, n_files));
    SEXP cov_matrix = PROTECT(allocMatrix(INTSXP, result->n_sites, n_files));

    double *beta_ptr = REAL(beta_matrix);
    int *cov_ptr = INTEGER(cov_matrix);

    // Initialize with NA/0
    for (size_t i = 0; i < (size_t)result->n_sites * n_files; i++) {
        beta_ptr[i] = R_NaN;
        cov_ptr[i] = 0;
    }

    // Copy data (row-major to column-major)
    for (int32_t local_idx = 0; local_idx < result->n_sites; local_idx++) {
        for (int s = 0; s < n_files; s++) {
            size_t src_idx = (size_t)local_idx * n_files + s;
            size_t dst_idx = local_idx + (size_t)s * result->n_sites;

            if (get_bit(result->has_data_bits, src_idx)) {
                beta_ptr[dst_idx] = result->beta_values[src_idx];
                cov_ptr[dst_idx] = result->cov_values[src_idx];
            }
        }
    }

    free_chunk_result(result);

    // Return list with matrices and metadata
    SEXP result_list = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result_list, 0, beta_matrix);
    SET_VECTOR_ELT(result_list, 1, cov_matrix);
    SET_VECTOR_ELT(result_list, 2, ScalarInteger(chunks->chunks[chunk_idx].start_site_idx));

    SEXP names = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("cov"));
    SET_STRING_ELT(names, 2, mkChar("start_idx"));

    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(4);
    return result_list;
}

// ============================================================================
// STREAMING HDF5 SUPPORT: CLEANUP
// ============================================================================

SEXP read_modkit_v2_cleanup_c(SEXP coord_ptr, SEXP chunks_ptr) {
    coordinate_set_t *coord_set = (coordinate_set_t*)R_ExternalPtrAddr(coord_ptr);
    chunk_list_t *chunks = (chunk_list_t*)R_ExternalPtrAddr(chunks_ptr);

    if (coord_set) {
        free_coordinate_set(coord_set);
        R_ClearExternalPtr(coord_ptr);
    }

    if (chunks) {
        free_chunk_list(chunks);
        R_ClearExternalPtr(chunks_ptr);
    }

    return R_NilValue;
}