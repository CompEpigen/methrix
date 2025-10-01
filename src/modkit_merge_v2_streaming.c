#include "modkit_merge_v2.h"
#include <string.h>
#include <stdlib.h>

// ============================================================================
// STREAMING ARCHITECTURE: INTERVAL-BASED PROCESSING
// ============================================================================
// Instead of discovering ALL sites then chunking, we process genome interval
// by interval. For each interval:
//   1. Discover sites in that interval (from all files)
//   2. Process the interval
//   3. Return results to R
//   4. Free memory
// This keeps memory bounded regardless of dataset size.
// ============================================================================

// Structure to hold one interval's result
typedef struct {
    char chr[MAX_CHR_LEN];
    int32_t start_pos;         // Genomic start
    int32_t end_pos;           // Genomic end
    int32_t n_sites;           // Number of sites discovered

    // Coordinate arrays
    int32_t *positions;        // Site positions
    char *strands;             // Site strands

    // Data matrices (row-major: site × sample)
    double *beta_values;
    int32_t *cov_values;
} interval_result_t;

// Free interval result
static void free_interval_result(interval_result_t *result) {
    if (!result) return;
    if (result->positions) free(result->positions);
    if (result->strands) free(result->strands);
    if (result->beta_values) free(result->beta_values);
    if (result->cov_values) free(result->cov_values);
    free(result);
}

// ============================================================================
// PROCESS ONE GENOMIC INTERVAL
// ============================================================================

interval_result_t* process_genomic_interval(
    const char *chr,
    int32_t interval_start,
    int32_t interval_end,
    char **input_files,
    int n_files,
    discovery_config_t *config
) {

    // Step 1: Discover coordinates in this interval using existing function
    // We'll use the coordinate discovery machinery but with tabix region query

    coordinate_set_t *coord_set = (coordinate_set_t*)malloc(sizeof(coordinate_set_t));
    if (!coord_set) return NULL;

    coord_set->capacity = INITIAL_SITES_CAPACITY;
    coord_set->sites = (site_key_t*)malloc(coord_set->capacity * sizeof(site_key_t));
    if (!coord_set->sites) {
        free(coord_set);
        return NULL;
    }
    coord_set->n_sites = 0;

    // Initialize chromosome index fields (required for free_coordinate_set)
    coord_set->chr_start_idx = NULL;
    coord_set->chr_end_idx = NULL;
    coord_set->chr_names = NULL;
    coord_set->n_chromosomes = 0;

    // Discover sites in this interval from all files
    site_hash_table_t *ht = create_site_hash_table(DISCOVERY_HASH_SIZE);
    if (!ht) {
        free(coord_set->sites);
        free(coord_set);
        return NULL;
    }

    size_t duplicates_skipped = 0;
    size_t total_records = 0;

    for (int f = 0; f < n_files; f++) {
        // Open file
        htsFile *fp = hts_open(input_files[f], "r");
        if (!fp) continue;

        tbx_t *idx = tbx_index_load(input_files[f]);
        if (!idx) {
            hts_close(fp);
            continue;
        }

        // Query this specific interval using region string
        char region[MAX_CHR_LEN + 64];
        snprintf(region, sizeof(region), "%s:%d-%d", chr, interval_start, interval_end);

        hts_itr_t *iter = tbx_itr_querys(idx, region);
        if (!iter) {
            tbx_destroy(idx);
            hts_close(fp);
            continue;
        }

        kstring_t line = {0, 0, NULL};

        while (tbx_itr_next(fp, idx, iter, &line) >= 0) {
            total_records++;

            // Parse bedMethyl line
            site_key_t site;
            site_counts_t counts;

            if (parse_bedmethyl_line(line.s, &site, &counts) != 0) {
                continue;
            }

            // Check if site is for target modification
            if (site.mod_code != config->target_mod) {
                continue;
            }

            // Check if position is within our interval
            if (site.pos < interval_start || site.pos >= interval_end) {
                continue;
            }

            // Check quality filters
            if (!passes_filters(&counts, config)) {
                continue;
            }

            // Combine strands if requested (for CpG modifications)
            if (config->combine_strands &&
                (config->target_mod == 'm' || config->target_mod == 'h' ||
                 config->target_mod == 'c' || config->target_mod == 'g')) {

                // For CpG modifications: normalize to + strand
                // Minus strand G at pos N+1 -> Plus strand C at pos N
                if (site.strand == '-') {
                    site.pos -= 1;
                }
                site.strand = '+';
            }

            // Try to insert (automatically deduplicates)
            int inserted = insert_unique_site(ht, &site);
            if (!inserted) {
                duplicates_skipped++;
            }
        }

        free(line.s);
        hts_itr_destroy(iter);
        tbx_destroy(idx);
        hts_close(fp);
    }


    // Extract sites from hash table into coordinate set
    coord_set->n_sites = 0;
    for (int i = 0; i < ht->size; i++) {
        hash_site_entry_t *entry = ht->buckets[i];
        while (entry) {
            if (coord_set->n_sites >= coord_set->capacity) {
                coord_set->capacity *= 2;
                site_key_t *new_sites = (site_key_t*)realloc(
                    coord_set->sites,
                    coord_set->capacity * sizeof(site_key_t)
                );
                if (!new_sites) {
                    free_site_hash_table(ht);
                    free_coordinate_set(coord_set);
                    return NULL;
                }
                coord_set->sites = new_sites;
            }

            // Copy site and ensure chr is set correctly
            coord_set->sites[coord_set->n_sites] = entry->site;
            // Ensure chr field is set (it should already be from parsing, but be safe)
            strncpy(coord_set->sites[coord_set->n_sites].chr, chr, MAX_CHR_LEN - 1);
            coord_set->sites[coord_set->n_sites].chr[MAX_CHR_LEN - 1] = '\0';
            coord_set->n_sites++;

            entry = entry->next;
        }
    }

    free_site_hash_table(ht);

    // If no sites found, return empty result
    if (coord_set->n_sites == 0) {
        free_coordinate_set(coord_set);

        interval_result_t *result = (interval_result_t*)calloc(1, sizeof(interval_result_t));
        strncpy(result->chr, chr, MAX_CHR_LEN - 1);
        result->start_pos = interval_start;
        result->end_pos = interval_end;
        result->n_sites = 0;
        return result;
    }

    // Sort sites by position, strand, mod_code
    qsort(coord_set->sites, coord_set->n_sites, sizeof(site_key_t),
          (int (*)(const void*, const void*))compare_sites);

    // Step 2: Process the interval (aggregate data from all files)
    // Create a chunk structure for this interval
    genomic_chunk_t chunk;
    strncpy(chunk.chr, chr, MAX_CHR_LEN - 1);
    chunk.start_site_idx = 0;
    chunk.end_site_idx = coord_set->n_sites;
    chunk.start_pos = interval_start;
    chunk.end_pos = interval_end;
    chunk.chunk_id = 0;

    // Process the chunk
    chunk_result_t *chunk_result = process_chunk(&chunk, coord_set, input_files, n_files, config);

    if (!chunk_result) {
        free_coordinate_set(coord_set);
        return NULL;
    }

    // Step 3: Package results for R
    interval_result_t *result = (interval_result_t*)malloc(sizeof(interval_result_t));
    if (!result) {
        free_chunk_result(chunk_result);
        free_coordinate_set(coord_set);
        return NULL;
    }

    strncpy(result->chr, chr, MAX_CHR_LEN - 1);
    result->start_pos = interval_start;
    result->end_pos = interval_end;
    result->n_sites = coord_set->n_sites;

    // Copy coordinates
    result->positions = (int32_t*)malloc(coord_set->n_sites * sizeof(int32_t));
    result->strands = (char*)malloc(coord_set->n_sites * sizeof(char));

    for (size_t i = 0; i < coord_set->n_sites; i++) {
        result->positions[i] = coord_set->sites[i].pos;
        result->strands[i] = coord_set->sites[i].strand;
    }

    // Copy data matrices
    size_t matrix_size = (size_t)coord_set->n_sites * n_files;
    result->beta_values = (double*)malloc(matrix_size * sizeof(double));
    result->cov_values = (int32_t*)malloc(matrix_size * sizeof(int32_t));

    memcpy(result->beta_values, chunk_result->beta_values, matrix_size * sizeof(double));
    memcpy(result->cov_values, chunk_result->cov_values, matrix_size * sizeof(int32_t));

    // Cleanup
    free_chunk_result(chunk_result);
    free_coordinate_set(coord_set);


    return result;
}

// ============================================================================
// R INTERFACE: PROCESS ONE INTERVAL
// ============================================================================

SEXP read_modkit_v2_process_interval_c(
    SEXP chr, SEXP interval_start, SEXP interval_end,
    SEXP files, SEXP target_mod, SEXP min_coverage,
    SEXP quality_filter, SEXP combine_strands, SEXP verbose
) {

    // Parse arguments
    const char *chr_str = CHAR(STRING_ELT(chr, 0));
    int32_t start_pos = asInteger(interval_start);
    int32_t end_pos = asInteger(interval_end);

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

    // Process the interval
    interval_result_t *result = process_genomic_interval(
        chr_str, start_pos, end_pos, input_files, n_files, &config
    );

    if (!result) {
        error("Failed to process interval %s:%d-%d", chr_str, start_pos, end_pos);
    }

    // Convert to R objects
    if (result->n_sites == 0) {
        // Return empty result
        SEXP result_list = PROTECT(allocVector(VECSXP, 6));
        SET_VECTOR_ELT(result_list, 0, ScalarString(mkChar(result->chr)));
        SET_VECTOR_ELT(result_list, 1, ScalarInteger(result->start_pos));
        SET_VECTOR_ELT(result_list, 2, ScalarInteger(result->end_pos));
        SET_VECTOR_ELT(result_list, 3, ScalarInteger(0));  // n_sites
        SET_VECTOR_ELT(result_list, 4, R_NilValue);  // No beta matrix
        SET_VECTOR_ELT(result_list, 5, R_NilValue);  // No cov matrix

        SEXP names = PROTECT(allocVector(STRSXP, 6));
        SET_STRING_ELT(names, 0, mkChar("chr"));
        SET_STRING_ELT(names, 1, mkChar("start_pos"));
        SET_STRING_ELT(names, 2, mkChar("end_pos"));
        SET_STRING_ELT(names, 3, mkChar("n_sites"));
        SET_STRING_ELT(names, 4, mkChar("beta"));
        SET_STRING_ELT(names, 5, mkChar("cov"));
        setAttrib(result_list, R_NamesSymbol, names);

        free_interval_result(result);
        UNPROTECT(2);
        return result_list;
    }

    // Create position and strand vectors
    SEXP pos_vec = PROTECT(allocVector(INTSXP, result->n_sites));
    SEXP strand_vec = PROTECT(allocVector(STRSXP, result->n_sites));

    for (int32_t i = 0; i < result->n_sites; i++) {
        INTEGER(pos_vec)[i] = result->positions[i];
        char strand_str[2] = {result->strands[i], '\0'};
        SET_STRING_ELT(strand_vec, i, mkChar(strand_str));
    }

    // Create matrices (convert from row-major to column-major for R)
    SEXP beta_matrix = PROTECT(allocMatrix(REALSXP, result->n_sites, n_files));
    SEXP cov_matrix = PROTECT(allocMatrix(INTSXP, result->n_sites, n_files));

    double *beta_ptr = REAL(beta_matrix);
    int *cov_ptr = INTEGER(cov_matrix);

    for (int32_t i = 0; i < result->n_sites; i++) {
        for (int j = 0; j < n_files; j++) {
            // Source: row-major (i * n_files + j)
            size_t src_idx = (size_t)i * n_files + j;
            // Dest: column-major (i + j * n_sites)
            size_t dst_idx = i + (size_t)j * result->n_sites;

            beta_ptr[dst_idx] = result->beta_values[src_idx];
            cov_ptr[dst_idx] = result->cov_values[src_idx];
        }
    }

    // Create result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 8));
    SET_VECTOR_ELT(result_list, 0, ScalarString(mkChar(result->chr)));
    SET_VECTOR_ELT(result_list, 1, ScalarInteger(result->start_pos));
    SET_VECTOR_ELT(result_list, 2, ScalarInteger(result->end_pos));
    SET_VECTOR_ELT(result_list, 3, ScalarInteger(result->n_sites));
    SET_VECTOR_ELT(result_list, 4, pos_vec);
    SET_VECTOR_ELT(result_list, 5, strand_vec);
    SET_VECTOR_ELT(result_list, 6, beta_matrix);
    SET_VECTOR_ELT(result_list, 7, cov_matrix);

    SEXP names = PROTECT(allocVector(STRSXP, 8));
    SET_STRING_ELT(names, 0, mkChar("chr"));
    SET_STRING_ELT(names, 1, mkChar("start_pos"));
    SET_STRING_ELT(names, 2, mkChar("end_pos"));
    SET_STRING_ELT(names, 3, mkChar("n_sites"));
    SET_STRING_ELT(names, 4, mkChar("positions"));
    SET_STRING_ELT(names, 5, mkChar("strand"));
    SET_STRING_ELT(names, 6, mkChar("beta"));
    SET_STRING_ELT(names, 7, mkChar("cov"));
    setAttrib(result_list, R_NamesSymbol, names);

    free_interval_result(result);
    UNPROTECT(6);  // pos_vec, strand_vec, beta_matrix, cov_matrix, result_list, names
    return result_list;
}

// ============================================================================
// R INTERFACE: CREATE INTERVAL LIST FROM CHROMOSOME SIZES
// ============================================================================

SEXP read_modkit_v2_create_intervals_c(
    SEXP chrom_sizes_file, SEXP files, SEXP interval_size, SEXP verbose
) {

    const char *chrom_sizes_path = CHAR(STRING_ELT(chrom_sizes_file, 0));
    int32_t interval_bp = asInteger(interval_size);
    int verb = asLogical(verbose);

    int n_files = length(files);
    char **input_files = (char**)R_alloc(n_files, sizeof(char*));
    for (int i = 0; i < n_files; i++) {
        input_files[i] = (char*)CHAR(STRING_ELT(files, i));
    }

    if (verb) {
        Rprintf("Creating genomic intervals for streaming processing...\n");
        Rprintf("  Interval size: %d bp\n", interval_bp);
    }

    // Read chromosome sizes
    chrom_sizes_t *all_sizes = read_chrom_sizes(chrom_sizes_path);
    if (!all_sizes) {
        error("Failed to read chromosome sizes file: %s", chrom_sizes_path);
    }

    // Filter to chromosomes present in tabix files
    chrom_sizes_t *filtered_sizes = filter_chrom_sizes_by_tabix(
        all_sizes, input_files, n_files, verb
    );
    free_chrom_sizes(all_sizes);

    if (!filtered_sizes) {
        error("No chromosomes match between sizes file and input files");
    }

    if (verb) {
        Rprintf("  Chromosomes to process: %d\n", filtered_sizes->n_chroms);
    }

    // Count total intervals needed
    int total_intervals = 0;
    for (int i = 0; i < filtered_sizes->n_chroms; i++) {
        int32_t chr_len = filtered_sizes->chroms[i].length;
        int n_intervals = (chr_len + interval_bp - 1) / interval_bp;  // ceiling division
        total_intervals += n_intervals;
    }

    if (verb) {
        Rprintf("  Total intervals: %d\n", total_intervals);
    }

    // Allocate R vectors for intervals
    SEXP chr_vec = PROTECT(allocVector(STRSXP, total_intervals));
    SEXP start_vec = PROTECT(allocVector(INTSXP, total_intervals));
    SEXP end_vec = PROTECT(allocVector(INTSXP, total_intervals));
    SEXP chr_len_vec = PROTECT(allocVector(INTSXP, total_intervals));

    // Generate intervals
    int interval_idx = 0;
    for (int chr_idx = 0; chr_idx < filtered_sizes->n_chroms; chr_idx++) {
        const char *chr_name = filtered_sizes->chroms[chr_idx].name;
        int32_t chr_len = filtered_sizes->chroms[chr_idx].length;

        if (verb) {
            Rprintf("  %s: length=%d bp, ", chr_name, chr_len);
        }

        int chr_intervals = 0;
        for (int32_t start = 0; start < chr_len; start += interval_bp) {
            int32_t end = start + interval_bp;
            if (end > chr_len) {
                end = chr_len;
            }

            SET_STRING_ELT(chr_vec, interval_idx, mkChar(chr_name));
            INTEGER(start_vec)[interval_idx] = start;
            INTEGER(end_vec)[interval_idx] = end;
            INTEGER(chr_len_vec)[interval_idx] = chr_len;

            interval_idx++;
            chr_intervals++;
        }

        if (verb) {
            Rprintf("intervals=%d\n", chr_intervals);
        }
    }

    free_chrom_sizes(filtered_sizes);

    // Create data.frame
    SEXP result_df = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result_df, 0, chr_vec);
    SET_VECTOR_ELT(result_df, 1, start_vec);
    SET_VECTOR_ELT(result_df, 2, end_vec);
    SET_VECTOR_ELT(result_df, 3, chr_len_vec);

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("chr"));
    SET_STRING_ELT(names, 1, mkChar("start"));
    SET_STRING_ELT(names, 2, mkChar("end"));
    SET_STRING_ELT(names, 3, mkChar("chr_length"));
    setAttrib(result_df, R_NamesSymbol, names);

    // Set row names
    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    INTEGER(rownames)[0] = NA_INTEGER;
    INTEGER(rownames)[1] = -total_intervals;
    setAttrib(result_df, R_RowNamesSymbol, rownames);

    // Set class to data.frame
    SEXP df_class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(df_class, 0, mkChar("data.frame"));
    setAttrib(result_df, R_ClassSymbol, df_class);

    if (verb) {
        Rprintf("Interval creation complete\n");
    }

    UNPROTECT(8);
    return result_df;
}