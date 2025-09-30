#include "modkit_merge_v2.h"
#include <stdlib.h>
#include <string.h>

// ============================================================================
// CHUNKING STRATEGY
// ============================================================================

chunk_list_t* create_chunks(coordinate_set_t *coord_set,
                           size_t sites_per_chunk,
                           int32_t max_genomic_span) {

    if (!coord_set || coord_set->n_sites == 0) {
        return NULL;
    }

    // Estimate number of chunks
    size_t estimated_chunks = (coord_set->n_sites + sites_per_chunk - 1) / sites_per_chunk;
    chunk_list_t *chunk_list = (chunk_list_t*)malloc(sizeof(chunk_list_t));
    if (!chunk_list) return NULL;

    chunk_list->chunks = (genomic_chunk_t*)malloc(estimated_chunks * 2 * sizeof(genomic_chunk_t));
    if (!chunk_list->chunks) {
        free(chunk_list);
        return NULL;
    }

    chunk_list->n_chunks = 0;
    chunk_list->sites_per_chunk = sites_per_chunk;

    int chunk_id = 0;

    // Create chunks within each chromosome
    for (int chr_idx = 0; chr_idx < coord_set->n_chromosomes; chr_idx++) {
        const char *chr = coord_set->chr_names[chr_idx];
        int32_t chr_start = coord_set->chr_start_idx[chr_idx];
        int32_t chr_end = coord_set->chr_end_idx[chr_idx];
        int32_t chr_n_sites = chr_end - chr_start;

        if (chr_n_sites == 0) continue;

        int32_t chunk_start_idx = chr_start;

        while (chunk_start_idx < chr_end) {
            genomic_chunk_t chunk;
            strncpy(chunk.chr, chr, MAX_CHR_LEN - 1);
            chunk.chr[MAX_CHR_LEN - 1] = '\0';

            chunk.start_site_idx = chunk_start_idx;
            chunk.start_pos = coord_set->sites[chunk_start_idx].pos;

            // Determine chunk end
            int32_t chunk_end_idx = chunk_start_idx + sites_per_chunk;
            if (chunk_end_idx > chr_end) {
                chunk_end_idx = chr_end;
            }

            // Check genomic span constraint
            int32_t genomic_span = coord_set->sites[chunk_end_idx - 1].pos - chunk.start_pos;
            if (genomic_span > max_genomic_span && (chunk_end_idx - chunk_start_idx) > 1000) {
                // Chunk is too large genomically - reduce size
                // Binary search for appropriate end point
                int32_t left = chunk_start_idx + 1000;  // Minimum 1000 sites
                int32_t right = chunk_end_idx;

                while (left < right) {
                    int32_t mid = left + (right - left) / 2;
                    int32_t span = coord_set->sites[mid - 1].pos - chunk.start_pos;

                    if (span <= max_genomic_span) {
                        left = mid + 1;
                    } else {
                        right = mid;
                    }
                }

                chunk_end_idx = left;
            }

            chunk.end_site_idx = chunk_end_idx;
            chunk.end_pos = coord_set->sites[chunk_end_idx - 1].pos + 1;  // Exclusive end
            chunk.chunk_id = chunk_id++;

            // Add chunk to list
            chunk_list->chunks[chunk_list->n_chunks++] = chunk;

            // Move to next chunk
            chunk_start_idx = chunk_end_idx;
        }
    }

    // Shrink allocation to actual size
    chunk_list->chunks = (genomic_chunk_t*)realloc(chunk_list->chunks,
                                                   chunk_list->n_chunks * sizeof(genomic_chunk_t));

    return chunk_list;
}

void free_chunk_list(chunk_list_t *chunks) {
    if (!chunks) return;
    if (chunks->chunks) free(chunks->chunks);
    free(chunks);
}

// ============================================================================
// CHUNK RESULT MANAGEMENT
// ============================================================================

chunk_result_t* alloc_chunk_result(int32_t n_sites, int32_t n_samples) {
    chunk_result_t *result = (chunk_result_t*)malloc(sizeof(chunk_result_t));
    if (!result) return NULL;

    result->n_sites = n_sites;
    result->n_samples = n_samples;

    size_t matrix_size = (size_t)n_sites * (size_t)n_samples;

    result->beta_values = (double*)malloc(matrix_size * sizeof(double));
    result->cov_values = (int32_t*)malloc(matrix_size * sizeof(int32_t));

    // Bitmap for has_data: 1 bit per cell
    size_t bitmap_size = (matrix_size + 7) / 8;
    result->has_data_bits = (uint8_t*)calloc(bitmap_size, sizeof(uint8_t));

    if (!result->beta_values || !result->cov_values || !result->has_data_bits) {
        free_chunk_result(result);
        return NULL;
    }

    // Initialize: beta = NA, cov = 0, has_data = false
    for (size_t i = 0; i < matrix_size; i++) {
        result->beta_values[i] = R_NaN;
        result->cov_values[i] = 0;
    }

    return result;
}

void free_chunk_result(chunk_result_t *result) {
    if (!result) return;
    if (result->beta_values) free(result->beta_values);
    if (result->cov_values) free(result->cov_values);
    if (result->has_data_bits) free(result->has_data_bits);
    free(result);
}

// ============================================================================
// GENOMIC INTERVAL CHUNKING (ModKit-style)
// ============================================================================
// TODO: This is currently a stub that calls the old create_chunks()
// Will be properly implemented to use genomic intervals from chrom_sizes
// ============================================================================

chunk_list_t* create_chunks_by_interval(coordinate_set_t *coord_set,
                                       chrom_sizes_t *chrom_sizes,
                                       int32_t interval_bp) {
    // For now, use old chunking with reasonable defaults
    // This allows compilation while we implement HashMap discovery
    size_t sites_per_chunk = 50000;
    int32_t max_genomic_span = interval_bp;

    if (coord_set->n_sites < 100000) {
        sites_per_chunk = coord_set->n_sites / 4;
        if (sites_per_chunk < 10000) sites_per_chunk = coord_set->n_sites;
    }

    return create_chunks(coord_set, sites_per_chunk, max_genomic_span);
}