#include "modkit_merge_v2.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// Simple hash table for aggregation within a chunk
// Key: (position, strand, mod_code, sample_idx) - automatically deduplicates sites
// Value: aggregated counts

#define HASH_TABLE_SIZE 1048576  // 1M entries

// Hash key: encodes site identity and sample index
typedef struct {
    int32_t pos;               // Genomic position
    char strand;               // Strand character
    char mod_code;             // Modification code
    int file_idx;              // Sample/file index
} hash_key_full_t;

typedef struct hash_entry {
    hash_key_full_t key;       // Full key for deduplication
    site_counts_t counts;      // Aggregated counts
    struct hash_entry *next;   // Chain for collisions
} hash_entry_t;

typedef struct {
    hash_entry_t **buckets;
    int size;
} hash_table_t;

// ============================================================================
// HASH TABLE FOR AGGREGATION
// ============================================================================

static hash_table_t* create_hash_table(int size) {
    hash_table_t *ht = (hash_table_t*)malloc(sizeof(hash_table_t));
    if (!ht) return NULL;

    ht->buckets = (hash_entry_t**)calloc(size, sizeof(hash_entry_t*));
    if (!ht->buckets) {
        free(ht);
        return NULL;
    }

    ht->size = size;
    return ht;
}

static void free_hash_table(hash_table_t *ht) {
    if (!ht) return;

    for (int i = 0; i < ht->size; i++) {
        hash_entry_t *entry = ht->buckets[i];
        while (entry) {
            hash_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
    }

    free(ht->buckets);
    free(ht);
}

// Hash function for full key
static inline uint32_t hash_function(hash_key_full_t *key, int table_size) {
    // Combine position, strand, mod_code, file_idx into hash
    uint32_t h = (uint32_t)key->pos;
    h = h * 31 + (uint32_t)key->strand;
    h = h * 31 + (uint32_t)key->mod_code;
    h = h * 31 + (uint32_t)key->file_idx;

    // Mix bits
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = (h >> 16) ^ h;
    return h % table_size;
}

// Compare two full keys for equality
static inline int keys_equal(hash_key_full_t *a, hash_key_full_t *b) {
    return a->pos == b->pos &&
           a->strand == b->strand &&
           a->mod_code == b->mod_code &&
           a->file_idx == b->file_idx;
}

static site_counts_t* get_or_create(hash_table_t *ht, hash_key_full_t *key) {
    uint32_t bucket = hash_function(key, ht->size);
    hash_entry_t *entry = ht->buckets[bucket];

    // Search in chain
    while (entry) {
        if (keys_equal(&entry->key, key)) {
            return &entry->counts;
        }
        entry = entry->next;
    }

    // Create new entry
    hash_entry_t *new_entry = (hash_entry_t*)calloc(1, sizeof(hash_entry_t));
    if (!new_entry) return NULL;

    new_entry->key = *key;  // Copy full key
    new_entry->next = ht->buckets[bucket];
    ht->buckets[bucket] = new_entry;

    return &new_entry->counts;
}

// ============================================================================
// CHUNK PROCESSING
// ============================================================================

chunk_result_t* process_chunk(genomic_chunk_t *chunk,
                              coordinate_set_t *coord_set,
                              char **input_files,
                              int n_files,
                              discovery_config_t *config) {

    int32_t chunk_n_sites = chunk->end_site_idx - chunk->start_site_idx;

    if (config->verbose && chunk->chunk_id % 10 == 0) {
        Rprintf("  Processing chunk %d: %s:%d-%d (%d sites)\n",
                chunk->chunk_id, chunk->chr, chunk->start_pos, chunk->end_pos,
                chunk_n_sites);
    }

    if (config->verbose) {
        Rprintf("    coord_set=%p, coord_set->sites=%p, coord_set->n_sites=%zu\n",
                (void*)coord_set, (void*)coord_set->sites, coord_set->n_sites);
        Rprintf("    chunk->start_site_idx=%d, chunk->end_site_idx=%d\n",
                chunk->start_site_idx, chunk->end_site_idx);
    }

    // Allocate result structure
    if (config->verbose) {
        Rprintf("    Allocating result for %d sites, %d files...\n", chunk_n_sites, n_files);
    }
    chunk_result_t *result = alloc_chunk_result(chunk_n_sites, n_files);
    if (!result) {
        Rprintf("Error: Failed to allocate chunk result\n");
        return NULL;
    }
    if (config->verbose) {
        Rprintf("    Result allocated: %p\n", (void*)result);
    }

    // Create hash table for aggregation
    if (config->verbose) {
        Rprintf("    Creating hash table...\n");
    }
    hash_table_t *agg_table = create_hash_table(HASH_TABLE_SIZE);
    if (!agg_table) {
        free_chunk_result(result);
        return NULL;
    }
    if (config->verbose) {
        Rprintf("    Hash table created: %p\n", (void*)agg_table);
    }

    // Process each file (sample)
    for (int file_idx = 0; file_idx < n_files; file_idx++) {
        const char *filename = input_files[file_idx];

        if (config->verbose) {
            Rprintf("    Processing file %d: %s\n", file_idx, filename);
        }

        // Open file and index
        htsFile *fp = hts_open(filename, "r");
        if (!fp) {
            Rprintf("Warning: Cannot open %s for chunk processing\n", filename);
            continue;
        }

        tbx_t *idx = tbx_index_load(filename);
        if (!idx) {
            Rprintf("Warning: Cannot load index for %s\n", filename);
            hts_close(fp);
            continue;
        }

        // Query the specific genomic region
        char region[256];
        snprintf(region, sizeof(region), "%s:%d-%d", chunk->chr, chunk->start_pos, chunk->end_pos);

        if (config->verbose) {
            Rprintf("      Querying region: %s\n", region);
        }

        hts_itr_t *iter = tbx_itr_querys(idx, region);
        if (!iter) {
            // No data in this region for this file
            if (config->verbose) {
                Rprintf("      No data in region\n");
            }
            tbx_destroy(idx);
            hts_close(fp);
            continue;
        }

        if (config->verbose) {
            Rprintf("      Iterator created, reading records...\n");
        }

        // Read records in this region
        kstring_t line = {0, 0, 0};
        while (tbx_itr_next(fp, idx, iter, &line) >= 0) {
            site_key_t site;
            site_counts_t counts;

            if (parse_bedmethyl_line(line.s, &site, &counts) != 0) {
                continue;
            }

            // Filter by modification code
            if (site.mod_code != config->target_mod) {
                continue;
            }

            // Apply strand normalization if needed
            if (config->combine_strands &&
                (config->target_mod == 'm' || config->target_mod == 'h' ||
                 config->target_mod == 'c' || config->target_mod == 'g')) {
                if (site.strand == '-') {
                    site.pos -= 1;
                }
                site.strand = '+';
            }

            // Find this site in the chunk's coordinate range
            int local_idx = binary_search_site(coord_set->sites,
                                              chunk->start_site_idx,
                                              chunk->end_site_idx,
                                              &site);

            if (local_idx < 0) {
                // Site not in discovered coordinates (shouldn't happen)
                continue;
            }

            // Create hash key based on site identity and sample
            // This automatically deduplicates if the same site appears multiple times
            // (e.g., at chunk boundaries or from discovery duplicates)
            hash_key_full_t hash_key;
            hash_key.pos = site.pos;
            hash_key.strand = site.strand;
            hash_key.mod_code = site.mod_code;
            hash_key.file_idx = file_idx;

            // Get or create entry in hash table
            // If this exact (pos, strand, mod_code, file_idx) was seen before,
            // we get the existing entry and aggregate counts
            site_counts_t *agg_counts = get_or_create(agg_table, &hash_key);
            if (!agg_counts) {
                Rprintf("Error: Hash table full\n");
                break;
            }

            // Aggregate counts (for strand merging)
            agg_counts->n_mod += counts.n_mod;
            agg_counts->n_canonical += counts.n_canonical;
            agg_counts->valid_cov += counts.valid_cov;
        }

        if (line.s) free(line.s);
        tbx_itr_destroy(iter);
        tbx_destroy(idx);
        hts_close(fp);
    }

    // Finalize: compute beta values from aggregated counts
    // Map hash table entries back to matrix positions using coordinate lookup
    if (config->verbose) {
        Rprintf("    Finalizing: mapping hash entries to matrix...\n");
        Rprintf("    Hash table size=%d, start_idx=%d, end_idx=%d\n",
                agg_table->size, chunk->start_site_idx, chunk->end_site_idx);
    }

    int entries_processed = 0;
    int max_entries = chunk_n_sites * n_files * 2;  // Safety limit: 2x expected max
    int loop_broken = 0;
    for (int i = 0; i < agg_table->size && !loop_broken; i++) {
        if (!agg_table->buckets) {
            Rprintf("ERROR: buckets is NULL at i=%d\n", i);
            break;
        }
        hash_entry_t *entry = agg_table->buckets[i];
        while (entry) {
            // Sanity check: entry pointer should look valid
            uintptr_t entry_addr = (uintptr_t)entry;
            if (entry_addr < 0x1000 || (entry_addr & 0x7) != 0) {
                Rprintf("ERROR: Corrupted entry pointer %p in bucket %d (after %d entries), stopping\n",
                        entry, i, entries_processed);
                loop_broken = 1;
                break;
            }

            entries_processed++;
            if (entries_processed > max_entries) {
                Rprintf("ERROR: Too many entries (%d > %d), likely infinite loop!\n",
                        entries_processed, max_entries);
                loop_broken = 1;
                break;
            }
            if (config->verbose && entries_processed % 10000 == 0) {
                Rprintf("      Processed %d entries... (bucket %d)\n", entries_processed, i);
            }

            // Extract site identity from hash key
            site_key_t lookup_site;
            lookup_site.pos = entry->key.pos;
            lookup_site.strand = entry->key.strand;
            lookup_site.mod_code = entry->key.mod_code;
            // Note: chr is not in hash key, but we know it's chunk->chr
            strncpy(lookup_site.chr, chunk->chr, MAX_CHR_LEN - 1);
            lookup_site.chr[MAX_CHR_LEN - 1] = '\0';

            // Find this site in the coordinate set
            int local_idx = binary_search_site(coord_set->sites,
                                              chunk->start_site_idx,
                                              chunk->end_site_idx,
                                              &lookup_site);

            if (local_idx >= 0) {
                // Convert to chunk-local index
                int32_t chunk_local_idx = local_idx - chunk->start_site_idx;
                int file_idx = entry->key.file_idx;

                // Bounds check
                if (chunk_local_idx < 0 || chunk_local_idx >= chunk_n_sites) {
                    Rprintf("ERROR: chunk_local_idx=%d out of range [0, %d)\n",
                            chunk_local_idx, chunk_n_sites);
                    entry = entry->next;
                    continue;
                }
                if (file_idx < 0 || file_idx >= n_files) {
                    Rprintf("ERROR: file_idx=%d out of range [0, %d)\n",
                            file_idx, n_files);
                    entry = entry->next;
                    continue;
                }

                // Compute matrix index (row-major: local_idx * n_samples + file_idx)
                size_t matrix_idx = (size_t)chunk_local_idx * (size_t)n_files + (size_t)file_idx;

                // Compute beta
                result->beta_values[matrix_idx] = compute_beta(&entry->counts);
                result->cov_values[matrix_idx] = entry->counts.valid_cov;

                // Mark as having data
                set_bit(result->has_data_bits, matrix_idx);
            }

            // Safely advance to next entry with corruption detection
            if (!entry) {
                // Already NULL, loop will end
                break;
            } else if (!entry->next) {
                // No next entry, end of chain
                entry = NULL;
            } else {
                // Check if next pointer looks valid before dereferencing
                uintptr_t next_addr = (uintptr_t)entry->next;
                if (next_addr < 0x1000 || (next_addr & 0x7) != 0) {
                    // Looks like a corrupted pointer (NULL-ish or misaligned)
                    Rprintf("WARNING: Corrupted next pointer %p in bucket %d (entry %d), stopping chain\n",
                            entry->next, i, entries_processed);
                    entry = NULL;
                } else {
                    entry = entry->next;
                }
            }
        }
    }

    if (config->verbose) {
        Rprintf("    Finalization complete, processed %d entries total\n", entries_processed);
        Rprintf("    Freeing hash table...\n");
    }

    free_hash_table(agg_table);

    if (config->verbose) {
        Rprintf("    Hash table freed, returning result\n");
    }

    // Periodic interrupt check
    if (chunk->chunk_id % 10 == 0) {
        R_CheckUserInterrupt();
    }

    return result;
}