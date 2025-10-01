#include "modkit_merge_v2.h"
#include <string.h>
#include <stdlib.h>

// strdup replacement for portability
static char* my_strdup(const char *s) {
    size_t len = strlen(s) + 1;
    char *new_str = (char*)malloc(len);
    if (new_str) {
        memcpy(new_str, s, len);
    }
    return new_str;
}

// ============================================================================
// HASH TABLE FOR COORDINATE DISCOVERY (ModKit-style)
// ============================================================================
// Key: (position, strand) - NOT including mod_code
// This guarantees automatic deduplication of sites
// ============================================================================
// Note: hash_site_entry_t, site_hash_table_t, and site_hash_key_t
// are all defined in modkit_merge_v2.h

site_hash_table_t* create_site_hash_table(int size) {
    site_hash_table_t *ht = (site_hash_table_t*)malloc(sizeof(site_hash_table_t));
    if (!ht) return NULL;

    ht->buckets = (hash_site_entry_t**)calloc(size, sizeof(hash_site_entry_t*));
    if (!ht->buckets) {
        free(ht);
        return NULL;
    }

    ht->size = size;
    ht->n_entries = 0;
    return ht;
}

void free_site_hash_table(site_hash_table_t *ht) {
    if (!ht) return;

    for (int i = 0; i < ht->size; i++) {
        hash_site_entry_t *entry = ht->buckets[i];
        while (entry) {
            hash_site_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
    }

    free(ht->buckets);
    free(ht);
}

// Hash function for (pos, strand)
static inline uint32_t hash_site_key(site_hash_key_t *key, int table_size) {
    uint32_t h = (uint32_t)key->pos;
    h = h * 31 + (uint32_t)key->strand;

    // Mix bits
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = (h >> 16) ^ h;
    return h % table_size;
}

// Insert or skip if exists (returns 1 if new, 0 if duplicate)
int insert_unique_site(site_hash_table_t *ht, const site_key_t *site) {
    site_hash_key_t key;
    key.pos = site->pos;
    key.strand = site->strand;

    uint32_t bucket = hash_site_key(&key, ht->size);
    hash_site_entry_t *entry = ht->buckets[bucket];

    // Search for existing entry with same (pos, strand)
    while (entry) {
        if (entry->key.pos == key.pos && entry->key.strand == key.strand) {
            // Duplicate found - ignore
            return 0;
        }
        entry = entry->next;
    }

    // New site - insert
    hash_site_entry_t *new_entry = (hash_site_entry_t*)malloc(sizeof(hash_site_entry_t));
    if (!new_entry) return -1;

    new_entry->key = key;
    new_entry->site = *site;
    new_entry->next = ht->buckets[bucket];
    ht->buckets[bucket] = new_entry;
    ht->n_entries++;

    return 1;
}

// ============================================================================
// SITE COMPARISON: Since tabix files are sorted, we leverage this!
// ============================================================================

int compare_sites(const site_key_t *a, const site_key_t *b) {
    // Compare: chromosome -> position -> strand -> mod_code
    int chr_cmp = strcmp(a->chr, b->chr);
    if (chr_cmp != 0) return chr_cmp;

    if (a->pos != b->pos) {
        return (a->pos < b->pos) ? -1 : 1;
    }

    if (a->strand != b->strand) {
        return (a->strand < b->strand) ? -1 : 1;
    }

    if (a->mod_code != b->mod_code) {
        return (a->mod_code < b->mod_code) ? -1 : 1;
    }

    return 0;
}

// ============================================================================
// COORDINATE SET MANAGEMENT
// ============================================================================

static coordinate_set_t* init_coordinate_set(size_t initial_capacity) {
    coordinate_set_t *cs = (coordinate_set_t*)malloc(sizeof(coordinate_set_t));
    if (!cs) return NULL;

    cs->sites = (site_key_t*)malloc(initial_capacity * sizeof(site_key_t));
    if (!cs->sites) {
        free(cs);
        return NULL;
    }

    cs->n_sites = 0;
    cs->capacity = initial_capacity;
    cs->n_chromosomes = 0;
    cs->chr_start_idx = NULL;
    cs->chr_end_idx = NULL;
    cs->chr_names = NULL;

    return cs;
}

static int grow_coordinate_set(coordinate_set_t *cs) {
    size_t new_capacity = cs->capacity * 2;

    // Check for overflow
    if (new_capacity < cs->capacity) {
        Rprintf("Error: Capacity overflow in coordinate set\n");
        return -1;
    }

    site_key_t *new_sites = (site_key_t*)realloc(cs->sites,
                                                  new_capacity * sizeof(site_key_t));
    if (!new_sites) {
        Rprintf("Error: Failed to reallocate coordinate set (tried %zu -> %zu)\n",
                cs->capacity, new_capacity);
        return -1;
    }

    cs->sites = new_sites;
    cs->capacity = new_capacity;
    return 0;
}

static int add_site_to_set(coordinate_set_t *cs, const site_key_t *site) {
    if (cs->n_sites >= cs->capacity) {
        if (grow_coordinate_set(cs) != 0) {
            return -1;
        }
    }

    cs->sites[cs->n_sites] = *site;
    cs->n_sites++;
    return 0;
}

void free_coordinate_set(coordinate_set_t *coord_set) {
    if (!coord_set) return;

    if (coord_set->sites) free(coord_set->sites);
    if (coord_set->chr_start_idx) free(coord_set->chr_start_idx);
    if (coord_set->chr_end_idx) free(coord_set->chr_end_idx);

    if (coord_set->chr_names) {
        for (int i = 0; i < coord_set->n_chromosomes; i++) {
            if (coord_set->chr_names[i]) free(coord_set->chr_names[i]);
        }
        free(coord_set->chr_names);
    }

    free(coord_set);
}

// Build chromosome index after all sites discovered
static int build_chromosome_index(coordinate_set_t *cs) {
    if (cs->n_sites == 0) return 0;

    // Count unique chromosomes
    int n_chr = 1;
    for (size_t i = 1; i < cs->n_sites; i++) {
        if (strcmp(cs->sites[i].chr, cs->sites[i-1].chr) != 0) {
            n_chr++;
        }
    }

    // Allocate chromosome arrays
    cs->chr_names = (char**)malloc(n_chr * sizeof(char*));
    cs->chr_start_idx = (int32_t*)malloc(n_chr * sizeof(int32_t));
    cs->chr_end_idx = (int32_t*)malloc(n_chr * sizeof(int32_t));

    if (!cs->chr_names || !cs->chr_start_idx || !cs->chr_end_idx) {
        return -1;
    }

    // Fill chromosome index
    int chr_idx = 0;
    cs->chr_names[0] = my_strdup(cs->sites[0].chr);
    cs->chr_start_idx[0] = 0;

    for (size_t i = 1; i < cs->n_sites; i++) {
        if (strcmp(cs->sites[i].chr, cs->sites[i-1].chr) != 0) {
            // End of previous chromosome
            cs->chr_end_idx[chr_idx] = i;

            // Start of new chromosome
            chr_idx++;
            cs->chr_names[chr_idx] = my_strdup(cs->sites[i].chr);
            cs->chr_start_idx[chr_idx] = i;
        }
    }
    cs->chr_end_idx[chr_idx] = cs->n_sites;
    cs->n_chromosomes = n_chr;

    return 0;
}

// ============================================================================
// DISCOVERY STREAM MANAGEMENT
// ============================================================================

static int init_discovery_stream(discovery_stream_t *stream,
                                 const char *filename,
                                 int file_idx) {
    memset(stream, 0, sizeof(discovery_stream_t));
    stream->file_idx = file_idx;
    stream->filename = my_strdup(filename);
    stream->has_site = 0;

    // Open file and index
    stream->fp = hts_open(filename, "r");
    if (!stream->fp) {
        Rprintf("Error: Cannot open file %s\n", filename);
        return -1;
    }

    stream->idx = tbx_index_load(filename);
    if (!stream->idx) {
        Rprintf("Error: Cannot load tabix index for %s\n", filename);
        hts_close(stream->fp);
        return -1;
    }

    // No iterator yet - will be created per-chromosome
    stream->iter = NULL;
    stream->line.s = NULL;
    stream->line.l = 0;
    stream->line.m = 0;

    return 0;
}

static void cleanup_discovery_stream(discovery_stream_t *stream) {
    if (stream->line.s) free(stream->line.s);
    if (stream->iter) tbx_itr_destroy(stream->iter);
    if (stream->idx) tbx_destroy(stream->idx);
    if (stream->fp) hts_close(stream->fp);
    if (stream->filename) free(stream->filename);
}

// Query next chromosome in stream
static int query_chromosome(discovery_stream_t *stream, const char *chr) {
    // Destroy previous iterator if exists
    if (stream->iter) {
        tbx_itr_destroy(stream->iter);
        stream->iter = NULL;
    }

    // Create iterator for this chromosome
    stream->iter = tbx_itr_querys(stream->idx, chr);
    if (!stream->iter) {
        // No data for this chromosome in this file (OK)
        stream->has_site = 0;
        return 0;
    }

    stream->has_site = 0;
    return 0;
}

// Advance stream to next valid site
static int advance_stream(discovery_stream_t *stream,
                         discovery_config_t *config) {
    if (!stream->iter) {
        stream->has_site = 0;
        return 0;
    }

    while (1) {
        stream->line.l = 0;
        int ret = tbx_itr_next(stream->fp, stream->idx, stream->iter, &stream->line);

        if (ret < 0) {
            // End of iterator
            stream->has_site = 0;
            return 0;
        }

        // Parse line
        site_counts_t counts;
        if (parse_bedmethyl_line(stream->line.s, &stream->current_site, &counts) != 0) {
            continue;  // Skip invalid lines
        }

        // Filter by modification code
        if (stream->current_site.mod_code != config->target_mod) {
            continue;
        }

        // Filter by quality/coverage
        if (!passes_filters(&counts, config)) {
            continue;
        }

        // Apply strand normalization if combining strands
        if (config->combine_strands &&
            (config->target_mod == 'm' || config->target_mod == 'h' ||
             config->target_mod == 'c' || config->target_mod == 'g')) {

            // For CpG modifications: normalize to + strand
            // Minus strand G at pos N+1 -> Plus strand C at pos N
            if (stream->current_site.strand == '-') {
                stream->current_site.pos -= 1;
            }
            stream->current_site.strand = '+';
        }

        // Valid site found
        stream->has_site = 1;
        return 0;
    }
}

// ============================================================================
// HASHMAP-BASED COORDINATE DISCOVERY (ModKit-style)
// ============================================================================

// Discover coordinates for a single chromosome using HashMap (ModKit-style)
static int discover_chromosome_sites(coordinate_set_t *cs,
                                     discovery_stream_t *streams,
                                     int n_streams,
                                     const char *chr,
                                     discovery_config_t *config) {

    if (config->verbose) {
        Rprintf("  Discovering sites on %s...\n", chr);
    }

    // Create hash table for this chromosome
    site_hash_table_t *ht = create_site_hash_table(DISCOVERY_HASH_SIZE);
    if (!ht) {
        Rprintf("Error: Failed to create hash table\n");
        return -1;
    }

    size_t total_records = 0;
    size_t duplicates_skipped = 0;

    // Process all files for this chromosome
    for (int file_idx = 0; file_idx < n_streams; file_idx++) {
        // Query this file for the chromosome
        query_chromosome(&streams[file_idx], chr);

        // Read all sites from this file
        while (1) {
            advance_stream(&streams[file_idx], config);

            if (!streams[file_idx].has_site) {
                // End of data for this file
                break;
            }

            total_records++;

            // Insert into hash table (automatically deduplicates by pos+strand)
            int result = insert_unique_site(ht, &streams[file_idx].current_site);
            if (result < 0) {
                Rprintf("Error: Hash table insertion failed\n");
                free_site_hash_table(ht);
                return -1;
            } else if (result == 0) {
                // Duplicate
                duplicates_skipped++;
            }

            // Periodic interrupt check
            if (total_records % 50000 == 0) {
                R_CheckUserInterrupt();
            }
        }
    }

    if (config->verbose) {
        Rprintf("    Total records: %zu\n", total_records);
        if (duplicates_skipped > 0) {
            Rprintf("    Duplicates skipped: %zu\n", duplicates_skipped);
        }
        Rprintf("    Unique sites: %zu\n", ht->n_entries);
    }

    // Transfer sites from hash table to coordinate set
    size_t n_unique = ht->n_entries;  // Save before freeing
    size_t chr_start_idx = cs->n_sites;  // Starting index for this chromosome

    for (int i = 0; i < ht->size; i++) {
        hash_site_entry_t *entry = ht->buckets[i];
        while (entry) {
            if (add_site_to_set(cs, &entry->site) != 0) {
                Rprintf("Error: Failed to add site to coordinate set\n");
                free_site_hash_table(ht);
                return -1;
            }
            entry = entry->next;
        }
    }

    free_site_hash_table(ht);

    // Sort the sites we just added (they come from hash table in random order)
    // Sort by position, strand, mod_code for consistent ordering
    qsort(&cs->sites[chr_start_idx], n_unique, sizeof(site_key_t),
          (int(*)(const void*, const void*))compare_sites);

    return 0;
}

// ============================================================================
// MAIN COORDINATE DISCOVERY
// ============================================================================

coordinate_set_t* discover_coordinates(char **input_files,
                                      int n_files,
                                      chrom_sizes_t *chrom_sizes,
                                      discovery_config_t *config) {

    if (config->verbose) {
        Rprintf("Phase 1: Coordinate Discovery\n");
        Rprintf("  Files: %d\n", n_files);
        Rprintf("  Target modification: %c\n", config->target_mod);
        Rprintf("  Chromosomes to process: %d\n", chrom_sizes->n_chroms);
    }

    // Initialize coordinate set
    coordinate_set_t *cs = init_coordinate_set(INITIAL_SITES_CAPACITY);
    if (!cs) {
        Rprintf("Error: Failed to initialize coordinate set\n");
        return NULL;
    }

    // Initialize discovery streams for all files
    discovery_stream_t *streams = (discovery_stream_t*)calloc(n_files,
                                                              sizeof(discovery_stream_t));
    if (!streams) {
        free_coordinate_set(cs);
        return NULL;
    }

    for (int i = 0; i < n_files; i++) {
        if (init_discovery_stream(&streams[i], input_files[i], i) != 0) {
            // Cleanup and error
            for (int j = 0; j < i; j++) {
                cleanup_discovery_stream(&streams[j]);
            }
            free(streams);
            free_coordinate_set(cs);
            return NULL;
        }
    }

    // Process each chromosome from chrom_sizes (already filtered to union of all files)
    // Since tabix files are sorted, we process chromosome-by-chromosome
    for (int c = 0; c < chrom_sizes->n_chroms; c++) {
        const char *chr = chrom_sizes->chroms[c].name;

        if (config->verbose) {
            Rprintf("  Processing %s (length: %d bp)...\n",
                    chr, chrom_sizes->chroms[c].length);
        }

        if (discover_chromosome_sites(cs, streams, n_files, chr, config) != 0) {
            goto cleanup;
        }
    }

    // Build chromosome index for fast lookup
    // Note: No post-discovery deduplication needed - HashMap guarantees uniqueness!
    if (build_chromosome_index(cs) != 0) {
        Rprintf("Error: Failed to build chromosome index\n");
        goto cleanup;
    }

    if (config->verbose) {
        Rprintf("  Total unique sites: %zu\n", cs->n_sites);
        Rprintf("  Chromosomes: %d\n", cs->n_chromosomes);
    }

    // Cleanup streams
    for (int i = 0; i < n_files; i++) {
        cleanup_discovery_stream(&streams[i]);
    }
    free(streams);

    return cs;

cleanup:
    for (int i = 0; i < n_files; i++) {
        cleanup_discovery_stream(&streams[i]);
    }
    free(streams);
    free_coordinate_set(cs);
    return NULL;
}