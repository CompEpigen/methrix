#ifndef MODKIT_MERGE_V2_H
#define MODKIT_MERGE_V2_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>
#include <stdint.h>
#include <htslib/tbx.h>
#include <htslib/hts.h>

// Maximum dimensions
#define MAX_CHR_LEN 64
#define MAX_FILES 1000
#define INITIAL_SITES_CAPACITY 1000000
#define MAX_CHROMOSOMES 1000

// ============================================================================
// CHROMOSOME SIZES
// ============================================================================

// Chromosome size record
typedef struct {
    char name[MAX_CHR_LEN];
    int32_t length;         // Chromosome length in bp
    int32_t tid;            // Tabix chromosome ID (-1 if not in tabix)
} chrom_size_t;

// Chromosome sizes from file
typedef struct {
    chrom_size_t *chroms;
    int n_chroms;
    int capacity;
} chrom_sizes_t;

// Reference cache structure
typedef struct {
    char chr[MAX_CHR_LEN];
    int32_t start;
    int32_t end;
    char *sequence;
    int valid;
} ref_cache_t;

// Function prototypes for chromosome sizes
chrom_sizes_t* read_chrom_sizes(const char *filename);
void free_chrom_sizes(chrom_sizes_t *sizes);
chrom_sizes_t* filter_chrom_sizes_by_tabix(chrom_sizes_t *all_sizes,
                                           char **input_files,
                                           int n_files,
                                           int verbose);

// ============================================================================
// PHASE 1: COORDINATE DISCOVERY
// ============================================================================

// Site key: uniquely identifies a genomic position
typedef struct {
    char chr[MAX_CHR_LEN];
    int32_t pos;          // 0-based position (bedMethyl standard)
    char strand;          // '+', '-', or '.' for both
    char mod_code;        // 'm', 'h', 'a', etc.
} site_key_t;

// Coordinate set: sorted list of all discovered sites
typedef struct {
    site_key_t *sites;    // Array of sites (sorted)
    size_t n_sites;       // Current count
    size_t capacity;      // Allocated capacity

    // Chromosome boundary indices for fast lookup
    int32_t *chr_start_idx;  // Start index for each chromosome
    int32_t *chr_end_idx;    // End index (exclusive) for each chromosome
    int n_chromosomes;       // Number of distinct chromosomes
    char **chr_names;        // Chromosome names
} coordinate_set_t;

// Per-file tabix stream for coordinate discovery
typedef struct {
    htsFile *fp;
    tbx_t *idx;
    hts_itr_t *iter;
    kstring_t line;

    site_key_t current_site;  // Current site from this file
    int has_site;             // 1 if current_site is valid
    int file_idx;             // Index of this file
    char *filename;
} discovery_stream_t;

// Configuration for coordinate discovery
typedef struct {
    char target_mod;          // Target modification code
    int combine_strands;      // Combine + and - strands?
    int min_coverage;         // Minimum coverage filter
    double quality_filter;    // Maximum failure rate
    int verbose;              // Print progress?
} discovery_config_t;

// Function prototypes for coordinate discovery
coordinate_set_t* discover_coordinates(
    char **input_files,
    int n_files,
    chrom_sizes_t *chrom_sizes,
    discovery_config_t *config
);

void free_coordinate_set(coordinate_set_t *coord_set);

int compare_sites(const site_key_t *a, const site_key_t *b);

// Hash table structures (used by both discovery.c and streaming.c)
// Forward declare the key type
typedef struct {
    int32_t pos;
    char strand;
} site_hash_key_t;

typedef struct hash_site_entry {
    site_hash_key_t key;
    site_key_t site;
    struct hash_site_entry *next;
} hash_site_entry_t;

typedef struct site_hash_table {
    hash_site_entry_t **buckets;
    int size;
    size_t n_entries;
} site_hash_table_t;

#define DISCOVERY_HASH_SIZE 2097152

site_hash_table_t* create_site_hash_table(int size);
void free_site_hash_table(site_hash_table_t *ht);
int insert_unique_site(site_hash_table_t *ht, const site_key_t *site);

// Function prototypes for context extraction (after coordinate_set_t definition)
int init_reference_v2(const char *ref_fasta_path);
void cleanup_reference_v2();
char get_ref_base_cached_v2(const char *chr, int32_t pos);
char* get_context_v2(const char *chr, int32_t pos, char strand);
int extract_context_for_sites(coordinate_set_t *coord_set,
                              char **ref_base_out,
                              char **context_out,
                              int verbose);

// R interface for standalone context extraction
SEXP add_context_c(SEXP chr_vec, SEXP pos_vec, SEXP strand_vec, SEXP ref_fasta, SEXP verbose);

// R interface for streaming HDF5 support (old chunk-based)
SEXP read_modkit_v2_init_chunks_c(SEXP files, SEXP chrom_sizes_file, SEXP target_mod,
                                   SEXP interval_size, SEXP min_coverage,
                                   SEXP quality_filter, SEXP combine_strands, SEXP ref_fasta,
                                   SEXP verbose);
SEXP read_modkit_v2_process_chunk_c(SEXP coord_ptr, SEXP chunks_ptr, SEXP chunk_id,
                                     SEXP files, SEXP target_mod, SEXP min_coverage,
                                     SEXP quality_filter, SEXP combine_strands, SEXP verbose);
SEXP read_modkit_v2_cleanup_c(SEXP coord_ptr, SEXP chunks_ptr);

// R interface for true streaming (interval-based discovery)
SEXP read_modkit_v2_process_interval_c(SEXP chr, SEXP interval_start, SEXP interval_end,
                                        SEXP files, SEXP target_mod, SEXP min_coverage,
                                        SEXP quality_filter, SEXP combine_strands, SEXP verbose);

// ============================================================================
// PHASE 2: CHUNKED PROCESSING
// ============================================================================

// Genomic chunk: contiguous range of sites to process together
typedef struct {
    char chr[MAX_CHR_LEN];
    int32_t start_site_idx;   // Index into coordinate_set
    int32_t end_site_idx;     // Exclusive end index
    int32_t start_pos;        // Genomic position (for tabix query)
    int32_t end_pos;          // Genomic position
    int chunk_id;             // Sequential chunk ID
} genomic_chunk_t;

// Chunk list
typedef struct {
    genomic_chunk_t *chunks;
    int n_chunks;
    size_t sites_per_chunk;   // Target sites per chunk
} chunk_list_t;

// Chunk result: aggregated data for one chunk
typedef struct {
    double *beta_values;      // [n_sites × n_samples]
    int32_t *cov_values;      // [n_sites × n_samples]
    uint8_t *has_data_bits;   // Bitmap: (n_sites × n_samples + 7) / 8 bytes
    int32_t n_sites;          // Sites in this chunk
    int32_t n_samples;        // Number of samples
} chunk_result_t;

// Count data for aggregation
typedef struct {
    int32_t n_mod;            // Modified base calls
    int32_t n_canonical;      // Canonical base calls
    int32_t valid_cov;        // Valid coverage
} site_counts_t;

// Function prototypes for chunking
chunk_list_t* create_chunks(
    coordinate_set_t *coord_set,
    size_t sites_per_chunk,
    int32_t max_genomic_span
);

chunk_list_t* create_chunks_by_interval(
    coordinate_set_t *coord_set,
    chrom_sizes_t *chrom_sizes,
    int32_t interval_bp
);

void free_chunk_list(chunk_list_t *chunks);

chunk_result_t* alloc_chunk_result(int32_t n_sites, int32_t n_samples);

chunk_result_t* process_chunk(
    genomic_chunk_t *chunk,
    coordinate_set_t *coord_set,
    char **input_files,
    int n_files,
    discovery_config_t *config
);

void free_chunk_result(chunk_result_t *result);

// ============================================================================
// PHASE 3: MATRIX DENSIFICATION
// ============================================================================

// R result structure
typedef struct {
    SEXP beta_matrix;         // R matrix (REALSXP)
    SEXP cov_matrix;          // R matrix (INTSXP)
    SEXP chr_vec;             // Chromosome names
    SEXP pos_vec;             // Positions
    SEXP strand_vec;          // Strands
    SEXP mod_code_vec;        // Modification codes
    int n_sites;
    int n_samples;
} r_result_t;

r_result_t* create_r_result(
    coordinate_set_t *coord_set,
    int n_samples
);

void populate_r_matrices(
    r_result_t *r_result,
    chunk_list_t *chunks,
    coordinate_set_t *coord_set
);

// ============================================================================
// UTILITIES
// ============================================================================

// Parse bedMethyl line into site_key and counts
int parse_bedmethyl_line(
    const char *line,
    site_key_t *site,
    site_counts_t *counts
);

// Binary search in sorted site array
int binary_search_site(
    const site_key_t *sites,
    int start_idx,
    int end_idx,
    const site_key_t *target
);

// Check if counts pass quality filters
int passes_filters(
    site_counts_t *counts,
    discovery_config_t *config
);

// Compute beta value from counts
double compute_beta(site_counts_t *counts);

// Bitmap operations
static inline void set_bit(uint8_t *bits, size_t idx) {
    bits[idx / 8] |= (1 << (idx % 8));
}

static inline int get_bit(const uint8_t *bits, size_t idx) {
    return (bits[idx / 8] >> (idx % 8)) & 1;
}

#endif // MODKIT_MERGE_V2_H