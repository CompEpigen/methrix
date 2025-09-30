#ifndef MODKIT_MERGE_IMPROVED_H
#define MODKIT_MERGE_IMPROVED_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>

// Include HTS libraries
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

// For R interface - include after standard headers
#ifndef R_R_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#endif

#define MAX_FILES 1000
#define MAX_CHR_LEN 64
#define MAX_BINS 100000
#define DEFAULT_BIN_SIZE 1000000  // 1MB bins
#define MAX_CHROMOSOMES 1000

// ModKit bedMethyl record structure
typedef struct {
  char chr[MAX_CHR_LEN];
  int32_t start;
  int32_t end;
  char mod_code;
  char strand;
  int32_t nvalid_cov;
  double frac_mod;
  int32_t nmod;
  int32_t ncanonical;
  int32_t nother_mod;
  int32_t nfail;
  int32_t nnocall;
} modkit_record_t;

// Coordinate structure
typedef struct {
  char chr[MAX_CHR_LEN];
  int32_t start;
  char strand;
} coord_t;

// Chromosome information with data range
typedef struct {
  char name[MAX_CHR_LEN];
  int32_t min_pos;
  int32_t max_pos;
  int has_data;
} chromosome_t;

// NEW: Chromosome length structure for fixed binning
typedef struct {
  char name[MAX_CHR_LEN];
  int64_t length;               // True chromosome length from .fai
} chromosome_length_t;

// NEW: Coordinate set for union operations
typedef struct {
  coord_t *coords;
  int n_coords;
  int capacity;
} coordinate_set_t;

// Genomic bin for processing
typedef struct {
  char chr[MAX_CHR_LEN];
  int32_t start;
  int32_t end;
  int bin_id;
} genomic_bin_t;

// Stream state for each input file
typedef struct {
  int file_idx;
  char *filename;
  htsFile *hts_file;
  tbx_t *tbx_handle;
  hts_itr_t *iterator;
  kstring_t line_buffer;
  
  int has_current_record;
  modkit_record_t current_record;
  
  int has_output_data;
  double output_beta;
  int32_t output_cov;

  // NEW: For strand aggregation
  int has_pending_record;
  modkit_record_t pending_record;
} file_stream_t;

// NEW: Structure for aggregating strand data
typedef struct {
  int32_t total_meth_calls;     // Sum of nmod from both strands
  int32_t total_canonical_calls; // Sum of ncanonical from both strands
  int32_t total_coverage;       // Sum of nvalid_cov from both strands
  double combined_beta;         // Calculated combined beta value
  int has_plus_strand;          // Flag: has plus strand data
  int has_minus_strand;         // Flag: has minus strand data
  int record_count;             // Number of records aggregated
} strand_aggregate_t;

// Reference sequence cache for performance
typedef struct {
  char chr[MAX_CHR_LEN];
  int32_t start;
  int32_t end;
  char *sequence;
  int valid;
} ref_cache_t;

// Configuration structure
typedef struct {
  char **input_files;
  int n_files;
  char *region;
  char target_mod;
  double quality_filter;
  int min_coverage;
  char *output_file;
  int combine_strands;
  int verbose;
  int bin_size;
  int n_threads;
  
  // Reference FASTA options
  char *ref_fasta;
  int include_context;
  int validate_cytosines;
  
  // Auto-discovered data
  chromosome_t *chromosomes;
  int n_chromosomes;
  genomic_bin_t *bins;
  int n_bins;

  // NEW: Fixed binning support
  int use_fixed_bins;           // Flag for fixed vs data-driven bins
  chromosome_length_t *chr_lengths;  // True chromosome lengths
  int n_chr_lengths;            // Number of chromosomes with known lengths
} config_t;

// Global reference handle and cache
extern faidx_t *ref_fai;
extern ref_cache_t ref_cache;
extern const int CACHE_SIZE;

// Core algorithm functions
int discover_chromosomes(config_t *config);
int discover_chromosome_from_file(const char *filename, chromosome_t *chromosomes, int *n_chromosomes);
void merge_chromosome_ranges(chromosome_t *chromosomes, int *n_chromosomes, 
                             const char *chr_name, int32_t min_pos, int32_t max_pos);
int generate_bins(config_t *config);

// File stream management
int init_file_stream_for_bin(file_stream_t *stream, const char *filename, int file_idx,
                             genomic_bin_t *bin, config_t *config);
int advance_stream(file_stream_t *stream, config_t *config);
void cleanup_streams(file_stream_t *streams, int n_files);

// Coordinate processing
coord_t find_minimum_coordinate(file_stream_t *streams, int n_files);
coord_t find_minimum_coordinate_with_merging(file_stream_t *streams, int n_files, config_t *config);
int coordinate_compare(const coord_t *a, const coord_t *b);
int coordinates_equal(const coord_t *a, const coord_t *b);
int any_stream_has_data(file_stream_t *streams, int n_files);
void collect_coordinate_data(file_stream_t *streams, int n_files, const coord_t *target_coord, config_t *config);
void reset_output_data(file_stream_t *streams, int n_files);

// NEW: Strand merging functions
coord_t normalize_cpg_coordinate(const modkit_record_t *record, config_t *config);
int cpg_coordinates_equal(const coord_t *a, const coord_t *b);
void init_strand_aggregate(strand_aggregate_t *agg);
void add_record_to_aggregate(strand_aggregate_t *agg, const modkit_record_t *record);
void finalize_strand_aggregate(strand_aggregate_t *agg);
void collect_coordinate_data_with_merging(file_stream_t *streams, int n_files, const coord_t *target_coord, config_t *config);

// NEW: Fixed binning functions
int extract_chromosomes_from_bedmethyl_files(config_t *config);
int extract_chromosome_lengths_from_fai(config_t *config);
int generate_fixed_genomic_bins(config_t *config);
coordinate_set_t collect_all_coordinates_in_bin(char **input_files, int n_files, genomic_bin_t *bin, config_t *config);
int process_bin_with_complete_union(genomic_bin_t *bin, char **input_files, int n_files, config_t *config);

void init_coordinate_set(coordinate_set_t *coord_set, int initial_capacity);
void add_coordinate_to_set(coordinate_set_t *coord_set, const coord_t *coord);
void sort_coordinate_set(coordinate_set_t *coord_set);
void cleanup_coordinate_set(coordinate_set_t *coord_set);
int query_file_for_coordinate(const char *filename, const coord_t *target_coord,
                             file_stream_t *result_stream, config_t *config);

// Record processing
int parse_modkit_record(char *line, modkit_record_t *record);
int apply_record_filters(modkit_record_t *record, config_t *config);
void apply_coordinate_transforms(modkit_record_t *record, config_t *config);

// Reference genome functions
int init_reference(config_t *config);
void cleanup_reference();
char get_ref_base_cached(const char *chr, int32_t pos);
char* get_context(const char *chr, int32_t pos, char strand);
char complement_base(char base);
void reverse_complement_inplace(char *seq, int len);

// Utility functions
char* format_region_string(genomic_bin_t *bin);
void cleanup_config(config_t *config);

#endif // MODKIT_MERGE_IMPROVED_H