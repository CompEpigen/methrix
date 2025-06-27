#ifndef MODKIT_R_INTERFACE_H
#define MODKIT_R_INTERFACE_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>
#include <stdint.h>

// Include the main header after R headers to avoid conflicts
#include "modkit_merge_improved.h"

// Result structure for returning data to R
typedef struct {
  // Coordinate vectors
  char **chr_names;      // chromosome names
  int32_t *positions;    // start positions
  char *strands;         // strand characters
  int n_sites;           // number of sites
  
  // NEW: Reference context fields
  char *ref_bases;       // Reference base at each position
  char **contexts;       // CG/CHG/CHH context strings
  int include_context;   // Flag indicating context data available
  
  // Data matrices (stored as column-major for R)
  double *beta_matrix;   // methylation values
  int32_t *cov_matrix;   // coverage values
  int n_samples;         // number of samples
  
  // Sample names
  char **sample_names;
  
  // Chromosome summary for metadata
  char **chr_list;
  int32_t *chr_counts;
  int n_chromosomes;
  
  // Processing info
  int total_sites_processed;
  char **processing_notes;
  int n_notes;
  
} modkit_r_result_t;

// Function prototypes
SEXP read_modkit_c(SEXP files, SEXP target_mod, SEXP bin_size, 
                   SEXP n_cores, SEXP quality_filter, SEXP min_coverage,
                   SEXP combine_strands, SEXP ref_fasta, SEXP verbose);

// Helper functions
modkit_r_result_t* init_modkit_result(int n_sites, int n_samples);
void cleanup_modkit_result(modkit_r_result_t *result);
SEXP convert_modkit_result_to_r(modkit_r_result_t *result);
config_t parse_r_config(SEXP files, SEXP target_mod, SEXP bin_size, 
                        SEXP n_cores, SEXP quality_filter, SEXP min_coverage,
                        SEXP combine_strands, SEXP ref_fasta, SEXP verbose);

// NEW: Context extraction function prototypes
void accumulate_coordinate_row(const coord_t *coord, file_stream_t *streams, int n_files,
                               config_t *config);
int process_all_bins_to_memory(config_t *config);
int process_single_bin_to_memory(genomic_bin_t *bin, char **input_files, int n_files,
                                 config_t *config);
void reallocate_result_if_needed(void);
int estimate_total_sites(config_t *config);

#endif // MODKIT_R_INTERFACE_H