#include "modkit_merge_v2.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

// ============================================================================
// BEDMETHYL PARSING
// ============================================================================

int parse_bedmethyl_line(const char *line,
                        site_key_t *site,
                        site_counts_t *counts) {
    // bedMethyl format:
    // chr start end mod_code score strand thickStart thickEnd itemRgb
    // coverage percent_modified n_mod n_canonical n_other_mod n_delete n_fail n_diff n_nocall

    char chr[MAX_CHR_LEN];
    int32_t start, end;
    char mod_code_str[32];
    int score;
    char strand_char;
    int thick_start, thick_end;
    char item_rgb[32];
    int coverage;
    float percent_modified;
    int n_mod, n_canonical, n_other, n_delete, n_fail, n_diff, n_nocall;

    int ret = sscanf(line, "%63s %d %d %31s %d %c %d %d %31s %d %f %d %d %d %d %d %d %d",
                    chr, &start, &end, mod_code_str, &score, &strand_char,
                    &thick_start, &thick_end, item_rgb,
                    &coverage, &percent_modified,
                    &n_mod, &n_canonical, &n_other,
                    &n_delete, &n_fail, &n_diff, &n_nocall);

    if (ret < 18) {
        return -1;  // Parse failed
    }

    // Fill site key
    strncpy(site->chr, chr, MAX_CHR_LEN - 1);
    site->chr[MAX_CHR_LEN - 1] = '\0';
    site->pos = start;  // bedMethyl is 0-based
    site->strand = strand_char;

    // Parse modification code
    // Handle comma-separated codes like "m,h" - take first
    char *comma = strchr(mod_code_str, ',');
    if (comma) *comma = '\0';

    if (strlen(mod_code_str) > 0) {
        site->mod_code = mod_code_str[0];
    } else {
        return -1;
    }

    // Fill counts
    counts->n_mod = n_mod;
    counts->n_canonical = n_canonical;
    counts->valid_cov = coverage;

    return 0;
}

// ============================================================================
// FILTERING
// ============================================================================

int passes_filters(site_counts_t *counts, discovery_config_t *config) {
    // Coverage filter
    if (counts->valid_cov < config->min_coverage) {
        return 0;
    }

    // Quality filter (failure rate)
    if (config->quality_filter > 0.0) {
        // Note: We don't have n_fail from parse yet (would need to parse more fields)
        // For now, skip quality filtering in discovery phase
        // Can add back if needed
    }

    return 1;
}

// ============================================================================
// BETA COMPUTATION
// ============================================================================

double compute_beta(site_counts_t *counts) {
    int total = counts->n_mod + counts->n_canonical;
    if (total == 0) {
        return R_NaN;
    }
    return (double)counts->n_mod / (double)total;
}

// ============================================================================
// BINARY SEARCH
// ============================================================================

int binary_search_site(const site_key_t *sites,
                      int start_idx,
                      int end_idx,
                      const site_key_t *target) {
    int left = start_idx;
    int right = end_idx - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        int cmp = compare_sites(&sites[mid], target);

        if (cmp == 0) {
            return mid;  // Found
        } else if (cmp < 0) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return -1;  // Not found
}

// ============================================================================
// CHROMOSOME LOOKUP IN COORDINATE SET
// ============================================================================

int get_chromosome_bounds(coordinate_set_t *cs,
                         const char *chr,
                         int32_t *start_idx,
                         int32_t *end_idx) {
    // Linear search through chromosome names (small number)
    for (int i = 0; i < cs->n_chromosomes; i++) {
        if (strcmp(cs->chr_names[i], chr) == 0) {
            *start_idx = cs->chr_start_idx[i];
            *end_idx = cs->chr_end_idx[i];
            return 0;
        }
    }
    return -1;  // Chromosome not found
}