#include "modkit_merge_v2.h"
#include <ctype.h>
#include <htslib/faidx.h>

// ============================================================================
// REFERENCE FASTA AND CONTEXT EXTRACTION
// ============================================================================
// Reusable functions for extracting reference base and sequence context
// Used by both read_modkit_v2() and add_context()
// ============================================================================

// Global reference handle and cache for performance
static faidx_t *ref_fai_v2 = NULL;
static ref_cache_t ref_cache_v2 = {0};
static const int CACHE_SIZE_V2 = 50000;  // 50KB cache chunks

// ============================================================================
// REFERENCE INITIALIZATION
// ============================================================================

int init_reference_v2(const char *ref_fasta_path) {
    if (!ref_fasta_path || strlen(ref_fasta_path) == 0) {
        return 0;  // No reference provided
    }

    ref_fai_v2 = fai_load(ref_fasta_path);
    if (!ref_fai_v2) {
        Rprintf("Error: Cannot load FASTA index for %s\n", ref_fasta_path);
        Rprintf("Make sure %s.fai exists or run: samtools faidx %s\n",
                ref_fasta_path, ref_fasta_path);
        return -1;
    }

    // Initialize cache
    ref_cache_v2.valid = 0;
    ref_cache_v2.sequence = NULL;

    return 0;
}

void cleanup_reference_v2() {
    if (ref_fai_v2) {
        fai_destroy(ref_fai_v2);
        ref_fai_v2 = NULL;
    }

    if (ref_cache_v2.sequence) {
        free(ref_cache_v2.sequence);
        ref_cache_v2.sequence = NULL;
    }
    ref_cache_v2.valid = 0;
}

// ============================================================================
// CACHED REFERENCE BASE RETRIEVAL
// ============================================================================

char get_ref_base_cached_v2(const char *chr, int32_t pos) {
    if (!ref_fai_v2) {
        return 'N';
    }

    // Check if position is in cache
    if (ref_cache_v2.valid &&
        strcmp(ref_cache_v2.chr, chr) == 0 &&
        pos >= ref_cache_v2.start &&
        pos < ref_cache_v2.end) {
        return toupper(ref_cache_v2.sequence[pos - ref_cache_v2.start]);
    }

    // Update cache with new chunk
    int cache_start = (pos / CACHE_SIZE_V2) * CACHE_SIZE_V2;
    int cache_end = cache_start + CACHE_SIZE_V2;

    if (ref_cache_v2.sequence) {
        free(ref_cache_v2.sequence);
    }

    int len;
    ref_cache_v2.sequence = faidx_fetch_seq(ref_fai_v2, chr, cache_start, cache_end - 1, &len);
    if (!ref_cache_v2.sequence) {
        ref_cache_v2.valid = 0;
        return 'N';
    }

    strncpy(ref_cache_v2.chr, chr, MAX_CHR_LEN - 1);
    ref_cache_v2.chr[MAX_CHR_LEN - 1] = '\0';
    ref_cache_v2.start = cache_start;
    ref_cache_v2.end = cache_start + len;
    ref_cache_v2.valid = 1;

    // Return the requested base
    if (pos >= ref_cache_v2.start && pos < ref_cache_v2.end) {
        return toupper(ref_cache_v2.sequence[pos - ref_cache_v2.start]);
    }

    return 'N';
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

static char complement_base(char base) {
    switch (toupper(base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

static void reverse_complement_inplace(char *seq, int len) {
    // Reverse the sequence
    for (int i = 0; i < len / 2; i++) {
        char temp = seq[i];
        seq[i] = seq[len - 1 - i];
        seq[len - 1 - i] = temp;
    }

    // Complement each base
    for (int i = 0; i < len; i++) {
        seq[i] = complement_base(seq[i]);
    }
}

// ============================================================================
// CONTEXT EXTRACTION
// ============================================================================

char* get_context_v2(const char *chr, int32_t pos, char strand) {
    if (!ref_fai_v2) {
        char *result = (char*)malloc(4);
        strcpy(result, "NNN");
        return result;
    }

    // Get 3-base context (position + 2 downstream bases)
    int len;
    char *seq = faidx_fetch_seq(ref_fai_v2, chr, pos, pos + 2, &len);
    if (!seq || len < 3) {
        if (seq) free(seq);
        char *result = (char*)malloc(4);
        strcpy(result, "NNN");
        return result;
    }

    // Convert to uppercase
    for (int i = 0; i < len; i++) {
        seq[i] = toupper(seq[i]);
    }

    // For minus strand, get reverse complement
    if (strand == '-') {
        reverse_complement_inplace(seq, 3);
    }

    // Determine context
    char *context = (char*)malloc(4);
    if (len >= 3) {
        if (seq[0] == 'C' && seq[1] == 'G') {
            strcpy(context, "CG");
        } else if (seq[0] == 'C' && seq[2] == 'G') {
            strcpy(context, "CHG");
        } else if (seq[0] == 'C') {
            strcpy(context, "CHH");
        } else {
            strcpy(context, "CNN");
        }
    } else {
        strcpy(context, "NNN");
    }

    free(seq);
    return context;
}

// ============================================================================
// BATCH CONTEXT EXTRACTION FOR COORDINATE SET
// ============================================================================

int extract_context_for_sites(coordinate_set_t *coord_set,
                              char **ref_base_out,
                              char **context_out,
                              int verbose) {

    if (!coord_set || coord_set->n_sites == 0) {
        return -1;
    }

    if (!ref_fai_v2) {
        Rprintf("Error: Reference FASTA not initialized\n");
        return -1;
    }

    if (verbose) {
        Rprintf("\nPhase 5: Extracting Reference Context\n");
        Rprintf("  Total sites: %zu\n", coord_set->n_sites);
    }

    // Allocate output arrays
    *ref_base_out = (char*)malloc(coord_set->n_sites * sizeof(char));
    *context_out = (char*)malloc(coord_set->n_sites * 4 * sizeof(char));  // 3 chars + null per site

    if (!*ref_base_out || !*context_out) {
        if (*ref_base_out) free(*ref_base_out);
        if (*context_out) free(*context_out);
        return -1;
    }

    // Extract context for each site
    const char *current_chr = NULL;
    size_t sites_processed = 0;

    for (size_t i = 0; i < coord_set->n_sites; i++) {
        site_key_t *site = &coord_set->sites[i];

        // Progress reporting per chromosome
        if (current_chr == NULL || strcmp(current_chr, site->chr) != 0) {
            if (verbose && current_chr != NULL) {
                Rprintf("    Processed %zu sites on %s\n", sites_processed, current_chr);
                sites_processed = 0;
            }
            current_chr = site->chr;
            if (verbose) {
                Rprintf("  Processing %s...\n", current_chr);
            }
        }

        // Get reference base
        (*ref_base_out)[i] = get_ref_base_cached_v2(site->chr, site->pos);

        // Get context
        char *ctx = get_context_v2(site->chr, site->pos, site->strand);
        strncpy(&(*context_out)[i * 4], ctx, 3);
        (*context_out)[i * 4 + 3] = '\0';
        free(ctx);

        sites_processed++;

        // Periodic interrupt check
        if (i % 10000 == 0 && i > 0) {
            R_CheckUserInterrupt();
        }
    }

    if (verbose) {
        Rprintf("    Processed %zu sites on %s\n", sites_processed, current_chr);
        Rprintf("  Context extraction complete\n");
    }

    return 0;
}

// ============================================================================
// R INTERFACE FOR STANDALONE CONTEXT EXTRACTION
// ============================================================================

SEXP add_context_c(SEXP chr_vec, SEXP pos_vec, SEXP strand_vec, SEXP ref_fasta, SEXP verbose) {
    // Parse inputs
    int n_sites = length(chr_vec);
    const char *ref_fasta_path = CHAR(STRING_ELT(ref_fasta, 0));
    int verb = asLogical(verbose);

    if (n_sites == 0) {
        error("No sites provided");
    }

    // Initialize reference
    if (init_reference_v2(ref_fasta_path) != 0) {
        error("Failed to initialize reference FASTA: %s", ref_fasta_path);
    }

    if (verb) {
        Rprintf("Extracting context for %d sites from %s\n", n_sites, ref_fasta_path);
    }

    // Allocate output vectors
    SEXP ref_base_vec = PROTECT(allocVector(STRSXP, n_sites));
    SEXP context_vec = PROTECT(allocVector(STRSXP, n_sites));

    // Extract context for each site
    const char *current_chr = NULL;
    int sites_processed = 0;

    for (int i = 0; i < n_sites; i++) {
        const char *chr = CHAR(STRING_ELT(chr_vec, i));
        int32_t pos = INTEGER(pos_vec)[i];
        const char *strand_str = CHAR(STRING_ELT(strand_vec, i));
        char strand = strand_str[0];

        // Progress reporting per chromosome
        if (current_chr == NULL || strcmp(current_chr, chr) != 0) {
            if (verb && current_chr != NULL) {
                Rprintf("  Processed %d sites on %s\n", sites_processed, current_chr);
                sites_processed = 0;
            }
            current_chr = chr;
            if (verb) {
                Rprintf("Processing %s...\n", chr);
            }
        }

        // Get reference base
        char ref_base = get_ref_base_cached_v2(chr, pos);
        char ref_str[2] = {ref_base, '\0'};
        SET_STRING_ELT(ref_base_vec, i, mkChar(ref_str));

        // Get context
        char *ctx = get_context_v2(chr, pos, strand);
        SET_STRING_ELT(context_vec, i, mkChar(ctx));
        free(ctx);

        sites_processed++;

        // Periodic interrupt check
        if (i % 10000 == 0 && i > 0) {
            R_CheckUserInterrupt();
        }
    }

    if (verb) {
        Rprintf("  Processed %d sites on %s\n", sites_processed, current_chr);
        Rprintf("Context extraction complete\n");
    }

    // Cleanup
    cleanup_reference_v2();

    // Create result list
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, ref_base_vec);
    SET_VECTOR_ELT(result, 1, context_vec);

    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("ref_base"));
    SET_STRING_ELT(names, 1, mkChar("context"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(4);
    return result;
}