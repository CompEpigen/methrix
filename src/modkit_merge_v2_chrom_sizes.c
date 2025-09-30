#include "modkit_merge_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// ============================================================================
// CHROMOSOME SIZES FILE PARSER
// ============================================================================
// Reads TSV file with format: <chrom><TAB><length>
// Also compatible with .fai files (uses first 2 columns only)
// ============================================================================

chrom_sizes_t* read_chrom_sizes(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        Rprintf("Error: Cannot open chromosome sizes file: %s\n", filename);
        return NULL;
    }

    chrom_sizes_t *sizes = (chrom_sizes_t*)malloc(sizeof(chrom_sizes_t));
    if (!sizes) {
        fclose(fp);
        return NULL;
    }

    sizes->capacity = MAX_CHROMOSOMES;
    sizes->chroms = (chrom_size_t*)malloc(sizes->capacity * sizeof(chrom_size_t));
    if (!sizes->chroms) {
        free(sizes);
        fclose(fp);
        return NULL;
    }

    sizes->n_chroms = 0;

    char line[1024];
    int line_num = 0;

    while (fgets(line, sizeof(line), fp)) {
        line_num++;

        // Skip empty lines and comments
        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0') {
            continue;
        }

        // Parse line: chrom<TAB>length
        char chrom_name[MAX_CHR_LEN];
        long long length_ll;

        // Use scanf to parse
        int n_parsed = sscanf(line, "%63s\t%lld", chrom_name, &length_ll);

        if (n_parsed != 2) {
            Rprintf("Warning: Skipping malformed line %d in %s\n", line_num, filename);
            continue;
        }

        if (length_ll <= 0 || length_ll > INT32_MAX) {
            Rprintf("Warning: Invalid chromosome length on line %d: %lld\n",
                    line_num, length_ll);
            continue;
        }

        // Check capacity
        if (sizes->n_chroms >= sizes->capacity) {
            Rprintf("Error: Too many chromosomes (max %d)\n", MAX_CHROMOSOMES);
            break;
        }

        // Store chromosome
        chrom_size_t *chr = &sizes->chroms[sizes->n_chroms];
        strncpy(chr->name, chrom_name, MAX_CHR_LEN - 1);
        chr->name[MAX_CHR_LEN - 1] = '\0';
        chr->length = (int32_t)length_ll;
        chr->tid = -1;  // Will be filled later when matching with tabix

        sizes->n_chroms++;
    }

    fclose(fp);

    if (sizes->n_chroms == 0) {
        Rprintf("Error: No valid chromosomes found in %s\n", filename);
        free(sizes->chroms);
        free(sizes);
        return NULL;
    }

    return sizes;
}

void free_chrom_sizes(chrom_sizes_t *sizes) {
    if (!sizes) return;
    if (sizes->chroms) free(sizes->chroms);
    free(sizes);
}

// ============================================================================
// CHROMOSOME UNION HELPER
// ============================================================================
// Scans ALL input files to build union of chromosomes present in tabix indexes
// Filters chromosome sizes to only those present in at least one file
// ============================================================================

typedef struct {
    char name[MAX_CHR_LEN];
    int present;
} chrom_union_t;

int find_or_add_chrom(chrom_union_t *union_set, int *n_union, int capacity,
                      const char *chr_name) {
    // Search for existing
    for (int i = 0; i < *n_union; i++) {
        if (strcmp(union_set[i].name, chr_name) == 0) {
            union_set[i].present = 1;
            return i;
        }
    }

    // Add new
    if (*n_union >= capacity) {
        return -1;  // Overflow
    }

    strncpy(union_set[*n_union].name, chr_name, MAX_CHR_LEN - 1);
    union_set[*n_union].name[MAX_CHR_LEN - 1] = '\0';
    union_set[*n_union].present = 1;
    (*n_union)++;
    return *n_union - 1;
}

chrom_sizes_t* filter_chrom_sizes_by_tabix(chrom_sizes_t *all_sizes,
                                           char **input_files,
                                           int n_files,
                                           int verbose) {
    // Build union of chromosomes from ALL tabix files
    chrom_union_t *tabix_union = (chrom_union_t*)calloc(MAX_CHROMOSOMES,
                                                         sizeof(chrom_union_t));
    if (!tabix_union) {
        return NULL;
    }

    int n_tabix_chroms = 0;

    if (verbose) {
        Rprintf("Scanning tabix files for chromosome union...\n");
    }

    for (int i = 0; i < n_files; i++) {
        htsFile *fp = hts_open(input_files[i], "r");
        if (!fp) {
            Rprintf("Warning: Cannot open %s for chromosome scan\n", input_files[i]);
            continue;
        }

        tbx_t *idx = tbx_index_load(input_files[i]);
        if (!idx) {
            Rprintf("Warning: Cannot load tabix index for %s\n", input_files[i]);
            hts_close(fp);
            continue;
        }

        // Get chromosome names from this file
        int n_seqs = 0;
        const char **seqnames = tbx_seqnames(idx, &n_seqs);

        if (seqnames) {
            for (int j = 0; j < n_seqs; j++) {
                find_or_add_chrom(tabix_union, &n_tabix_chroms,
                                 MAX_CHROMOSOMES, seqnames[j]);
            }
            free(seqnames);
        }

        tbx_destroy(idx);
        hts_close(fp);
    }

    if (verbose) {
        Rprintf("  Found %d unique chromosome(s) across all files\n", n_tabix_chroms);
    }

    // Filter chromosome sizes to only those in tabix union
    chrom_sizes_t *filtered = (chrom_sizes_t*)malloc(sizeof(chrom_sizes_t));
    if (!filtered) {
        free(tabix_union);
        return NULL;
    }

    filtered->capacity = n_tabix_chroms;
    filtered->chroms = (chrom_size_t*)malloc(n_tabix_chroms * sizeof(chrom_size_t));
    if (!filtered->chroms) {
        free(filtered);
        free(tabix_union);
        return NULL;
    }

    filtered->n_chroms = 0;

    // Match sizes with tabix chroms
    for (int i = 0; i < n_tabix_chroms; i++) {
        const char *tabix_chr = tabix_union[i].name;

        // Find this chromosome in sizes file
        int found = 0;
        for (int j = 0; j < all_sizes->n_chroms; j++) {
            if (strcmp(all_sizes->chroms[j].name, tabix_chr) == 0) {
                // Match - include this chromosome
                filtered->chroms[filtered->n_chroms] = all_sizes->chroms[j];
                filtered->chroms[filtered->n_chroms].tid = i;
                filtered->n_chroms++;
                found = 1;
                break;
            }
        }

        if (!found && verbose) {
            Rprintf("  Warning: Chromosome '%s' in tabix but not in sizes file (skipped)\n",
                    tabix_chr);
        }
    }

    free(tabix_union);

    if (filtered->n_chroms == 0) {
        Rprintf("Error: No chromosomes match between sizes file and tabix indexes\n");
        free(filtered->chroms);
        free(filtered);
        return NULL;
    }

    if (verbose) {
        Rprintf("  Processing %d chromosome(s) with known sizes\n", filtered->n_chroms);
    }

    return filtered;
}