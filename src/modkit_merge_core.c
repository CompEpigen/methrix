#include "modkit_merge_improved.h"

// Global reference handle and cache
faidx_t *ref_fai = NULL;
ref_cache_t ref_cache = {0};
const int CACHE_SIZE = 50000;  // 50KB cache chunks

// Auto-discovery functions
int discover_chromosomes(config_t *config) {
  config->chromosomes = (chromosome_t*)R_alloc(MAX_CHROMOSOMES, sizeof(chromosome_t));
  config->n_chromosomes = 0;
  
  if (config->verbose) {
    Rprintf("Auto-discovering chromosomes from %d files...\n", config->n_files);
  }
  
  // Process each input file to discover chromosome ranges
  for (int i = 0; i < config->n_files; i++) {
    if (config->verbose) {
      Rprintf("  Scanning %s...\n", config->input_files[i]);
    }
    
    if (discover_chromosome_from_file(config->input_files[i], 
                                      config->chromosomes, &config->n_chromosomes) != 0) {
      if (config->verbose) {
        Rprintf("Warning: Failed to read index from %s\n", config->input_files[i]);
      }
      continue;
    }
  }
  
  if (config->n_chromosomes == 0) {
    error("No chromosomes discovered from input files");
  }
  
  // Sort chromosomes by name for consistent processing
  for (int i = 0; i < config->n_chromosomes - 1; i++) {
    for (int j = i + 1; j < config->n_chromosomes; j++) {
      if (strcmp(config->chromosomes[i].name, config->chromosomes[j].name) > 0) {
        chromosome_t temp = config->chromosomes[i];
        config->chromosomes[i] = config->chromosomes[j];
        config->chromosomes[j] = temp;
      }
    }
  }
  
  if (config->verbose) {
    Rprintf("Discovered %d chromosomes with data\n", config->n_chromosomes);
    for (int i = 0; i < config->n_chromosomes; i++) {
      Rprintf("  Chr %s: positions %d-%d\n", config->chromosomes[i].name,
              config->chromosomes[i].min_pos, config->chromosomes[i].max_pos);
    }
    Rprintf("  WARNING: Files may have different starting positions within this range\n");
  }
  
  return 0;
}

int discover_chromosome_from_file(const char *filename, chromosome_t *chromosomes, int *n_chromosomes) {
  htsFile *fp = hts_open(filename, "r");
  if (!fp) {
    return 1;
  }
  
  tbx_t *tbx = tbx_index_load(filename);
  if (!tbx) {
    hts_close(fp);
    return 1;
  }
  
  // Try common chromosome names to discover what's available
  const char *common_chrs[] = {
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM",
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y", "MT", "M"
  };
  
  for (size_t i = 0; i < sizeof(common_chrs) / sizeof(common_chrs[0]); i++) {
    hts_itr_t *iter = tbx_itr_querys(tbx, common_chrs[i]);
    if (!iter) {
      continue;  // No data for this chromosome
    }
    
    kstring_t line = {0, 0, 0};
    int32_t min_pos = INT32_MAX;
    int32_t max_pos = 0;
    int has_data = 0;
    int record_count = 0;
    
    // Sample some records to find the range
    while (tbx_itr_next(fp, tbx, iter, &line) >= 0 && record_count < 500) {
      char *tab1 = strchr(line.s, '\t');
      if (tab1) {
        int32_t pos = atoi(tab1 + 1);
        if (pos < min_pos) min_pos = pos;
        if (pos > max_pos) max_pos = pos;
        has_data = 1;
      }
      record_count++;
    }
    
    if (line.s) {
      free(line.s);
    }
    tbx_itr_destroy(iter);
    
    // If we found data, add/merge this chromosome
    if (has_data) {
      merge_chromosome_ranges(chromosomes, n_chromosomes, common_chrs[i], min_pos, max_pos);
    }
  }
  
  tbx_destroy(tbx);
  hts_close(fp);
  
  return 0;
}

void merge_chromosome_ranges(chromosome_t *chromosomes, int *n_chromosomes, 
                             const char *chr_name, int32_t min_pos, int32_t max_pos) {
  
  if (!chromosomes || !n_chromosomes || !chr_name) {
    return;
  }
  
  // Look for existing chromosome entry
  for (int i = 0; i < *n_chromosomes; i++) {
    if (strcmp(chromosomes[i].name, chr_name) == 0) {
      // Update existing entry with expanded range
      if (min_pos < chromosomes[i].min_pos) {
        chromosomes[i].min_pos = min_pos;
      }
      if (max_pos > chromosomes[i].max_pos) {
        chromosomes[i].max_pos = max_pos;
      }
      return;
    }
  }
  
  // Add new chromosome entry
  if (*n_chromosomes < MAX_CHROMOSOMES) {
    strncpy(chromosomes[*n_chromosomes].name, chr_name, MAX_CHR_LEN - 1);
    chromosomes[*n_chromosomes].name[MAX_CHR_LEN - 1] = '\0';
    chromosomes[*n_chromosomes].min_pos = min_pos;
    chromosomes[*n_chromosomes].max_pos = max_pos;
    chromosomes[*n_chromosomes].has_data = 1;
    (*n_chromosomes)++;
  } else {
    if (ref_fai) {  // Only print warning if we're in verbose mode
      Rprintf("Warning: Maximum chromosomes reached (%d), ignoring %s\n", 
              MAX_CHROMOSOMES, chr_name);
    }
  }
}

int generate_bins(config_t *config) {
  config->bins = (genomic_bin_t*)R_alloc(MAX_BINS, sizeof(genomic_bin_t));
  config->n_bins = 0;
  
  for (int chr_idx = 0; chr_idx < config->n_chromosomes; chr_idx++) {
    chromosome_t *chr = &config->chromosomes[chr_idx];
    
    // Start from the minimum position with data, rounded down to bin boundary
    int32_t start_pos = (chr->min_pos / config->bin_size) * config->bin_size;
    // End at the maximum position with data, rounded up to bin boundary  
    int32_t end_pos = ((chr->max_pos / config->bin_size) + 1) * config->bin_size;
    
    for (int32_t pos = start_pos; pos < end_pos; pos += config->bin_size) {
      if (config->n_bins >= MAX_BINS) {
        Rprintf("Warning: Maximum bins reached (%d), truncating\n", MAX_BINS);
        goto done;
      }
      
      genomic_bin_t *bin = &config->bins[config->n_bins];
      strncpy(bin->chr, chr->name, MAX_CHR_LEN - 1);
      bin->chr[MAX_CHR_LEN - 1] = '\0';
      bin->start = pos;
      bin->end = pos + config->bin_size;
      bin->bin_id = config->n_bins;
      
      config->n_bins++;
    }
  }
  
  done:
    if (config->verbose) {
      Rprintf("Generated %d bins of size %d covering data regions\n",
              config->n_bins, config->bin_size);
      for (int i = 0; i < config->n_bins; i++) {
        Rprintf("  Bin %d: %s:%d-%d\n", i, config->bins[i].chr,
                config->bins[i].start, config->bins[i].end);
      }
    }

    return 0;
}

// File stream management
int init_file_stream_for_bin(file_stream_t *stream, const char *filename, int file_idx,
                             genomic_bin_t *bin, config_t *config) {
  
  stream->file_idx = file_idx;
  stream->filename = (char*)R_alloc(strlen(filename) + 1, sizeof(char));
  strcpy(stream->filename, filename);
  stream->has_current_record = 0;
  stream->has_output_data = 0;
  stream->line_buffer = (kstring_t){0, 0, 0};
  
  stream->hts_file = hts_open(filename, "r");
  if (!stream->hts_file) {
    return 1;
  }
  
  stream->tbx_handle = tbx_index_load(filename);
  if (!stream->tbx_handle) {
    return 1;
  }
  
  // Create region string for this bin
  char *region = format_region_string(bin);
  // DEBUG: Region querying debug commented out
  // if (config->verbose) {
  //   Rprintf("    File %d: Querying region %s\n", file_idx, region);
  // }
  stream->iterator = tbx_itr_querys(stream->tbx_handle, region);
  free(region);
  
  if (!stream->iterator) {
    return 0;  // No data in this region
  }
  
  // Read first record
  stream->has_current_record = advance_stream(stream, config);

  // DEBUG: First record debug commented out
  // if (config->verbose && stream->has_current_record) {
  //   Rprintf("    File %d: First record at %s:%d:%c\n", file_idx,
  //           stream->current_record.chr, stream->current_record.start, stream->current_record.strand);
  // }

  return 0;
}

char* format_region_string(genomic_bin_t *bin) {
  char *region = (char*)malloc(256);
  snprintf(region, 256, "%s:%d-%d", bin->chr, bin->start, bin->end);
  return region;
}

int advance_stream(file_stream_t *stream, config_t *config) {
  if (!stream->iterator) {
    return 0;
  }
  
  while (tbx_itr_next(stream->hts_file, stream->tbx_handle, 
                      stream->iterator, &stream->line_buffer) >= 0) {
    
    modkit_record_t record;
    if (parse_modkit_record(stream->line_buffer.s, &record) == 0) {
      if (!apply_record_filters(&record, config)) {
        continue;
      }
      
      apply_coordinate_transforms(&record, config);
      stream->current_record = record;
      return 1;
    }
  }
  
  return 0;
}

void cleanup_streams(file_stream_t *streams, int n_files) {
  for (int i = 0; i < n_files; i++) {
    if (streams[i].line_buffer.s) {
      free(streams[i].line_buffer.s);
    }
    if (streams[i].iterator) {
      tbx_itr_destroy(streams[i].iterator);
    }
    if (streams[i].tbx_handle) {
      tbx_destroy(streams[i].tbx_handle);
    }
    if (streams[i].hts_file) {
      hts_close(streams[i].hts_file);
    }
  }
}

// Coordinate processing functions
coord_t find_minimum_coordinate(file_stream_t *streams, int n_files) {
  coord_t min_coord = {{0}, 0, 0};
  int found_first = 0;

  for (int i = 0; i < n_files; i++) {
    if (streams[i].has_current_record) {
      coord_t coord;
      strncpy(coord.chr, streams[i].current_record.chr, MAX_CHR_LEN - 1);
      coord.chr[MAX_CHR_LEN - 1] = '\0';
      coord.start = streams[i].current_record.start;
      coord.strand = streams[i].current_record.strand;

      if (!found_first || coordinate_compare(&coord, &min_coord) < 0) {
        min_coord = coord;
        found_first = 1;
      }
    }
  }

  return min_coord;
}

// NEW: Find minimum coordinate with strand normalization for merging
coord_t find_minimum_coordinate_with_merging(file_stream_t *streams, int n_files, config_t *config) {
  coord_t min_coord = {{0}, 0, 0};
  int found_first = 0;

  for (int i = 0; i < n_files; i++) {
    if (streams[i].has_current_record) {
      coord_t coord = normalize_cpg_coordinate(&streams[i].current_record, config);

      if (!found_first || coordinate_compare(&coord, &min_coord) < 0) {
        min_coord = coord;
        found_first = 1;
      }
    }
  }

  return min_coord;
}

int coordinate_compare(const coord_t *a, const coord_t *b) {
  int chr_cmp = strcmp(a->chr, b->chr);
  if (chr_cmp != 0) return chr_cmp;
  
  if (a->start != b->start) {
    return (a->start < b->start) ? -1 : 1;
  }
  
  if (a->strand != b->strand) {
    return (a->strand < b->strand) ? -1 : 1;
  }
  
  return 0;
}

int coordinates_equal(const coord_t *a, const coord_t *b) {
  return (strcmp(a->chr, b->chr) == 0 && 
          a->start == b->start && 
          a->strand == b->strand);
}

int any_stream_has_data(file_stream_t *streams, int n_files) {
  for (int i = 0; i < n_files; i++) {
    if (streams[i].has_current_record) {
      return 1;
    }
  }
  return 0;
}

void collect_coordinate_data(file_stream_t *streams, int n_files, const coord_t *target_coord, config_t *config) {
  for (int i = 0; i < n_files; i++) {
    if (streams[i].has_current_record) {
      coord_t current_coord;
      strncpy(current_coord.chr, streams[i].current_record.chr, MAX_CHR_LEN - 1);
      current_coord.chr[MAX_CHR_LEN - 1] = '\0';
      current_coord.start = streams[i].current_record.start;
      current_coord.strand = streams[i].current_record.strand;
      
      if (coordinates_equal(&current_coord, target_coord)) {
        streams[i].has_output_data = 1;
        streams[i].output_beta = streams[i].current_record.frac_mod;
        streams[i].output_cov = streams[i].current_record.nvalid_cov;
        streams[i].has_current_record = advance_stream(&streams[i], config);
      }
    }
  }
}

void reset_output_data(file_stream_t *streams, int n_files) {
  for (int i = 0; i < n_files; i++) {
    streams[i].has_output_data = 0;
    streams[i].output_beta = R_NaN;
    streams[i].output_cov = 0;
  }
}

// Record processing functions
int parse_modkit_record(char *line, modkit_record_t *record) {
  char *fields[20];
  int n_fields = 0;
  
  char *token = strtok(line, "\t\n");
  while (token && n_fields < 20) {
    fields[n_fields++] = token;
    token = strtok(NULL, "\t\n");
  }
  
  if (n_fields < 18) {
    return -1;
  }
  
  strncpy(record->chr, fields[0], MAX_CHR_LEN - 1);
  record->chr[MAX_CHR_LEN - 1] = '\0';
  record->start = atoi(fields[1]);
  record->end = atoi(fields[2]);
  record->mod_code = fields[3][0];
  record->strand = fields[5][0];
  record->nvalid_cov = atoi(fields[9]);
  record->frac_mod = atof(fields[10]);
  record->nmod = atoi(fields[11]);
  record->ncanonical = atoi(fields[12]);
  record->nother_mod = atoi(fields[13]);
  record->nfail = atoi(fields[15]);
  record->nnocall = atoi(fields[17]);
  
  return 0;
}

int apply_record_filters(modkit_record_t *record, config_t *config) {
  if (record->mod_code != config->target_mod) {
    return 0;
  }
  
  if (record->nvalid_cov < config->min_coverage) {
    return 0;
  }
  
  if (config->quality_filter > 0) {
    int total_calls = record->nmod + record->ncanonical + 
      record->nother_mod + record->nfail;
    if (total_calls > 0) {
      double fail_rate = (double)record->nfail / total_calls;
      if (fail_rate > config->quality_filter) {
        return 0;
      }
    }
  }
  
  return 1;
}

void apply_coordinate_transforms(modkit_record_t *record, config_t *config) {
  // NO coordinate conversion - keep original 0-based coordinates
  // NOTE: Strand merging is now handled in the new aggregation functions
  // This function is kept for backwards compatibility but does nothing
  // when combine_strands is enabled, since coordinate normalization
  // happens in normalize_cpg_coordinate() instead
}

// Reference FASTA handling functions
int init_reference(config_t *config) {
  if (!config->ref_fasta) {
    return 0;
  }
  
  ref_fai = fai_load(config->ref_fasta);
  if (!ref_fai) {
    error("Cannot load FASTA index for %s. Make sure %s.fai exists or run: samtools faidx %s", 
          config->ref_fasta, config->ref_fasta, config->ref_fasta);
  }
  
  // Initialize cache
  ref_cache.valid = 0;
  ref_cache.sequence = NULL;
  
  return 0;
}

void cleanup_reference() {
  if (ref_fai) {
    fai_destroy(ref_fai);
    ref_fai = NULL;
  }
  
  if (ref_cache.sequence) {
    free(ref_cache.sequence);
    ref_cache.sequence = NULL;
  }
  ref_cache.valid = 0;
}

char get_ref_base_cached(const char *chr, int32_t pos) {
  if (!ref_fai) {
    return 'N';
  }
  
  // Check if position is in cache
  if (ref_cache.valid && 
      strcmp(ref_cache.chr, chr) == 0 &&
      pos >= ref_cache.start && 
      pos < ref_cache.end) {
    return toupper(ref_cache.sequence[pos - ref_cache.start]);
  }
  
  // Update cache with new chunk
  int cache_start = (pos / CACHE_SIZE) * CACHE_SIZE;
  int cache_end = cache_start + CACHE_SIZE;
  
  if (ref_cache.sequence) {
    free(ref_cache.sequence);
  }
  
  int len;
  ref_cache.sequence = faidx_fetch_seq(ref_fai, chr, cache_start, cache_end - 1, &len);
  if (!ref_cache.sequence) {
    ref_cache.valid = 0;
    return 'N';
  }
  
  strncpy(ref_cache.chr, chr, MAX_CHR_LEN - 1);
  ref_cache.chr[MAX_CHR_LEN - 1] = '\0';
  ref_cache.start = cache_start;
  ref_cache.end = cache_start + len;
  ref_cache.valid = 1;
  
  // Return the requested base
  if (pos >= ref_cache.start && pos < ref_cache.end) {
    return toupper(ref_cache.sequence[pos - ref_cache.start]);
  }
  
  return 'N';
}

char* get_context(const char *chr, int32_t pos, char strand) {
  if (!ref_fai) {
    char *result = (char*)malloc(4);
    strcpy(result, "NNN");
    return result;
  }
  
  // Get 3-base context (position + 2 downstream bases)
  int len;
  char *seq = faidx_fetch_seq(ref_fai, chr, pos, pos + 2, &len);
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

char complement_base(char base) {
  switch (toupper(base)) {
  case 'A': return 'T';
  case 'T': return 'A';
  case 'C': return 'G';
  case 'G': return 'C';
  default: return 'N';
  }
}

void reverse_complement_inplace(char *seq, int len) {
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

// NEW: Fixed binning implementation
int extract_chromosomes_from_bedmethyl_files(config_t *config) {
  if (!config->input_files || config->n_files == 0) {
    Rprintf("No input files provided\n");
    return 1;
  }

  // Use a simple array to collect unique chromosome names
  char unique_chromosomes[MAX_CHROMOSOMES][MAX_CHR_LEN];
  int n_unique_chromosomes = 0;

  Rprintf("Discovering chromosomes from %d bedMethyl files...\n", config->n_files);

  // For each input file, use tabix to list chromosomes
  for (int file_idx = 0; file_idx < config->n_files; file_idx++) {
    htsFile *fp = hts_open(config->input_files[file_idx], "r");
    if (!fp) {
      Rprintf("Warning: Cannot open file %s\n", config->input_files[file_idx]);
      continue;
    }

    tbx_t *tbx = tbx_index_load(config->input_files[file_idx]);
    if (!tbx) {
      Rprintf("Warning: Cannot load tabix index for %s\n", config->input_files[file_idx]);
      hts_close(fp);
      continue;
    }

    // Get sequence names from tabix index
    int n_seqs;
    const char **seq_names = tbx_seqnames(tbx, &n_seqs);

    for (int i = 0; i < n_seqs; i++) {
      // Check if this chromosome is already in our list
      int already_exists = 0;
      for (int j = 0; j < n_unique_chromosomes; j++) {
        if (strcmp(unique_chromosomes[j], seq_names[i]) == 0) {
          already_exists = 1;
          break;
        }
      }

      // Add new chromosome if not already present
      if (!already_exists && n_unique_chromosomes < MAX_CHROMOSOMES) {
        strncpy(unique_chromosomes[n_unique_chromosomes], seq_names[i], MAX_CHR_LEN - 1);
        unique_chromosomes[n_unique_chromosomes][MAX_CHR_LEN - 1] = '\0';
        n_unique_chromosomes++;
      }
    }

    free(seq_names);
    tbx_destroy(tbx);
    hts_close(fp);
  }

  if (n_unique_chromosomes == 0) {
    Rprintf("Error: No chromosomes found in input files\n");
    return 1;
  }

  // Store discovered chromosomes in config for later filtering
  config->chromosomes = (chromosome_t*)R_alloc(n_unique_chromosomes, sizeof(chromosome_t));
  config->n_chromosomes = n_unique_chromosomes;

  for (int i = 0; i < n_unique_chromosomes; i++) {
    strncpy(config->chromosomes[i].name, unique_chromosomes[i], MAX_CHR_LEN - 1);
    config->chromosomes[i].name[MAX_CHR_LEN - 1] = '\0';
    config->chromosomes[i].has_data = 1; // We know these have data
    config->chromosomes[i].min_pos = 0;
    config->chromosomes[i].max_pos = 0;
  }

  Rprintf("Found %d unique chromosomes: ", n_unique_chromosomes);
  for (int i = 0; i < n_unique_chromosomes; i++) {
    Rprintf("%s%s", config->chromosomes[i].name, (i < n_unique_chromosomes - 1) ? ", " : "\n");
  }

  return 0;
}

int extract_chromosome_lengths_from_fai(config_t *config) {
  if (!config->ref_fasta) {
    Rprintf("Error: No reference FASTA provided for fixed binning\n");
    return 1;
  }

  if (!config->chromosomes || config->n_chromosomes == 0) {
    Rprintf("Error: No chromosomes discovered from bedMethyl files\n");
    return 1;
  }

  char fai_path[2048];
  snprintf(fai_path, sizeof(fai_path), "%s.fai", config->ref_fasta);

  FILE *fai_file = fopen(fai_path, "r");
  if (!fai_file) {
    Rprintf("Error: Cannot open FASTA index file: %s\n", fai_path);
    return 1;
  }

  // Allocate memory for chromosome lengths (only for discovered chromosomes)
  config->chr_lengths = (chromosome_length_t*)R_alloc(config->n_chromosomes, sizeof(chromosome_length_t));
  config->n_chr_lengths = 0;

  Rprintf("Extracting chromosome lengths for %d discovered chromosomes...\n", config->n_chromosomes);

  char line[1024];
  int chromosomes_found = 0;
  int chromosomes_missing = 0;

  // Parse .fai file: chr_name \t length \t offset \t line_bases \t line_width
  while (fgets(line, sizeof(line), fai_file)) {
    char chr_name[MAX_CHR_LEN];
    int64_t chr_length;

    if (sscanf(line, "%63s\t%ld", chr_name, &chr_length) == 2) {
      // Check if this chromosome was found in bedMethyl files
      for (int i = 0; i < config->n_chromosomes; i++) {
        if (strcmp(config->chromosomes[i].name, chr_name) == 0) {
          // Found a match - store the length
          strncpy(config->chr_lengths[config->n_chr_lengths].name, chr_name, MAX_CHR_LEN - 1);
          config->chr_lengths[config->n_chr_lengths].name[MAX_CHR_LEN - 1] = '\0';
          config->chr_lengths[config->n_chr_lengths].length = chr_length;
          config->n_chr_lengths++;
          chromosomes_found++;
          break;
        }
      }
    }
  }

  fclose(fai_file);

  // Check for missing chromosomes
  chromosomes_missing = config->n_chromosomes - chromosomes_found;

  if (chromosomes_missing > 0) {
    Rprintf("Warning: %d chromosomes from bedMethyl files not found in FASTA index:\n", chromosomes_missing);
    for (int i = 0; i < config->n_chromosomes; i++) {
      int found = 0;
      for (int j = 0; j < config->n_chr_lengths; j++) {
        if (strcmp(config->chromosomes[i].name, config->chr_lengths[j].name) == 0) {
          found = 1;
          break;
        }
      }
      if (!found) {
        Rprintf("  %s (skipping)\n", config->chromosomes[i].name);
      }
    }
  }

  Rprintf("Found lengths for %d/%d chromosomes in FASTA index\n", chromosomes_found, config->n_chromosomes);

  return 0;
}

int generate_fixed_genomic_bins(config_t *config) {
  config->bins = (genomic_bin_t*)R_alloc(MAX_BINS, sizeof(genomic_bin_t));
  config->n_bins = 0;

  if (config->verbose) {
    Rprintf("Generating fixed genomic bins (bin_size=%d)\n", config->bin_size);
  }

  for (int chr_idx = 0; chr_idx < config->n_chr_lengths; chr_idx++) {
    chromosome_length_t *chr = &config->chr_lengths[chr_idx];

    // Generate fixed bins across entire chromosome
    for (int64_t pos = 1; pos <= chr->length; pos += config->bin_size) {
      if (config->n_bins >= MAX_BINS) {
        Rprintf("Warning: Maximum bins reached (%d), truncating\n", MAX_BINS);
        goto done;
      }

      genomic_bin_t *bin = &config->bins[config->n_bins];
      strncpy(bin->chr, chr->name, MAX_CHR_LEN - 1);
      bin->chr[MAX_CHR_LEN - 1] = '\0';
      bin->start = (int32_t)pos;
      bin->end = (int32_t)(pos + config->bin_size - 1 < chr->length ?
                           pos + config->bin_size - 1 : chr->length);
      bin->bin_id = config->n_bins;

      config->n_bins++;
    }

    if (config->verbose) {
      int n_bins_for_chr = (chr->length + config->bin_size - 1) / config->bin_size;
      Rprintf("  %s: %d bins (length=%ld)\n", chr->name, n_bins_for_chr, chr->length);
    }
  }

  done:
    if (config->verbose) {
      Rprintf("Generated %d fixed genomic bins\n", config->n_bins);
    }

    return 0;
}

// Coordinate set management functions
void init_coordinate_set(coordinate_set_t *coord_set, int initial_capacity) {
  coord_set->coords = (coord_t*)R_alloc(initial_capacity, sizeof(coord_t));
  coord_set->n_coords = 0;
  coord_set->capacity = initial_capacity;
}

void add_coordinate_to_set(coordinate_set_t *coord_set, const coord_t *coord) {
  // Check if coordinate already exists (avoid duplicates)
  for (int i = 0; i < coord_set->n_coords; i++) {
    if (coordinates_equal(&coord_set->coords[i], coord)) {
      return; // Already exists
    }
  }

  // Expand capacity if needed
  if (coord_set->n_coords >= coord_set->capacity) {
    int new_capacity = coord_set->capacity * 2;
    coord_t *new_coords = (coord_t*)R_alloc(new_capacity, sizeof(coord_t));
    memcpy(new_coords, coord_set->coords, coord_set->n_coords * sizeof(coord_t));
    coord_set->coords = new_coords;
    coord_set->capacity = new_capacity;
  }

  // Add new coordinate
  coord_set->coords[coord_set->n_coords] = *coord;
  coord_set->n_coords++;
}

int coord_compare_for_sort(const void *a, const void *b) {
  const coord_t *coord_a = (const coord_t*)a;
  const coord_t *coord_b = (const coord_t*)b;
  return coordinate_compare(coord_a, coord_b);
}

void sort_coordinate_set(coordinate_set_t *coord_set) {
  qsort(coord_set->coords, coord_set->n_coords, sizeof(coord_t), coord_compare_for_sort);
}

void cleanup_coordinate_set(coordinate_set_t *coord_set) {
  // R_alloc'd memory is automatically cleaned up by R
  coord_set->coords = NULL;
  coord_set->n_coords = 0;
  coord_set->capacity = 0;
}

coordinate_set_t collect_all_coordinates_in_bin(char **input_files, int n_files, genomic_bin_t *bin, config_t *config) {
  coordinate_set_t coord_set;
  init_coordinate_set(&coord_set, 10000); // Start with 10K coordinates

  if (config->verbose) {
    Rprintf("Collecting all coordinates in bin %s:%d-%d\n", bin->chr, bin->start, bin->end);
  }

  // Collect coordinates from each file
  for (int file_idx = 0; file_idx < n_files; file_idx++) {
    file_stream_t stream;
    memset(&stream, 0, sizeof(file_stream_t));

    // Initialize stream for this bin
    if (init_file_stream_for_bin(&stream, input_files[file_idx], file_idx, bin, config) == 0) {
      int coord_count = 0;

      // Read all records in this bin
      while (stream.has_current_record) {
        coord_t coord;

        if (config->combine_strands &&
            (config->target_mod == 'm' || config->target_mod == 'h' ||
             config->target_mod == 'c' || config->target_mod == 'g')) {
          // Use normalized coordinates for strand merging
          coord = normalize_cpg_coordinate(&stream.current_record, config);
        } else {
          // Use original coordinates
          strncpy(coord.chr, stream.current_record.chr, MAX_CHR_LEN - 1);
          coord.chr[MAX_CHR_LEN - 1] = '\0';
          coord.start = stream.current_record.start;
          coord.strand = stream.current_record.strand;
        }

        // Add to coordinate set (duplicates will be ignored)
        add_coordinate_to_set(&coord_set, &coord);
        coord_count++;

        // Advance to next record
        stream.has_current_record = advance_stream(&stream, config);
      }

      // DEBUG: Coordinate counting debug commented out
      // if (config->verbose) {
      //   Rprintf("  File %d: %d coordinates\n", file_idx, coord_count);
      // }

      // Cleanup stream
      cleanup_streams(&stream, 1);
    }
  }

  // Sort coordinates for efficient processing
  sort_coordinate_set(&coord_set);

  if (config->verbose) {
    Rprintf("  Total unique coordinates: %d\n", coord_set.n_coords);
  }

  return coord_set;
}


int query_file_for_coordinate(const char *filename, const coord_t *target_coord,
                             file_stream_t *result_stream, config_t *config) {
  file_stream_t temp_stream;
  memset(&temp_stream, 0, sizeof(file_stream_t));

  // Create a small bin around the target coordinate for efficient querying
  genomic_bin_t query_bin;
  strncpy(query_bin.chr, target_coord->chr, MAX_CHR_LEN - 1);
  query_bin.chr[MAX_CHR_LEN - 1] = '\0';
  query_bin.start = target_coord->start;
  query_bin.end = target_coord->start + 1;

  // Initialize stream for this specific region
  if (init_file_stream_for_bin(&temp_stream, filename, 0, &query_bin, config) != 0) {
    return 0; // Failed to initialize
  }

  // Search for the exact coordinate
  while (temp_stream.has_current_record) {
    coord_t current_coord;

    if (config->combine_strands &&
        (config->target_mod == 'm' || config->target_mod == 'h' ||
         config->target_mod == 'c' || config->target_mod == 'g')) {
      current_coord = normalize_cpg_coordinate(&temp_stream.current_record, config);
    } else {
      strncpy(current_coord.chr, temp_stream.current_record.chr, MAX_CHR_LEN - 1);
      current_coord.chr[MAX_CHR_LEN - 1] = '\0';
      current_coord.start = temp_stream.current_record.start;
      current_coord.strand = temp_stream.current_record.strand;
    }

    if (coordinates_equal(&current_coord, target_coord)) {
      // Found matching coordinate - collect data
      if (config->combine_strands) {
        // Use strand aggregation logic
        strand_aggregate_t agg;
        init_strand_aggregate(&agg);
        add_record_to_aggregate(&agg, &temp_stream.current_record);

        // Check for additional records at same coordinate (strand pairs)
        temp_stream.has_current_record = advance_stream(&temp_stream, config);
        while (temp_stream.has_current_record) {
          coord_t next_coord = normalize_cpg_coordinate(&temp_stream.current_record, config);
          if (coordinates_equal(&next_coord, target_coord)) {
            add_record_to_aggregate(&agg, &temp_stream.current_record);
            temp_stream.has_current_record = advance_stream(&temp_stream, config);
          } else {
            break;
          }
        }

        finalize_strand_aggregate(&agg);
        result_stream->output_beta = agg.combined_beta;
        result_stream->output_cov = agg.total_coverage;
      } else {
        // Use direct values
        result_stream->output_beta = temp_stream.current_record.frac_mod;
        result_stream->output_cov = temp_stream.current_record.nvalid_cov;
      }

      cleanup_streams(&temp_stream, 1);
      return 1; // Found data
    }

    // Move to next record
    temp_stream.has_current_record = advance_stream(&temp_stream, config);
  }

  cleanup_streams(&temp_stream, 1);
  return 0; // No data found
}

// NEW: Strand merging implementation
coord_t normalize_cpg_coordinate(const modkit_record_t *record, config_t *config) {
  coord_t normalized_coord;
  strncpy(normalized_coord.chr, record->chr, MAX_CHR_LEN - 1);
  normalized_coord.chr[MAX_CHR_LEN - 1] = '\0';

  // For CpG strand merging:
  // Plus strand CpG at position N stays at N
  // Minus strand CpG at position N+1 maps to N (the corresponding C position)
  if (config->combine_strands &&
      (config->target_mod == 'm' || config->target_mod == 'h' ||
       config->target_mod == 'c' || config->target_mod == 'g')) {

    if (record->strand == '-') {
      // Minus strand G at pos N+1 represents CpG at pos N
      normalized_coord.start = record->start - 1;
    } else {
      // Plus strand C at pos N
      normalized_coord.start = record->start;
    }
    normalized_coord.strand = '+';  // All normalized to plus strand
  } else {
    // No strand merging - keep original coordinates
    normalized_coord.start = record->start;
    normalized_coord.strand = record->strand;
  }

  return normalized_coord;
}

int cpg_coordinates_equal(const coord_t *a, const coord_t *b) {
  return (strcmp(a->chr, b->chr) == 0 &&
          a->start == b->start &&
          a->strand == b->strand);
}

void init_strand_aggregate(strand_aggregate_t *agg) {
  agg->total_meth_calls = 0;
  agg->total_canonical_calls = 0;
  agg->total_coverage = 0;
  agg->combined_beta = 0.0;
  agg->has_plus_strand = 0;
  agg->has_minus_strand = 0;
  agg->record_count = 0;
}

void add_record_to_aggregate(strand_aggregate_t *agg, const modkit_record_t *record) {
  agg->total_meth_calls += record->nmod;
  agg->total_canonical_calls += record->ncanonical;
  agg->total_coverage += record->nvalid_cov;

  // Track which strands we have data for
  if (record->strand == '+') {
    agg->has_plus_strand = 1;
  } else if (record->strand == '-') {
    agg->has_minus_strand = 1;
  }

  agg->record_count++;
}

void finalize_strand_aggregate(strand_aggregate_t *agg) {
  // Calculate combined beta using the proper formula:
  // β_combined = (n_meth_plus + n_meth_minus) / (n_total_plus + n_total_minus)
  int32_t total_calls = agg->total_meth_calls + agg->total_canonical_calls;

  if (total_calls > 0) {
    agg->combined_beta = (double)agg->total_meth_calls / (double)total_calls;
  } else {
    agg->combined_beta = 0.0;
  }
}

void collect_coordinate_data_with_merging(file_stream_t *streams, int n_files, const coord_t *target_coord, config_t *config) {
  // Array to store aggregates for each file
  strand_aggregate_t *file_aggregates = (strand_aggregate_t*)R_alloc(n_files, sizeof(strand_aggregate_t));

  // Initialize aggregates
  for (int i = 0; i < n_files; i++) {
    init_strand_aggregate(&file_aggregates[i]);
    streams[i].has_output_data = 0;
  }

  // DEBUG: Debug output commented out
  // // DEBUG: if (config->verbose && target_coord->start == 10003) {
  //   Rprintf("DEBUG: Processing coordinate %s:%d:%c\n",
  //           target_coord->chr, target_coord->start, target_coord->strand);
  // }

  // Collect all records that match the target coordinate (potentially from multiple strands)
  for (int i = 0; i < n_files; i++) {
    // Check current record
    if (streams[i].has_current_record) {
      coord_t current_normalized = normalize_cpg_coordinate(&streams[i].current_record, config);

      // DEBUG: Position comparison debug commented out
      // if (config->verbose && target_coord->start == 10003) {
      //   Rprintf("    File %d: current_pos=%d:%c, normalized=%d:%c, match=%d\n",
      //           i, streams[i].current_record.start, streams[i].current_record.strand,
      //           current_normalized.start, current_normalized.strand,
      //           cpg_coordinates_equal(&current_normalized, target_coord));
      // }

      // Skip records that are before the target coordinate
      while (streams[i].has_current_record) {
        current_normalized = normalize_cpg_coordinate(&streams[i].current_record, config);
        int coord_cmp = coordinate_compare(&current_normalized, target_coord);

        if (coord_cmp < 0) {
          // Current record is before target, advance to next record
          // DEBUG: Record skipping debug commented out
          // if (config->verbose && target_coord->start == 10003) {
          //   Rprintf("      File %d: Skipping record at %d (before target %d)\n",
          //           i, current_normalized.start, target_coord->start);
          // }
          streams[i].has_current_record = advance_stream(&streams[i], config);
        } else {
          break; // Found record at or after target coordinate
        }
      }

      // Now check if we have a record that matches the target
      if (streams[i].has_current_record) {
        current_normalized = normalize_cpg_coordinate(&streams[i].current_record, config);

        if (cpg_coordinates_equal(&current_normalized, target_coord)) {
          add_record_to_aggregate(&file_aggregates[i], &streams[i].current_record);
          // DEBUG: Record addition debug commented out
          // if (config->verbose && target_coord->start == 10003) {
          //   Rprintf("      File %d: Added record cov=%d beta=%.3f\n",
          //           i, streams[i].current_record.nvalid_cov, streams[i].current_record.frac_mod);
          // }
          streams[i].has_current_record = advance_stream(&streams[i], config);

          // Check if there's another record at the same coordinate (strand pair)
          while (streams[i].has_current_record) {
            coord_t next_normalized = normalize_cpg_coordinate(&streams[i].current_record, config);
            if (cpg_coordinates_equal(&next_normalized, target_coord)) {
              add_record_to_aggregate(&file_aggregates[i], &streams[i].current_record);
              // DEBUG: Additional record debug commented out
              // if (config->verbose && target_coord->start == 10003) {
              //   Rprintf("      File %d: Added additional record cov=%d beta=%.3f\n",
              //           i, streams[i].current_record.nvalid_cov, streams[i].current_record.frac_mod);
              // }
              streams[i].has_current_record = advance_stream(&streams[i], config);
            } else {
              break;  // No more records at this coordinate
            }
          }
        }
      }
    }
  }

  // Finalize aggregates and set output data
  for (int i = 0; i < n_files; i++) {
    if (file_aggregates[i].record_count > 0) {
      finalize_strand_aggregate(&file_aggregates[i]);
      streams[i].has_output_data = 1;
      streams[i].output_beta = file_aggregates[i].combined_beta;
      streams[i].output_cov = file_aggregates[i].total_coverage;

      // DEBUG: Final result debug commented out
      // if (config->verbose && target_coord->start == 10003) {
      //   Rprintf("    File %d final: cov=%d beta=%.3f (from %d records)\n",
      //           i, streams[i].output_cov, streams[i].output_beta, file_aggregates[i].record_count);
      // }
    } else {
      // DEBUG: No data debug commented out
      // if (config->verbose && target_coord->start == 10003) {
      //   Rprintf("    File %d: No data\n", i);
      // }
    }
  }
}

// Cleanup function
void cleanup_config(config_t *config) {
  // R handles memory cleanup automatically for R_alloc'd memory
  // This function is mainly for completeness and any malloc'd memory
  if (config->region) {
    // Only free if it was malloc'd, not R_alloc'd
    // In R interface, this will be R_alloc'd, so no free needed
  }
  if (config->output_file) {
    // Only free if it was malloc'd, not R_alloc'd
    // In R interface, this will be R_alloc'd, so no free needed
  }
  if (config->ref_fasta) {
    // Only free if it was malloc'd, not R_alloc'd
    // In R interface, this will be R_alloc'd, so no free needed
  }
}