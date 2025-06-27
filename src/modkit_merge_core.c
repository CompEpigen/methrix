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
  stream->iterator = tbx_itr_querys(stream->tbx_handle, region);
  free(region);
  
  if (!stream->iterator) {
    return 0;  // No data in this region
  }
  
  // Read first record
  stream->has_current_record = advance_stream(stream, config);
  
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
  
  if (config->combine_strands && 
      (config->target_mod == 'm' || config->target_mod == 'h' || 
      config->target_mod == 'c' || config->target_mod == 'g')) {
    if (record->strand == '-') {
      record->start -= 1;
    }
    record->strand = '+';
  }
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