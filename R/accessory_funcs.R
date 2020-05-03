# Tiny function to check if object is h5
is_h5 = function(m) {
  return(m@metadata$is_h5)
}


get_source_idx = function(protocol = NULL) {
  if (protocol == "MethylcTools") {
    return(list(col_idx = c(chr = 1, start = 2, strand = 3, context = 4,
                            cnp_rate = 5, M = 6, U = 7),
                col_classes = c("character", "numeric",
                                "character", "character", "numeric", "integer",
                                "integer"),
                fix_missing = c("end := start+1", "cov := M+U",
                                "beta := M/cov")))
  } else {
    # Bismark and methyldackel have same output format
    return(list(col_idx = c(chr = 1, start = 2, end = 3, beta = 4,
                            M = 5, U = 6),
                col_classes = c("character", "numeric", "numeric",
                                "numeric", "numeric", "numeric"),
                fix_missing = c("cov := M+U",
                                "strand := '.'")))
  }
}

#--------------------------------------------------------------------------------------------------------------------------

# Parse custom indices and return missing info
parse_source_idx = function(chr = NULL, start = NULL, end = NULL, strand = NULL,
                            beta = NULL, n_meth = NULL, n_unmeth = NULL,
                            cov = NULL, beta_fract = FALSE,
                            verbose = TRUE) {
  
  # mandatory chr and start field
  if (is.null(chr) | is.null(start)) {
    stop("missing chromosome/start indices\nUse pipeline argument if the files are from Bismark, MethyDeckal, or MethylcTools",
         call. = FALSE)
  }
  
  # See if any indices are duplicated
  if (length(which(duplicated(c(chr, start, end, strand, beta, n_meth,
                                n_unmeth, cov)))) > 0) {
    stop("Duplicated indices.", call. = FALSE)
  }
  
  # Check maximum betavalues (Can be 1 or 100)
  fix_missing = vector()
  if (is.null(strand)) {
    fix_missing = "strand := '*'"
  }
  
  # If beta and cov are missing
  if (all(is.null(beta), is.null(cov))) {
    if (is.null(n_meth) | is.null(n_unmeth)) {
      stop("Missing beta or coverage values.\nU and M are not available either!",
           call. = FALSE)
    } else {
      if (verbose) {
        message("--Missing beta and coverage info. Estimating them from M and U values")
      }
      
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, M = n_meth,
                              U = n_unmeth,
                              cov = cov),
                  fix_missing = c(fix_missing, "cov := M+U",
                                  "beta := M/(M+U)")))
    }
  } else if (is.null(beta) & !is.null(cov)) {
    # If anyone of them is present (case-1: coverage available, estimate
    # beta)
    if (all(is.null(n_meth), is.null(n_unmeth))) {
      stop("Missing beta values but coverage info available.\nEither U or M are required for estimating beta values!",
           call. = FALSE)
    } else if (all(!is.null(n_meth), !is.null(n_unmeth))) {
      message("--Estimating beta values from M and U")
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, M = n_meth,
                              U = n_unmeth,  cov = cov),
                  fix_missing = c(fix_missing, "beta := M/(M+U")))
    } else if (!is.null(n_meth)) {
      # M available
      message("--Estimating beta values from M and coverage")
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, M = n_meth,
                              U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "beta := M/cov",
                                  "U := cov - M")))
    } else if (!is.null(n_unmeth)) {
      # U available
      message("--Estimating beta values from U and coverage")
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta,M = n_meth,
                              U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, paste0("beta := 1- (U/cov)"),
                                  "M := cov - U")))
    }
  } else if (!is.null(beta) & is.null(cov)) {
    # If anyone of them is present (case-2: beta available, estimate
    # coverage)
    if (all(is.null(n_meth), is.null(n_unmeth))) {
      stop("Missing coverage info but beta values are available.\nU and M are required for estimating coverage values!",
           call. = FALSE)
    } else {
      if (verbose) {
        message("--Estimating coverage from M and U")
      }
      
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, M = n_meth,
                              U = n_unmeth, cov = cov),
                  fix_missing = c(fix_missing, "cov := M+U")))
    }
  } else if (!is.null(beta) & !is.null(cov)) {
    # If both present (case-3: beta and coverage available, but missing M
    # and U)
    if (all(is.null(n_meth), is.null(n_unmeth))) {
      if (verbose) {
        message("--Estimating M and U from coverage and beta values")
      }
      
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, cov = cov),
                  fix_missing = c("M := as.integer(cov * beta)",
                                  "U := cov - M")))
    } else {
      if (verbose) {
        message("--All fields are present. Nice.")
      }
      
      return(list(col_idx = c(chr = chr, start = start, end = end,
                              strand = strand, beta = beta, cov = cov, M = n_meth,
                              U = n_unmeth),
                  fix_missing = NULL))
    }
  }
}



#--------------------------------------------------------------------------------------------------------------------------

# Read bedgraphs, and add missing info
read_bdg = function(bdg, col_list = NULL, genome = NULL, verbose = TRUE,
                    strand_collapse = FALSE, fill_cpgs = TRUE,
                    contigs = contigs, synced_coordinates = synced_coordinates,
                    file_uncovered = NULL, zero_based = TRUE) {
  
  chr <- M <- U <- . <- NULL
  message(paste0("-Processing:    ", basename(bdg)))
  bdg_dat = suppressWarnings(data.table::fread(file = bdg, sep = "\t",
                                               colClasses = col_list$col_classes,
                                               verbose = FALSE,
                                               showProgress = FALSE))
  colnames(bdg_dat)[col_list$col_idx] = names(col_list$col_idx)
  
  if ("beta" %in% colnames(bdg_dat)) {
    if (nrow(bdg_dat) < 1000) {
      sample_row_idx = 1:nrow(bdg_dat)
      max_beta = max(bdg_dat[, beta], na.rm = TRUE)
    } else {
      # Choose 1000 random beta values
      sample_row_idx = sample(x = seq_len(nrow(bdg_dat)), size = 1000,
                              replace = FALSE)
      max_beta = max(bdg_dat[sample_row_idx, beta], na.rm = TRUE)
      rm(sample_row_idx)
    }
    if (max_beta > 1) {
      bdg_dat[, `:=`(beta, beta/100)]
      if (verbose) {
        message("--Note:         Converted beta values from percent to fractions")
      }
      rm(max_beta)
    }
    gc(verbose = FALSE)
  }
  
  if (!is.null(col_list$fix_missing)) {
    for (cmd in col_list$fix_missing) {
      bdg_dat[, eval(parse(text = cmd))]
    }
  }
  
  bdg_dat[, `:=`(chr, as.character(chr))]
  bdg_dat[, `:=`(start, as.integer(start))]
  
  if (zero_based) {
    # Bring bedgraphs to 1-based cordinate
    bdg_dat[, `:=`(start, start + 1)]
  }
  
  # Check for contig prefixes and add them if necessary
  if (nrow(bdg_dat) < 1000) {
    sample_row_idx = sample(x = seq_len(nrow(bdg_dat)),
                            size = as.integer(nrow(bdg_dat)/2), replace = FALSE)
  } else {
    sample_row_idx = sample(x = seq_len(nrow(bdg_dat)),
                            size = 1000, replace = FALSE)
  }
  if (grepl(pattern = "chr", x = genome[1, chr]) != any(grepl(pattern = "chr",
                                                              x = bdg_dat[sample_row_idx, chr]))) {
    if (grepl(pattern = "chr", x = genome[1, chr])) {
      bdg_dat[, `:=`(chr, paste0("chr", chr))]
    } else if (grepl(pattern = "chr", x = bdg_dat[1, chr])) {
      bdg_dat[, `:=`(chr, paste0("chr", chr))]
    } else {
      stop("Prefix mismatch between provided CpGs and bedgraphs")
    }
  }
  
  if (!is.null(contigs)) {
    bdg_dat = bdg_dat[chr %in% as.character(contigs)]
  }
  
  if (synced_coordinates) {
    bdg_dat = bdg_dat[strand == "-", `:=`(start, start + 1L)]
  }
  
  data.table::setkey(x = bdg_dat, "chr", "start")
  data.table::setkey(x = genome, "chr", "start")
  
  missing_cpgs = genome[!bdg_dat[, list(chr, start)], on = c("chr", "start")]
  
  # Write missing CpGs to an op_dir
  if (!is.null(file_uncovered) && nrow(missing_cpgs) > 0) {
    fwrite(x = missing_cpgs, file = paste0(file_uncovered,
                                           gsub("\\.[[:alnum:]]+(\\.gz)?$",
                                                "", basename(bdg)), "_uncovered.bed"),
           sep = "\t", row.names = FALSE)
  }
  
  if (verbose) {
    if (nrow(missing_cpgs) > 0) {
      message(paste0("--CpGs missing: ", format(nrow(missing_cpgs), big.mark = ",")))
    }
    # message(paste0('Missing ', format(nrow(missing_cpgs), big.mark =
    # ','), ' reference CpGs from: ', basename(bdg)))
  }
  if (nrow(missing_cpgs) > 0) {
    missing_cpgs[, `:=`(width, NULL)][, `:=`(beta, NA)][, `:=`(cov, NA)][, `:=`(M, NA)][, `:=`(U, NA)]
    bdg_dat = data.table::rbindlist(list(bdg_dat, missing_cpgs), use.names = TRUE,
                                    fill = TRUE)
    data.table::setkey(x = bdg_dat, "chr", "start")
  }
  # Better than identical(); seems to take couple of seconds but this is
  # crucial to make sure everything is in order
  is_identical = all.equal(target = bdg_dat[, .(chr, start)],
                           current = genome[, .(chr, start)], ignore.row.order = FALSE)
  
  if (is(is_identical, "character")) {
    non_ref_cpgs = bdg_dat[!genome[, list(chr, start)], on = c("chr",
                                                               "start")]
    message(paste0("--Non ref CpGs: ", format(nrow(non_ref_cpgs), big.mark = ","),
                   " [Removing them]"))
    bdg_dat = bdg_dat[genome[, list(chr, start)], on = c("chr", "start")]
    data.table::setkey(x = bdg_dat, "chr", "start")
    is_identical = all.equal(target = bdg_dat[, .(chr, start)],
                             current = genome[, .(chr, start)],
                             ignore.row.order = FALSE)
    if (is(is_identical, "character")) {
      stop("Something went wrong with filling up of uncovered CpG sites.")
    }
  }
  # Re-assign strand info from genome (since some bedgraphs have no
  # strand info, yet cover CpGs from both strands. i,e MethylDackel)
  bdg_dat[, `:=`(strand, genome$strand)]
  
  if (strand_collapse) {
    # If strand information needs to collapsed, bring start position of
    # crick strand to previous base (on watson base) and estimate new M, U
    # and beta values
    if (!all(c("M", "U") %in% names(bdg_dat))) {
      stop("strand_collapse works only when M and U are available!")
    }
    
    bdg_dat[, `:=`(start, ifelse(strand == "-", yes = start - 1, no = start))]
    bdg_dat = bdg_dat[, .(M = sum(M, na.rm = TRUE), U = sum(U, na.rm = TRUE)),
                      .(chr, start)]
    bdg_dat[, `:=`(cov, M + U)]
    bdg_dat[, `:=`(beta, M/cov)]
    bdg_dat[, `:=`(strand, "*")]
  }
  
  # data.table::set(bdg_dat, which(is.nan(bdg_dat[,beta])), 'beta', NA)
  # If coverage is 0, convert corresponding beta as well as coverage
  # values to NA
  data.table::set(bdg_dat, which(bdg_dat[, cov] == 0), c("cov", "beta"),
                  NA)
  bdg_dat = bdg_dat[, .(chr, start, beta, cov, strand)]
  
  bdg_genome_stat = bdg_dat[!is.na(beta), .(mean_meth = mean(beta),
                                            median_meth = median(beta),
                                            mean_cov = mean(cov),
                                            median_cov = median(cov))]
  bdg_chr_stat = bdg_dat[!is.na(beta), .(mean_meth = mean(beta),
                                         median_meth = median(beta),
                                         mean_cov = mean(cov),
                                         median_cov = median(cov)), .(chr)]
  bdg_ncpg_stat = bdg_dat[!is.na(beta), .N, .(chr)]
  
  return(list(bdg = bdg_dat, genome_stat = bdg_genome_stat,
              chr_stat = bdg_chr_stat,
              ncpg = bdg_ncpg_stat))
}

#--------------------------------------------------------------------------------------------------------------------------

# Process samples in batches. Batches are processed in vectorized
# manner (ideal for large number of samples)
vect_code_batch <- function(files, col_idx, batch_size, col_data = NULL,
                            genome = NULL, strand_collapse = FALSE, thr = 1, contigs = contigs,
                            synced_coordinates, file_uncovered = NULL, zero_based = TRUE) {
  . <- NULL
  batches <- split(files, ceiling(seq_along(files)/batch_size))
  batches_samp_names <- split(rownames(col_data), ceiling(seq_along(rownames(col_data))/batch_size))
  
  beta_mat_final <- data.table::data.table()
  cov_mat_final <- data.table::data.table()
  genome_stat_final <- data.table::data.table()
  chr_stat_final <- data.table::data.table()
  ncpg_final <- data.table::data.table()
  
  for (i in seq_along(batches)) {
    # browser()
    message(paste0("-Batch:         ", i, "/", length(batches)))
    batch_files <- batches[[i]]
    samp_names <- batches_samp_names[[i]]
    if (grepl("Windows", Sys.getenv("OS"))) {
      if (thr > 1) {
        warning("Windows doesn't support parallel processing. Setting n_threads to 1.")
      }
      bdgs <- lapply(batch_files, read_bdg, col_list = col_idx, genome = genome,
                     strand_collapse = strand_collapse, contigs = contigs, synced_coordinates = synced_coordinates,
                     file_uncovered = file_uncovered, zero_based = zero_based)
    } else {
      bdgs <- parallel::mclapply(batch_files, read_bdg, col_list = col_idx,
                                 genome = genome, strand_collapse = strand_collapse, mc.cores = thr,
                                 contigs = contigs, synced_coordinates = synced_coordinates,
                                 file_uncovered = file_uncovered, zero_based = zero_based)
    }
    names(bdgs) <- samp_names
    
    if (i == 1) {
      cov_mat_final <- data.frame(lapply(bdgs, function(x) x$bdg[,
                                                                 .(cov)]), stringsAsFactors = FALSE)
      beta_mat_final <- data.frame(lapply(bdgs, function(x) x$bdg[,
                                                                  .(beta)]), stringsAsFactors = FALSE)
      colnames(cov_mat_final) <- colnames(beta_mat_final) <- samp_names
      
      genome_stat_final <- data.table::rbindlist(lapply(bdgs, function(x) x$genome_stat),
                                                 use.names = TRUE, fill = TRUE, idcol = "Sample_Name")
      chr_stat_final <- data.table::rbindlist(lapply(bdgs, function(x) x$chr_stat),
                                              use.names = TRUE, fill = TRUE, idcol = "Sample_Name")
      ncpg_final <- data.table::rbindlist(lapply(bdgs, function(x) x$ncpg),
                                          use.names = TRUE, fill = TRUE, idcol = "Sample_Name")
    } else {
      cov_mat <- data.frame(lapply(bdgs, function(x) x$bdg[, .(cov)]),
                            stringsAsFactors = FALSE)
      beta_mat <- data.frame(lapply(bdgs, function(x) x$bdg[, .(beta)]),
                             stringsAsFactors = FALSE)
      colnames(cov_mat) <- colnames(beta_mat) <- samp_names
      cov_mat_final <- cbind(cov_mat_final, cov_mat)
      beta_mat_final <- cbind(beta_mat_final, beta_mat)
      
      genome_stat_final <- rbind(genome_stat_final, data.table::rbindlist(lapply(bdgs,
                                                                                 function(x) x$genome_stat), use.names = TRUE, fill = TRUE,
                                                                          idcol = "Sample_Name"))
      chr_stat_final <- rbind(chr_stat_final, data.table::rbindlist(lapply(bdgs,
                                                                           function(x) x$chr_stat), use.names = TRUE, fill = TRUE,
                                                                    idcol = "Sample_Name"))
      ncpg_final <- rbind(ncpg_final, data.table::rbindlist(lapply(bdgs,
                                                                   function(x) x$ncpg), use.names = TRUE, fill = TRUE, idcol = "Sample_Name"))
      
      rm(cov_mat)
      rm(beta_mat)
      gc()
    }
  }
  gc()
  ncpg_final <- data.table::dcast(data = ncpg_final, chr ~ Sample_Name,
                                  value.var = "N")
  
  return(list(beta_matrix = data.table::setDT(beta_mat_final), cov_matrix = data.table::setDT(cov_mat_final),
              genome_stat = genome_stat_final, chr_stat = chr_stat_final, ncpg = ncpg_final))
}

#--------------------------------------------------------------------------------------------------------------------------


# Use for loop for sample-by-sample processing, memory efficient, uses
# HDF5Array
non_vect_code <- function(files, col_idx, coldata, verbose = TRUE, genome = NULL,
                          h5temp = NULL, h5 = FALSE, strand_collapse = FALSE, contigs = contigs,
                          synced_coordinates, file_uncovered = NULL, zero_based = TRUE) {
  
  Sample_Name <- . <- chr <- NULL
  if (strand_collapse) {
    dimension <- as.integer(nrow(genome)/2)
  } else {
    dimension <- as.integer(nrow(genome))
  }
  
  if (h5) {
    if (is.null(h5temp)) {
      h5temp <- tempdir()
    }
    sink_counter <- 1
    while (any(c(paste0("M_sink_", sink_counter, ".h5"), paste0("cov_sink_",
                                                                sink_counter, ".h5")) %in% dir(h5temp))) {
      sink_counter <- sink_counter + 1
      
    }
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, 1L))
    
    M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                             dimnames = NULL, type = "double",
                                             filepath = file.path(h5temp, paste0("M_sink_", sink_counter, ".h5")), name = "M", level = 6)
    cov_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                               dimnames = NULL, type = "integer",
                                               filepath = file.path(h5temp, paste0("cov_sink_", sink_counter, ".h5")), name = "cov",
                                               level = 6)
  } else {
    beta_mat <- data.table::data.table()
    cov_mat <- data.table::data.table()
  }
  
  if (h5) {
    for (i in seq_along(files)) {
      if (i == 1) {
        b <- read_bdg(bdg = files[i], col_list = col_idx, genome = genome,
                      strand_collapse = strand_collapse, contigs = contigs, synced_coordinates = synced_coordinates,
                      file_uncovered = file_uncovered, zero_based = zero_based)
        
        DelayedArray::write_block(block = as.matrix(b$bdg[, .(beta)]),
                                  viewport = grid[[i]], x = M_sink)
        DelayedArray::write_block(block = as.matrix(b$bdg[, .(cov)]),
                                  viewport = grid[[i]], x = cov_sink)
        genome_stat_final <- b$genome_stat[, `:=`(Sample_Name,
                                                  rownames(coldata)[i])]
        chr_stat_final <- b$chr_stat[, `:=`(Sample_Name, rownames(coldata)[i])]
        ncpg_final <- b$ncpg[, `:=`(Sample_Name, rownames(coldata)[i])]
        rm(b)
        gc()
      } else {
        b <- read_bdg(bdg = files[i], col_list = col_idx, genome = genome,
                      strand_collapse = strand_collapse, contigs = contigs, synced_coordinates = synced_coordinates,
                      file_uncovered = file_uncovered, zero_based = zero_based)
        
        DelayedArray::write_block(block = as.matrix(b$bdg[, .(beta)]),
                                  viewport = grid[[i]], x = M_sink)
        DelayedArray::write_block(block = as.matrix(b$bdg[, .(cov)]),
                                  viewport = grid[[i]], x = cov_sink)
        genome_stat_final <- rbind(genome_stat_final, b$genome_stat[,
                                                                    `:=`(Sample_Name, rownames(coldata)[i])])
        chr_stat_final <- rbind(chr_stat_final, b$chr_stat[, `:=`(Sample_Name,
                                                                  rownames(coldata)[i])])
        ncpg_final <- rbind(ncpg_final, b$ncpg[, `:=`(Sample_Name,
                                                      rownames(coldata)[i])])
        rm(b)
        gc()
      }
    }
    ncpg_final <- data.table::dcast(data = ncpg_final, chr ~ Sample_Name,
                                    value.var = "N")
    return(list(beta_matrix = as(M_sink, "HDF5Array"), cov_matrix = as(cov_sink, "HDF5Array"),
                genome_stat = genome_stat_final, chr_stat = chr_stat_final,
                ncpg = ncpg_final))
  } else {
    for (i in seq_along(files)) {
      if (i == 1) {
        b <- read_bdg(bdg = files[i], col_list = col_idx, genome = genome,
                      strand_collapse = strand_collapse, contigs = contigs,
                      synced_coordinates = synced_coordinates, file_uncovered = file_uncovered,
                      zero_based = zero_based)
        
        beta_mat <- b$bdg[, .(chr, start, beta)]
        cov_mat <- b$bdg[, .(chr, start, cov)]
        genome_stat_final <- b$genome_stat[, `:=`(Sample_Name,
                                                  rownames(coldata)[i])]
        chr_stat_final <- b$chr_stat[, `:=`(Sample_Name, rownames(coldata)[i])]
        ncpg_final <- b$ncpg[, `:=`(Sample_Name, rownames(coldata)[i])]
      } else {
        b <- read_bdg(bdg = files[i], col_list = col_idx, genome = genome,
                      strand_collapse = strand_collapse, contigs = contigs,
                      synced_coordinates = synced_coordinates, file_uncovered = file_uncovered,
                      zero_based = zero_based)
        
        beta_mat <- cbind(beta_mat, b$bdg[, .(beta)])
        cov_mat <- cbind(cov_mat, b$bdg[, .(cov)])
        genome_stat_final <- rbind(genome_stat_final, b$genome_stat[,
                                                                    `:=`(Sample_Name, rownames(coldata)[i])])
        chr_stat_final <- rbind(chr_stat_final, b$chr_stat[, `:=`(Sample_Name,
                                                                  rownames(coldata)[i])])
        ncpg_final <- rbind(ncpg_final, b$ncpg[, `:=`(Sample_Name,
                                                      rownames(coldata)[i])])
      }
      colnames(beta_mat)[ncol(beta_mat)] <- colnames(cov_mat)[ncol(cov_mat)] <- rownames(coldata)[i]
    }
    
    ncpg_final <- data.table::dcast(data = ncpg_final, chr ~ Sample_Name,
                                    value.var = "N")
    return(list(beta_matrix = beta_mat[, -(seq_len(2))], cov_matrix = cov_mat[, -(seq_len(2))],
                genome_stat = genome_stat_final, chr_stat = chr_stat_final,
                ncpg = ncpg_final))
  }
}

#--------------------------------------------------------------------------------------------------------------------------
# Parse genomic regions and convert them to key'd data.table
cast_ranges <- function(regions, set.key = TRUE) {
  chr <- . <- NULL
  if (is(regions, "GRanges")) {
    target_regions <- data.table::as.data.table(x = regions)
    target_regions[, `:=`(seqnames, as.character(seqnames))]
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    if (set.key){
    data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
    target_regions <- target_regions[, .(chr, start, end)]
  } else if (is(regions, "data.frame")) {
    if (all(c("chr", "start", "end") %in% colnames(regions))){
      regions <- regions[,c("chr", "start", "end")]
    } else {
      warning("Columns with names chr, start and end are not found. Assuming that the first three columns are chr, start and end.")
    }
    target_regions <- data.table::as.data.table(x = regions)
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    target_regions <- target_regions[, .(chr, start, end)]
    target_regions[, `:=`(chr, as.character(chr))]
    target_regions[, `:=`(start, as.numeric(as.character(start)))]
    target_regions[, `:=`(end, as.numeric(as.character(end)))]
    if (set.key){
    data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
  } else {
    stop("Invalid input class for regions. Must be a data.table or GRanges object")
  }
  
  target_regions
}

#--------------------------------------------------------------------------------------------------------------------------
# Get min/max/mean/median of a matrix
giveme_this <- function(mat, stat = "mean", na_rm = TRUE, ish5 = FALSE) {
  stat <- match.arg(arg = stat, choices = c("mean", "median", "min",
                                            "max", "sum"))
  
  if (ish5) {
    if (stat == "mean") {
      res <- DelayedMatrixStats::colMeans2(mat, na.rm = na_rm)
    } else if (stat == "median") {
      res <- DelayedMatrixStats::colMedians(mat, na.rm = na_rm)
    } else if (stat == "min") {
      res <- colMins(mat, na.rm = na_rm)
    } else if (stat == "max") {
      res <- colMaxs(mat, na.rm = na_rm)
    } else if (stat == "sum") {
      res <- DelayedMatrixStats::colSums2(mat, na.rm = na_rm)
    }
  } else {
    if (stat == "mean") {
      res <- matrixStats::colMeans2(mat, na.rm = na_rm)
    } else if (stat == "median") {
      res <- matrixStats::colMedians(mat, na.rm = na_rm)
    } else if (stat == "min") {
      res <- matrixStats::colMins(mat, na.rm = na_rm)
    } else if (stat == "max") {
      res <- matrixStats::colMaxs(mat, na.rm = na_rm)
    } else if (stat == "sum") {
      res <- matrixStats::colSums2(mat, na.rm = na_rm)
    }
  }
  
  res
}


#--------------------------------------------------------------------------------------------------------------------------
# Tiny script to get axis and limits
get_y_lims <- function(vec) {
  
  y_lims <- range(vec)
  y_at <- pretty(y_lims)
  
  if (y_at[1] > min(vec, na.rm = TRUE)) {
    y_at[1] <- min(vec, na.rm = TRUE)
  }
  if (y_at[length(y_at)] < max(vec, na.rm = TRUE)) {
    y_at[length(y_at)] <- max(vec, na.rm = TRUE)
  }
  y_lims <- range(y_at, na.rm = TRUE)
  
  list(y_lims = y_lims, y_at = y_at)
}