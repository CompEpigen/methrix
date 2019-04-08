
#Process all samples in one go (ideal for few number of samples)
merge_vect_code = function(files, col_idx, col_data, ideal = FALSE, genome = NULL){
  bdgs = lapply(files, read_bdg, col_list = col_idx, genome = genome, fill_cpgs = FALSE)
  names(bdgs) = rownames(col_data)
  if(ideal){
    cov_mat = data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE)
    beta_mat = data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE)
    colnames(cov_mat) = colnames(beta_mat) = rownames(col_data)
    cov_mat = cbind(bdgs[[1]][,.(chr, start)], cov_mat)
    beta_mat = cbind(bdgs[[1]][,.(chr, start)], beta_mat)
  }else{
    bdgs = data.table::rbindlist(l = bdgs, idcol = "sample")
    beta_mat = data.table::dcast(bdgs, chr + start ~ sample, value.var = "beta")
    cov_mat = data.table::dcast(bdgs, chr + start ~ sample, value.var = "cov", fill = 0)
  }
  rm(bdgs)
  gc()
  return(list(beta_matrix = beta_mat, cov_matrix = cov_mat))
}

#--------------------------------------------------------------------------------------------------------------------------

#Process samples in batches. Batches are processed in vectorized manner (ideal for large number of samples)
merge_vect_code_batch = function(files, col_idx, batch_size, ideal = FALSE, col_data = NULL, genome = NULL){
  batches = split(files, ceiling(seq_along(files)/batch_size))
  batches_samp_names = split(rownames(col_data), ceiling(seq_along(rownames(col_data))/batch_size))

  beta_mat_final = data.table::data.table()
  cov_mat_final = data.table::data.table()

  for(i in seq_along(batches)){
    message(paste0("Processing batch ",  i , " of ", length(batches)))
    batch_files = batches[[i]]
    samp_names = batches_samp_names[[i]]
    bdgs = lapply(batch_files, read_bdg, col_list = col_idx, genome = genome, fill_cpgs = FALSE)
    names(bdgs) = samp_names
    if(ideal){
      if(i == 1){
        cov_mat_final = cbind(bdgs[[1]][,.(chr, start)],
                              data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE))
        beta_mat_final = cbind(bdgs[[1]][,.(chr, start)],
                               data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE))
        colnames(cov_mat_final) = colnames(beta_mat_final) = c("chr", "start", samp_names)
      }else{
        cov_mat = data.frame(lapply(bdgs, function(x) x[,.(cov)]), stringsAsFactors = FALSE)
        beta_mat = data.frame(lapply(bdgs, function(x) x[,.(beta)]), stringsAsFactors = FALSE)
        colnames(cov_mat) = colnames(beta_mat) = samp_names
        cov_mat_final = cbind(cov_mat_final, cov_mat)
        beta_mat_final = cbind(beta_mat_final, beta_mat)
      }
    }else{
      bdgs = data.table::rbindlist(l = bdgs, idcol = "sample")
      if(i == 1){
        beta_mat_final = data.table::dcast(bdgs, chr + start ~ sample, value.var = "beta")
        cov_mat_final = data.table::dcast(bdgs, chr + start ~ sample, value.var = "cov", fill = 0)
      }else{
        beta_mat_final = merge(beta_mat_final,
                               data.table::dcast(bdgs, chr + start ~ sample, value.var = "beta"),
                               by = c("chr", "start"), all = TRUE)
        cov_mat_final = merge(cov_mat_final,
                              data.table::dcast(bdgs, chr + start ~ sample, value.var = "cov", fill = 0),
                              by = c("chr", "start"), all = TRUE)
      }
    }
  }

  return(list(beta_matrix = data.table::setDT(beta_mat_final), cov_matrix = data.table::setDT(cov_mat_final)))
}

#--------------------------------------------------------------------------------------------------------------------------

#Use for loop for sample-by-sample processing. (Slow af!)
merge_non_vect_code = function(files, col_idx, coldata, verbose = TRUE, ideal = FALSE, genome = NULL){
  beta_mat = data.table::data.table()
  cov_mat = data.table::data.table()

  for(i in seq_along(files)){
    if(verbose){
      message("Processing: ", files[i])
    }
    if(ideal){
      if(i == 1){
        b = read_bdg(bdg = files[i], col_list = col_idx, genome = genome, fill_cpgs = FALSE)
        beta_mat = b[,.(chr, start, beta)]
        cov_mat = b[,.(chr, start, cov)]
      }else{
        b = read_bdg(bdg = files[i], col_list = col_idx, genome = genome, fill_cpgs = FALSE)
        beta_mat = cbind(beta_mat, b[,.(beta)])
        cov_mat = cbind(cov_mat, b[,.(cov)])
      }
      colnames(beta_mat)[ncol(beta_mat)] = colnames(cov_mat)[ncol(cov_mat)] = rownames(coldata)[i]
    }else{
      if(i == 1){
        b = read_bdg(bdg = files[i], col_list = col_idx, genome = genome, fill_cpgs = FALSE)
        #b[, id := paste0(chr, ":", start)]
        beta_mat = b[,.(chr, start, beta)]
        cov_mat = b[,.(chr, start, cov)]
        colnames(beta_mat)[ncol(beta_mat)] = colnames(cov_mat)[ncol(cov_mat)] = rownames(coldata)[i]
      }else{
        b = read_bdg(bdg = files[i], col_list = col_idx, genome = genome, fill_cpgs = FALSE)
        #b[, id := paste0(chr, ":", start)]
        beta_mat = merge(beta_mat, b[,.(chr, start, beta)], by = c("chr", "start"), all = TRUE)
        cov_mat = merge(cov_mat, b[,.(chr, start, cov)], by = c("chr", "start"), all = TRUE)
        colnames(beta_mat)[ncol(beta_mat)] = colnames(cov_mat)[ncol(cov_mat)] = rownames(coldata)[i]
      }
    }
  }
  mat_list = list(beta_matrix = beta_mat, cov_matrix = cov_mat)
}
