#' Writes bedGraphs from methrix object
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved. If \code{NULL} creats a \code{tempdir}
#' @param n_thr Default 4.
#' @param rm_NA remove NAs
#' @param force forces to create files if they are existing
#' @param compress Whether to compress the output. Default TRUE
#' @examples
#' data("methrix_data")
#' write_bedgraphs(m = methrix_data, output_dir = "./temp")
#' @return writes bedgraph files to output
#' @export

write_bedgraphs = function(m, output_dir = NULL, rm_NA = TRUE, force = FALSE, n_thr = 4, compress = TRUE){


  if(!dir.exists(output_dir)){
    dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  mat = methrix::get_matrix(m = m, type = "M", add_loci = TRUE)
  mat[,end := start + 1]
  mat = mat[,c("chr", "start", "end", "strand", rownames(colData(x = m))), with = FALSE]

  if(is.null(output_dir)){
    output_dir = getwd()
  }

  cat("----------------------\n")
  cat("*Writing bedGraphs:\n")
  op_bdgs = lapply(seq_len(nrow(colData(m))), function(i){
                mat_i = mat[,c(1:3, i+4), with = FALSE]
                if(rm_NA){
                  mat_i = mat_i[complete.cases(mat_i),,]
                }

                if(compress){
                  op_bdg = paste0(output_dir, "/", rownames(colData(m))[i], ".bedGraph.gz")
                }else{
                  op_bdg = paste0(output_dir, "/", rownames(colData(m))[i], ".bedGraph")
                }

                if(file.exists(op_bdg)){
                  if(force){
                    cat(paste0("**Writing ", rownames(colData(m))[i], "\n"))
                    colnames(mat_i) = paste0("V", 1:ncol(mat_i))
                    data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
                  }else{
                    cat(paste0("**File ", basename(op_bdg), " already exists. Skipped re-writing\n"))
                  }
                }else{
                  cat(paste0("**Writing ", rownames(colData(m))[i], "\n"))
                  colnames(mat_i) = paste0("V", 1:ncol(mat_i))
                  data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", col.names = FALSE, nThread = n_thr)
                }
                op_bdg
              })
  cat("----------------------\n")
}
