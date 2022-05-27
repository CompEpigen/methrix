#' Writes bedGraphs from methrix object
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved.
#' If \code{NULL} creats a \code{tempdir}
#' @param n_thr Default 4.
#' @param rm_NA remove NAs
#' @param force forces to create files if they are existing
#' @param compress Whether to compress the output. Default TRUE
#' @param SeqStyle Default `UCSC` with `chr` prefix.
#' @param multiBed Default NULL. If provided a filename, a single bedGraph file with all samples included is generated.
#' @param metilene Default FALSE. If TRUE outputs bedgraphs in `metilene` format that can be directly used for DMR calling with `metilene`. This option works only when \code{multiBed = TRUE}.
#' @param phenoCol Default NULL. `condition` column from colData. Only applicable if \code{metilene = TRUE}
#' @param add_coverage Default FALSE. Whether to add a coverage column to the output. Only applicable if \code{multiBed = NULL}
#' @examples
#' data('methrix_data')
#' write_bedgraphs(m = methrix_data, output_dir = './temp')
#' #Export to metline format for DMR calling with metline
#' write_bedgraphs(m = methrix_data, output_dir = "./temp", rm_NA = FALSE, metilene = TRUE,multiBed = "metline_ip", phenoCol = "Condition")
#' @return writes bedgraph files to output
#' @export

write_bedgraphs <- function(m, output_dir = NULL, rm_NA = TRUE, force = FALSE, 
  n_thr = 4, compress = TRUE, SeqStyle="UCSC", multiBed = NULL, metilene = FALSE, 
  phenoCol = NULL, add_coverage = FALSE) {
      
  if (!dir.exists(output_dir)) {
        dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  } else if (is.null(output_dir)) {
    output_dir <- getwd()
  }
    
  mat_gr <- methrix::get_matrix(m = m, type = "M", add_loci = TRUE, in_granges = TRUE)
  GenomeInfoDb::seqlevelsStyle(mat_gr)<- SeqStyle
  mat <- as.data.table(mat_gr) 
  mat <- mat[, c("seqnames", "start", "end", "strand", rownames(colData(x = m))), 
             with = FALSE]


  if (!is.null(multiBed)) {    
    op_bdg = file.path(output_dir, paste0(
      multiBed, ".bedGraph", ifelse(compress, yes = ".gz", no = ""))
    )
    if (rm_NA) {
        mat = mat[complete.cases(mat),]
    }
    
    if (metilene) {
      if (is.null(phenoCol)) {
        stop("Please provide a value to phenoCol.")
      }
      else if (!phenoCol %in% colnames(colData(m))) {
        stop(phenoCol, 
          " is not a valid column name in colData().\n Available column name are: ",
          paste(colnames(colData(m)), collapse = ", "))
      }
      
      colnames(mat)[1:2] = c("chr", "pos")
      colnames(mat)[5:ncol(mat)] = paste(as.character(colData(m)[,phenoCol]),
        rownames(colData(m)), sep = "_")
      mat = mat[,c(1,2,5:ncol(mat)), with = FALSE]
    } else {
      colnames(mat)[1] = paste0("#", colnames(mat)[1])
    }

    message("*Writing ", op_bdg, " ")
    data.table::fwrite(x = mat, file = op_bdg, sep = "\t", na = ".", scipen = 7,
      nThread = n_thr, compress = "auto", showProgress = TRUE, quote = FALSE)
    
  } else {
    
    if (add_coverage){
      cov_gr <- methrix::get_matrix(m = m, type = "C", add_loci = TRUE, 
                                  in_granges = TRUE)
      GenomeInfoDb::seqlevelsStyle(cov_gr) <- SeqStyle
      cov_mat <- as.data.table(cov_gr)
      cov_mat <- cov_mat[, c("seqnames", "start", "end", "strand", 
        rownames(colData(x = m))), with = FALSE] 
    }
    op_bdgs <- lapply(seq_len(nrow(colData(m))), function(i) {
      mat_i <- mat[, c(seq_len(3), i + 4), with = FALSE]
      if (add_coverage) {
        mat_i <- cbind(mat_i, cov_mat[, c(i + 4), with = FALSE])
      }
      if (rm_NA) {
        mat_i <- mat_i[complete.cases(mat_i), , ]
      }
      
      op_bdg = file.path(output_dir, paste0(
        rownames(colData(m))[i], ".bedGraph", ifelse(compress, yes = ".gz", no = ""))
      )
      
      # bedGraph header line
      parameters <-c("color" = "255,0,0", 
                     "visibility"="full", 
                     "altColor" = "128,128,128",
                     "autoScale"="on", 
                     "viewLimits"="0:1", 
                     "windowingFunction"="mean")  
      parameters <- paste0(" ", paste(names(parameters), parameters, 
                                      sep = "=", collapse = " "))
      header <- data.table(paste0('track type=bedGraph name="', 
        rownames(colData(m))[i], '"', parameters))
      
      if (file.exists(op_bdg) & !force) {
        message(paste0("*File ", basename(op_bdg), " already exists. Skipped re-writing"))
        return(op_bdg)
      }

      message(paste0("*Writing ", rownames(colData(m))[i]))
      colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
      data.table::fwrite(x = header, file = op_bdg, sep = "\t", append=FALSE, quote=FALSE,
                         col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
      data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", append=TRUE,
                         col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
      return(op_bdg)
    })
    
  }
  
  message("----------------------")
}
