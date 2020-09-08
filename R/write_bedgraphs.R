#' Writes bedGraphs from methrix object
#' @param m \code{\link{methrix}} object
#' @param output_dir Output directory name where the files should be saved.
#' If \code{NULL} creats a \code{tempdir}
#' @param n_thr Default 4.
#' @param rm_NA remove NAs
#' @param force forces to create files if they are existing
#' @param compress Whether to compress the output. Default TRUE
#' @param SeqStyle Default `UCSC` with `chr` prefix.
#' @examples
#' data('methrix_data')
#' write_bedgraphs(m = methrix_data, output_dir = './temp')
#' @return writes bedgraph files to output
#' @export

write_bedgraphs <- function(m, output_dir = NULL, rm_NA = TRUE, force = FALSE, 
    n_thr = 4, compress = TRUE, SeqStyle="UCSC") {
    
    
    if (!dir.exists(output_dir)) {
        dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
  mat_gr <- methrix::get_matrix(m = m, type = "M", add_loci = TRUE, in_granges = TRUE)
  
  #change SeqlevelsStyle
  GenomeInfoDb::seqlevelsStyle(mat_gr)<- SeqStyle
  mat <- as.data.table(mat_gr) 
  mat <- mat[, c("seqnames", "start", "end", "strand", rownames(colData(x = m))), 
             with = FALSE]
    
    if (is.null(output_dir)) {
        output_dir <- getwd()
    }
   
  #parameters of trackline for bedgraph 
  parameters <-c("color" = "255,0,0", 
                 "visibility"="full", 
                 "altColor" = "128,128,128",
                 "autoScale"="on", 
                 "viewLimits"="0:1", 
                 "windowingFunction"="mean")  
  parameters <- paste0(" ", paste(names(parameters), parameters, 
                                  sep = "=", collapse = " "))
  
  
    message("----------------------")
    message("*Writing bedGraphs:")
    op_bdgs <- lapply(seq_len(nrow(colData(m))), function(i) {
        mat_i <- mat[, c(seq_len(3), i + 4), with = FALSE]
        if (rm_NA) {
            mat_i <- mat_i[complete.cases(mat_i), , ]
        }
        
        if (compress) {
            op_bdg <- paste0(output_dir, "/", rownames(colData(m))[i], 
                ".bedGraph.gz")
        } else {
            op_bdg <- paste0(output_dir, "/", rownames(colData(m))[i], 
                ".bedGraph")
        }
        
        #add sample name to trackline
        header <- data.table(paste0('track type=bedGraph name="', rownames(colData(m))[i], '"', parameters))
        
        if (file.exists(op_bdg)) {
            if (force) {
                message(paste0("**Writing ", rownames(colData(m))[i]))
                colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
                data.table::fwrite(x = header, file = op_bdg, sep = "\t", append=FALSE, quote=FALSE,
                                   col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
                data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", append=TRUE,
                  col.names = FALSE, nThread = n_thr, scipen = 7, compress = "auto")
            } else {
                message(paste0("**File ", basename(op_bdg), " already exists. Skipped re-writing"))
            }
        } else {
            message(paste0("**Writing ", rownames(colData(m))[i]))
            colnames(mat_i) <- paste0("V", seq_len(ncol(mat_i)))
            data.table::fwrite(x = header, file = op_bdg, sep = "\t", col.names = FALSE, append=FALSE, quote=FALSE,
                               nThread = n_thr)
            data.table::fwrite(x = mat_i, file = op_bdg, sep = "\t", col.names = FALSE, append=TRUE,
                nThread = n_thr)
        }
        op_bdg
    })
    message("----------------------")
}
