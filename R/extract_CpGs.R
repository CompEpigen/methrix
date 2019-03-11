#' Extracts all CpGs from a genome
#' @param ref_genome BSgenome object or name of the installed BSgenome package. Example: BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome installed.genomes getBSgenome seqnames
#' @export
#' @return a data.table of CpG's

extract_CPGs = function(ref_genome = NULL){

  gnoms_installed = BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = gnoms_installed)

  if(is.null(ref_genome)){
    if(nrow(gnoms_installed) == 0){
      stop("Could not find any installed BSgenomes.\nUse BSgenome::available.genomes() for options.")
    }else{
      message("Found following BSgenome installtions. Use the required one.")
      print(gnoms_installed)
      stop()
      #ref_genome = gnoms_installed[,pkgname][1]
    }
  }else{
    if(nrow(gnoms_installed[pkgname %in% ref_genome]) == 0){
      message(paste0("Could not find BSgenome "), ref_genome)
      if(nrow(gnoms_installed) == 0){
        stop("Could not find any installed BSgenomes either.\nUse BSgenome::available.genomes() for options.")
      }else{
        message("Found following BSgenome installtions. Correct ref_genome argument if necessary.")
        print(gnoms_installed)
        stop()
      }
    }
  }

  ref_genome = BSgenome::getBSgenome(genome = ref_genome)
  chrs = names(ref_genome)

  #Code borrwed from from: https://support.bioconductor.org/p/95239/
  message("Extracting CpGs..")
  cgs = lapply(chrs, function(x) start(matchPattern("CG", ref_genome[[x]])))
  cpgs = do.call(c, lapply(seq_along(chrs), function(x) GRanges(names(ref_genome)[x], IRanges(cgs[[x]], width = 2))))
  cpgs = data.table::as.data.table(as.data.frame(cpgs))
  message(paste0("Done.\nExtracted ", nrow(cpgs), " from ", length(chrs), " contigs."))

  return(cpgs)
}
