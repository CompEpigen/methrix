#' Extracts all CpGs from a genome
#' @param ref_genome BSgenome object or name of the installed BSgenome package. Example: BSgenome.Hsapiens.UCSC.hg19
#' @param bored Tell me a Chuck Norris joke while I wait for my CpG extraction. Default TRUE. jokes can be very explicit! you have been warned!
#' @importFrom BSgenome installed.genomes getBSgenome seqnames
#' @importFrom rjson fromJSON
#' @export
#' @return a data.table of CpG's

extract_CPGs = function(ref_genome = NULL, bored = TRUE){

  gnoms_installed = BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = gnoms_installed)

  if(is.null(ref_genome)){
    if(nrow(gnoms_installed) == 0){
      stop("Could not find any installed BSgenomes.\nUse BSgenome::available.genomes() for options.")
    }else{
      message("Found following BSgenome installtions. Use the required 'pkgname'.")
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
        message("Found following BSgenome installtions. Provide the correct 'pkgname'")
        print(gnoms_installed)
        stop()
      }
    }
  }

  if(bored){
    #Try to fetch a random joke from Internet Chuck Norris Database
    rec_url = paste("http://api.icndb.com/jokes/random/")
    #joke = rjson::fromJSON(file = rec_url)
    joke = suppressWarnings(try(expr = rjson::fromJSON(file = rec_url), silent = TRUE))
    if(class(joke) == "try-error"){
      joke = list(value = "fail")
    }
  }

  ref_genome = BSgenome::getBSgenome(genome = ref_genome)
  chrs = names(ref_genome)

  #Code borrwed from from: https://support.bioconductor.org/p/95239/
  message("Extracting CpGs..")
  if(bored){
    if(joke$type == "success"){
      message("Here is a Chuck Norris joke while you wait..")
      cat(joke$value$joke, sep = "\n")
    }
  }
  cgs = lapply(chrs, function(x) start(matchPattern("CG", ref_genome[[x]])))
  cpgs = do.call(c, lapply(seq_along(chrs), function(x) GRanges(names(ref_genome)[x], IRanges(cgs[[x]], width = 2))))
  cpgs = data.table::as.data.table(as.data.frame(cpgs, stringsAsFactors = FALSE))
  colnames(cpgs) = c("chr", "start", "end", "width", "strand")
  cpgs[, chr := as.character(chr)][, start := as.numeric(start)][, end := as.numeric(end)][, width := as.numeric(width)]
  data.table::setkey(x = cpgs, "chr", "start")
  message(paste0("Done.\nExtracted ", nrow(cpgs), " from ", length(chrs), " contigs."))

  return(cpgs)
}
