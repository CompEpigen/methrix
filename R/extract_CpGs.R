#' Extracts all CpGs from a genome
#' @param ref_genome BSgenome object or name of the installed BSgenome package. Example: BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome installed.genomes getBSgenome seqnames
#' @importFrom rjson fromJSON
#' @export
#' @return a list of data.table containing number of CpG's and contig lengths
#' @examples
#'\dontrun{
#' hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
#' }

extract_CPGs = function(ref_genome = NULL) {

    pkgname <- seqlengths <- chr <- NULL
    gnoms_installed = BSgenome::installed.genomes(splitNameParts = TRUE)
    data.table::setDT(x = gnoms_installed)

    if (is.null(ref_genome)) {
        if (nrow(gnoms_installed) == 0) {
            stop("Could not find any installed BSgenomes.\nUse BSgenome::available.genomes() for options.")
        } else {
            message("Found following BSgenome installations. Use the required 'pkgname'.")
            print(gnoms_installed)
            stop()
            # ref_genome = gnoms_installed[,pkgname][1]
        }
    } else {
        if (nrow(gnoms_installed[pkgname %in% ref_genome]) == 0) {
            message("Could not find BSgenome ", ref_genome)
            if (nrow(gnoms_installed) == 0) {
                stop("Could not find any installed BSgenomes either.\nUse BSgenome::available.genomes() for options.")
            } else {
                message("Found following BSgenome installations. Provide the correct 'pkgname'")
                print(gnoms_installed)
                stop()
            }
        }
    }
    requireNamespace(ref_genome, quietly = TRUE)


    ref_genome = BSgenome::getBSgenome(genome = ref_genome)
    ref_build = attributes(x = ref_genome)$provider_version
    chrom_sizes = data.table::data.table(contig = names(seqlengths(x = ref_genome)),
        length = seqlengths(x = ref_genome))
    chrs = names(ref_genome)

    message("-Extracting CpGs\n")

    # Code borrwed from from: https://support.bioconductor.org/p/95239/
    cgs = lapply(chrs, function(x) start(Biostrings::matchPattern("CG", ref_genome[[x]])))
    cpgs = do.call(c, lapply(seq_along(chrs), function(x) GenomicRanges::GRanges(names(ref_genome)[x],
        IRanges::IRanges(cgs[[x]], width = 2))))
    cpgs = data.table::as.data.table(as.data.frame(cpgs, stringsAsFactors = FALSE))
    colnames(cpgs) = c("chr", "start", "end", "width", "strand")
    cpgs[, `:=`(chr, as.character(chr))][, `:=`(start, as.numeric(start))][, `:=`(end,
        as.numeric(end))][, `:=`(width, as.numeric(width))]
    data.table::setkey(x = cpgs, "chr", "start")
    message(paste0("-Done. Extracted ", format(nrow(cpgs), big.mark = ","), " CpGs from ",
        length(chrs), " contigs.\n"))

    return(list(cpgs = cpgs, contig_lens = chrom_sizes, release_name = ref_build))
}
