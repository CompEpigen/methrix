#' Versatile BedGraph reader.
#' @details Reads BedGraph files and generates methylation and coverage matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files bedgraph files.
#' @param pipeline Default NULL. Can be 'Bismark' or 'MethylDeckal' or 'MethylcTools'.
#' If not known use idx arguments for manual column assignments.
#' @param zero_based Are bedgraph regions zero based ? Default TRUE
#' @param stranded Default FALSE
#' @param collapse_strands If TRUE collapses CpGs on different crick strand into watson. Deafult FALSE
#' @param ref_cpgs BSgenome object, or name of the installed BSgenome package, or an output from \code{\link{extract_CPGs}}.
#' Example: BSgenome.Hsapiens.UCSC.hg19
#' @param ref_build reference genome for bedgraphs. Default NULL. Only used for additional details. Doesnt affect in any way.
#' @param contigs contigs to restrict genomic CpGs to. Default all autosomes and allosomes - ignoring extra contigs.
#' @param vect To use vectorized code. Default FALSE. Set to TRUE if you don't have large number of BedGraph files.
#' @param vect_batch_size Default NULL. Process samples in batches. Applicable only when vect = TRUE
#' @param coldata An optional DataFrame describing the samples. Row names, if present, become the column names of the matrix.
#' If NULL, then a DataFrame will be created with basename of files used as the row names.
#' @param chr_idx column index for chromosome in bedgraph files
#' @param start_idx column index for start position in bedgraph files
#' @param end_idx column index for end position in bedgraph files
#' @param beta_idx column index for beta values in bedgraph files
#' @param M_idx column index for read counts supporting Methylation in bedgraph files
#' @param U_idx column index for read counts supporting Un-methylation in bedgraph files
#' @param strand_idx column index for strand information in bedgraph files
#' @param cov_idx column index for total-coverage in bedgraph files
#' @param synced_coordinates Are the start and end coordinates of a stranded
#' bedgraph are synchronized between + and - strands? Possible values: FALSE (default),
#' TRUE if the start coordinates are the start coordinates of the C on the plus strand.
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param h5 Should the coverage and methylation matrices be stored as 'HDF5Array'
#' @param h5_dir directory to store H5 based object
#' @param h5temp temporary directory to store hdf5
#' @param verbose Be little chatty ? Default TRUE.
#' @export
#' @return An object of class \code{\link{methrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import parallel
#' @import DelayedMatrixStats
#' @import SummarizedExperiment DelayedArray HDF5Array
#' @examples
#'\dontrun{
#'bdg_files = list.files(path = system.file('extdata', package = 'methrix'),
#'pattern = '*\\.bedGraph\\.gz$', full.names = TRUE)
#' hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
#' meth = methrix::read_bedgraphs( files = bdg_files, ref_cpgs = hg19_cpgs,
#' chr_idx = 1, start_idx = 2, M_idx = 3, U_idx = 4,
#' stranded = FALSE, zero_based = FALSE, collapse_strands = FALSE)
#'}
#'

read_bedgraphs <- function(files = NULL, pipeline = NULL, zero_based = TRUE,
    stranded = FALSE, collapse_strands = FALSE, ref_cpgs = NULL, ref_build = NULL,
    contigs = NULL, vect = FALSE, vect_batch_size = NULL, coldata = NULL,
    chr_idx = NULL, start_idx = NULL, end_idx = NULL, beta_idx = NULL,
    M_idx = NULL, U_idx = NULL, strand_idx = NULL, cov_idx = NULL, synced_coordinates = FALSE,
    n_threads = 1, h5 = FALSE, h5_dir = NULL, h5temp = NULL, verbose = TRUE) {

    contig <- chr <- genome_contigs <- . <- NULL
    if (is.null(files)) {
        stop("Missing input files.", call. = FALSE)
    }

    if (is.null(ref_cpgs)) {
        stop("Missing genome. Please provide a valid genome name", call. = FALSE)
    }

    if (collapse_strands) {
        if (!stranded) {
            stop("collapse_strands works only when stranded = TRUE ")
        }
    }

    if (vect && is.null(vect_batch_size)) {
        vect_batch_size <- length(files)
    }

    start_proc_time <- proc.time()

    # Final aim is to bring input data to the following order: chr start
    # end beta cov starnd <rest..>
    message(paste0("----------------------------", "\n"))
    if (is.null(pipeline)) {
        message(paste0("-Preset:        Custom \n"))
        col_idx <- parse_source_idx(chr = chr_idx, start = start_idx, end = end_idx,
            beta = beta_idx, cov = cov_idx, strand = strand_idx, n_meth = M_idx,
            n_unmeth = U_idx, verbose = verbose)
        col_idx$col_classes <- NULL
    } else {
        pipeline <- match.arg(arg = pipeline, choices = c("Bismark_cov",
            "MethylDackel", "MethylcTools"))
        if (verbose) {
            message(paste0("-Preset:        ", pipeline, "\n"))
        }

        col_idx <- get_source_idx(protocol = pipeline)

        if (any(pipeline %in% c("Bismark_cov"))) {
            if (zero_based) {
                message("*BismarkCov files are one based. You may want to re-run with zero_based=FALSE\n")
            }
        }
    }

    if (is.null(contigs)) {
        # Work with only main contrigs (either with chr prefix - UCSC style, or
        # ensemble style)
        contigs <- c(paste0("chr", c(seq_len(22), "X", "Y", "M")), seq_len(22),
            "X", "Y", "MT")
    }

    # Extract CpG's
    conig_lens <- NA
    if (is(ref_cpgs, "list") & all(names(ref_cpgs) == c("cpgs", "contig_lens",
        "release_name"))) {
        conig_lens <- ref_cpgs$contig_lens[contig %in% contigs]
        genome <- data.table::copy(x = ref_cpgs$cpgs)
    } else if (any(is(ref_cpgs[1], "character"), is(ref_cpgs[1], "BSgenome"))) {
        genome <- extract_CPGs(ref_genome = ref_cpgs)
        conig_lens <- genome$contig_lens
        genome <- genome$cpgs
    } else if (is(ref_cpgs[1], "data.frame")) {
        genome <- data.table::copy(x = ref_cpgs)
        data.table::setkey(x = genome, chr, start)
    } else {
        stop("Could not figure out the genome class.")
    }

    if (verbose) {
        message(paste0("-CpGs raw:      ", format(nrow(genome), big.mark = ","),
            "\n"))
    }
    genome <- genome[chr %in% as.character(contigs)]

    if (nrow(genome) == 0) {
        message("No more CpG's left after subsetting for contigs. It appears provided contig names do not match to the BSgenome.")
        message("Contigs provied:")
        print(contigs)
        message("Contigs from BSGenome:")
        print(genome_contigs)
        stop(call. = FALSE)
    }

    if (verbose) {
        message(paste0("-CpGs filtered: ", format(nrow(genome), big.mark = ","),
            "\n"))
    }

    # check it with the strand column
    if (stranded) {
        genome_plus <- data.table::copy(genome)
        genome_plus[, `:=`(strand, "+")]
        genome[, `:=`(start, start + 1)]
        genome[, `:=`(strand, "-")]
        genome <- data.table::rbindlist(list(genome, genome_plus), use.names = TRUE)
        data.table::setkeyv(genome, cols = c("chr", "start"))
        message(paste0("-CpGs stranded: ", format(nrow(genome), big.mark = ","),
            "\n"))
        rm(genome_plus)
        gc()
    }

    # Set colData
    if (is.null(coldata)) {
        coldata <- data.frame(row.names = unlist(data.table::tstrsplit(x = basename(files),
            split = "\\.", keep = 1)), stringsAsFactors = FALSE)
    } else {
        if (length(files) != nrow(coldata)) {
            stop("Number of samples in coldata does not match the number of input files.")
        }
        if (any(make.names(rownames(coldata), unique = TRUE) != rownames(coldata))) {
            modified <- which(make.names(rownames(coldata), unique = TRUE) !=
                rownames(coldata))
            message("The sample names contained a non-valid character or were duplicated.
                The following changes were made:\n")
            message(paste(paste(rownames(coldata)[modified], make.names(rownames(coldata),
                unique = TRUE)[modified], sep = " => "), collapse = " \n "),
                "\n")
            rownames(coldata) <- make.names(rownames(coldata), unique = TRUE)

        }

    }

    message(paste0("----------------------------", "\n"))

    # Summarize bedgraphs and create a matrix
    if (vect) {
        mat_list <- vect_code_batch(files = files, col_idx = col_idx, batch_size = vect_batch_size,
            col_data = coldata, genome = genome, strand_collapse = collapse_strands,
            thr = n_threads, contigs = contigs, synced_coordinates = synced_coordinates,
            zero_based = zero_based)
    } else {
        mat_list <- non_vect_code(files = files, col_idx = col_idx, coldata = coldata,
            strand_collapse = collapse_strands, verbose = verbose, genome = genome,
            h5 = h5, h5temp = h5temp, contigs = contigs, synced_coordinates = synced_coordinates,
            zero_based = zero_based)
    }

    if (nrow(mat_list$beta_matrix) != nrow(mat_list$cov_matrix)) {
        stop("Discrepancies in dimensions of coverage and beta value matrices.")
    }

    # Finally collapse ref CpGs strands
    if (collapse_strands) {
        genome <- genome[strand %in% "+"]
        genome[, `:=`(strand, "*")]
    }

    ref_cpgs_chr <- genome[, .N, chr]

    descriptive_stats <- list(genome_stat = mat_list$genome_stat, chr_stat = mat_list$chr_stat,
        n_cpgs_covered = mat_list$ncpg)

    if (is.null(ref_build)) {
        ref_build <- ref_cpgs$release_name
    }

    m_obj <- create_methrix(beta_mat = mat_list$beta_mat, cov_mat = mat_list$cov_matrix,
        cpg_loci = genome[, .(chr, start, strand)], is_hdf5 = h5, genome_name = ref_build,
        col_data = coldata, h5_dir = h5_dir, ref_cpg_dt = ref_cpgs_chr,
        chrom_sizes = conig_lens, desc = descriptive_stats)

    rm(genome)
    gc()

    message("-Finished in:  ", data.table::timetaken(start_proc_time), "\n")
    return(m_obj)
}
