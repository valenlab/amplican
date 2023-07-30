#' Automated analysis of CRISPR experiments.
#'
#' Main goals:
#' \enumerate{
#' \item Flexible pipeline for analysis of the CRISPR Mi-Seq or Hi-Seq data.
#' \item Compatible with GRanges and data.table style.
#' \item Precise quantification of mutation rates.
#' \item Prepare automatic reports as .Rmd files that are flexible
#' and open for manipulation.
#' \item Provide specialized plots for deletions, insertions, mismatches,
#' variants, heterogeneity of the reads.
#' }
#'
#' To learn more about amplican, start with the vignettes:
#' \code{browseVignettes(package = "amplican")}
#'
#' @docType package
#' @name amplican
#' @useDynLib amplican
#'
#' @import Rcpp ggthemes waffle knitr methods BiocGenerics Biostrings data.table
#' @importFrom Rcpp sourceCpp
#'
"_PACKAGE"
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0("Please consider supporting this software by citing:\n\n",
          "Labun et al. 2019\n",
          "Accurate analysis of genuine CRISPR editing events with ampliCan.\n",
          "Genome Res. 2019 Mar 8\n",
          "doi: 10.1101/gr.244293.118\n"))
}

amplicanPipe <- function(min_freq_default) {
  function(
    config, fastq_folder, results_folder, knit_reports = TRUE,
    write_alignments_format = "None", average_quality = 30,
    min_quality = 0, filter_n = FALSE, batch_size = 1e7, use_parallel = FALSE,
    scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
      match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA"),
    gap_opening = 25, gap_extension = 0, fastqfiles = 0.5,
    primer_mismatch = 2,
    donor_mismatch = 3, donor_strict = FALSE,
    PRIMER_DIMER = 30,
    event_filter = TRUE, cut_buffer = 5,
    promiscuous_consensus = TRUE, normalize = c("guideRNA", "Group"),
    min_freq = min_freq_default,
    continue = TRUE) {

    config <- normalizePath(config)
    fastq_folder <- normalizePath(fastq_folder)
    results_folder <- normalizePath(results_folder)

    message("Checking write access...")
    checkFileWriteAccess(results_folder)

    if (!continue) {
      message("continue is FALSE, removeing contents of results folder.")
      unlink(results_folder, recursive = TRUE)
      dir.create(results_folder, showWarnings = FALSE)
    }
    resultsFolder <- file.path(results_folder, "alignments")
    if (!dir.exists(resultsFolder)) {
      dir.create(resultsFolder)
    }
    rds_file <- file.path(resultsFolder, "AlignmentsExperimentSet.rds")

    if (file.exists(rds_file)) {
      message("Loading alignments...")
      aln <- readRDS(rds_file)
    } else {
      aln <- amplicanAlign(config = config,
                    fastq_folder = fastq_folder,
                    use_parallel = use_parallel,
                    average_quality = average_quality,
                    batch_size = batch_size,
                    scoring_matrix = scoring_matrix,
                    gap_opening = gap_opening,
                    gap_extension = gap_extension,
                    min_quality = min_quality,
                    filter_n = filter_n,
                    fastqfiles = fastqfiles,
                    primer_mismatch = primer_mismatch,
                    donor_mismatch = donor_mismatch,
                    donor_strict = donor_strict)
      message("Saving alignments...")
      saveRDS(aln, rds_file)
      message("Saved alignments.")
    }

    # save as other formats
    if (!"None" %in% write_alignments_format) {
      for (frmt in write_alignments_format) {
        aln_file_frmt <- file.path(resultsFolder,
                                   paste0("alignments.", frmt))
        if (!file.exists(aln_file_frmt)) {
          writeAlignments(aln, aln_file_frmt, frmt)
        }
      }
    }

    logFileName <- file.path(results_folder, "RunParameters.txt")
    if (!file.exists(logFileName)) {
      message("Saving parameters...")
      logFileConn <- file(logFileName, open = "at")
      writeLines(c(paste("Config file:        ", config),
                   paste("Average Quality:    ", average_quality),
                   paste("Minimum Quality:    ", min_quality),
                   paste("Filter N-reads:     ", filter_n),
                   paste("Batch size:         ", batch_size),
                   paste("Write Alignments:   ", toString(write_alignments_format)),
                   paste("Fastq files Mode:   ", fastqfiles),
                   paste("Gap Opening:        ", gap_opening),
                   paste("Gap Extension:      ", gap_extension),
                   paste("Consensus:          ", promiscuous_consensus),
                   paste("Normalize:          ", toString(normalize)),
                   paste("PRIMER DIMER buffer:", PRIMER_DIMER),
                   paste("Cut buffer:", cut_buffer),
                   "Scoring Matrix:"), logFileConn)
      utils::write.csv(scoring_matrix, logFileConn, quote = FALSE, row.names = TRUE)
      close(logFileConn)
    }


    un_file <- file.path(resultsFolder, "unassigned_reads.csv")
    if (!file.exists(un_file)) {
      message("Saving unassigned sequences...")
      unData <- unassignedData(aln)
      if (!is.null(unData)) data.table::fwrite(
        unData, file.path(resultsFolder, "unassigned_reads.csv"))
    }

    bd_file <- file.path(results_folder, "barcode_reads_filters.csv")
    if (!file.exists(bd_file)) {
      message("Saving barcode statistics...")
      data.table::fwrite(barcodeData(aln), bd_file)
    }
    cfgT <- experimentData(aln)

    re_file <- file.path(resultsFolder, "raw_events.csv")
    if (!file.exists(re_file)) {
      message("Translating alignments into events...")
      aln <- extractEvents(aln, use_parallel = use_parallel)
      message("Saving complete events - unfiltered...")
      data.table::fwrite(aln, re_file)
      message("Saved complete events - unfiltered.")
      aln <- data.table::as.data.table(aln)
    } else {
      message("Reading complete events - unfiltered.")
      aln <- fread(re_file)
    }

    seqnames <- read_id <- counts <- NULL
    if (dim(aln)[1] == 0) stop("There are no events.",
                               "Check whether you have correct primers in the config file.")

    efs_file <- file.path(resultsFolder, "events_filtered_shifted.csv")
    cs_file <- file.path(results_folder, "config_summary.csv")
    if (!file.exists(efs_file)) {
      aln$overlaps <- amplicanOverlap(aln, cfgT, cut_buffer = cut_buffer)
      aln$consensus <- if (fastqfiles <= 0.5) {
        amplicanConsensus(aln, cfgT, promiscuous = promiscuous_consensus)
      } else { TRUE }

      # filter events overlapping primers
      eOP <- findEOP(aln, cfgT)
      aln <- aln[!eOP, ]

      # find PRIMER DIMERS
      PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

      # summarize how many PRIMER DIMER reads per ID
      onlyPD <- aln[PD, ]
      onlyPD <- unique(onlyPD, by = c("seqnames", "read_id"))
      onlyPD <- data.table::as.data.table(onlyPD)
      summaryPD <- onlyPD[, list(counts  = sum(counts)), by = c("seqnames")]
      cfgT$PRIMER_DIMER <- 0
      cfgT$PRIMER_DIMER[match(summaryPD$seqnames, cfgT$ID)] <- summaryPD$counts

      # apply filter - remove all events that come from PD infected reads
      aln <- aln[!onlyPD, on = list(seqnames, read_id)]

      # alignment event filter
      cfgT$Low_Score <- 0
      if (event_filter) {
        for (i in seq_len(dim(cfgT)[1])) {
          aln_id <- aln[seqnames == cfgT$ID[i], ]
          if (dim(aln_id)[1] == 0 | cfgT$Donor[i] != "") next()
          onlyBR <- aln_id[findLQR(aln_id), ]
          onlyBR <- unique(onlyBR, by = "read_id")
          cfgT[i, "Low_Score"] <- sum(onlyBR$counts)
          aln <- aln[!(aln$seqnames == cfgT$ID[i] &
                         aln$read_id %in% onlyBR$read_id), ]
        }
      }
      cfgT$Reads_Filtered <- cfgT$Reads - cfgT$PRIMER_DIMER - cfgT$Low_Score

      if (donor_strict) {
        for (i in seq_len(dim(cfgT)[1])) {
          aln[seqnames == cfgT$ID[i], ]$readType <-
            is_hdr_strict_by_id(aln[seqnames == cfgT$ID[i], ],
                                get_seq(cfgT, cfgT$ID[i]),
                                get_seq(cfgT, cfgT$ID[i], "Donor"),
                                scoring_matrix, gap_opening,
                                gap_extension, donor_mismatch)

        }
      }

      # shift to relative (most left UPPER case is position 0)
      message("Shifting events as relative...")
      data.table::setDF(aln)
      aln <- data.frame(amplicanMap(aln, cfgT), stringsAsFactors = FALSE)
      message("Saving shifted events - filtered...")
      data.table::fwrite(aln, efs_file)
      message("Saved shifted events - filtered.")
      # revert guides to 5'-3'
      cfgT$guideRNA[cfgT$Direction] <- revComp(cfgT$guideRNA[cfgT$Direction])
      data.table::fwrite(cfgT, cs_file)
    } else {
      message("Reading shifted events - filtered.")
      aln <- fread(efs_file)
      cfgT <- fread(cs_file)
    }

    efsn_file <- file.path(resultsFolder,
                           "events_filtered_shifted_normalized.csv")
    if (!file.exists(efsn_file)) {
      message("Normalizing events...")
      # we remove all N as they are just noise from poor sequencing
      aln <- aln[aln$replacement != "N", ]
      aln <- amplicanNormalize(aln, cfgT, min_freq = min_freq, add = normalize)
      message("Saving normalized events...")
      data.table::fwrite(aln, efsn_file)
      message("Saved normalized events.")
    } else {
      message("Reading normalized events.")
      aln <- fread(efsn_file)
    }

    # summarize
    cfgT <- amplicanSummarize(aln[aln$consensus & aln$overlaps, ], cfgT)
    data.table::fwrite(
      cfgT[, c("ID", "Barcode", "Forward_Reads_File", "Reverse_Reads_File",
               "Group", "guideRNA", "Found_Guide", "Control", "Forward_Primer",
               "Reverse_Primer", "Direction", "Amplicon", "Donor", "fwdPrPosEnd",
               "rvePrPos", "Reads", "PRIMER_DIMER", "Low_Score",
               "Reads_Filtered", "Reads_Del", "Reads_In",
               "Reads_Edited", "Reads_Frameshifted", "HDR")], cs_file)

    # reports
    reportsFolder <- file.path(results_folder, "reports")
    if (!dir.exists(reportsFolder)) {
      dir.create(reportsFolder)
    } else {
      unlink(reportsFolder, recursive = T)
    }

    message(paste0("Making reports... \nDue to high quality ",
                   "figures, it is time consuming. Use .Rmd templates for ",
                   "more control."))
    amplicanReport(results_folder,
                   knit_reports = knit_reports,
                   cut_buffer = cut_buffer,
                   report_files = file.path(reportsFolder,
                                            c("id_report",
                                              "barcode_report",
                                              "group_report",
                                              "guide_report",
                                              "amplicon_report",
                                              "index")))
    message("Finished.")
    invisible(results_folder)
  }
}


#' Wraps main package functionality into one function.
#'
#' amplicanPipeline is convenient wrapper around all functionality of the
#' package with the most robust settings. It will generate all results in the
#' \code{result_folder} and also knit prepared reports into 'reports' folder.
#' @param results_folder (string) Where do you want to store results?
#' The package will create files in that folder so make sure you have writing
#' permissions.
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}.
#' Configuration file can contain additional columns, but first 11 columns
#' have to follow the example config specification.
#' @param fastq_folder (string) Path to FASTQ files. If not specified,
#' FASTQ files should be in the same directory as config file.
#' @param knit_reports (boolean) whether function should "knit" all
#' reports automatically for you (it is time consuming, be patient), when false
#' reports will be prepared, but not knitted
#' @param use_parallel (boolean) Set to TRUE, if you have registered
#' multicore back-end.
#' @param average_quality (numeric) The FASTQ file have a quality for each
#' nucleotide, depending on sequencing technology there exist many formats.
#' This package uses \code{\link[ShortRead]{readFastq}} to parse the reads.
#' If the average quality of the reads fall below value of
#' \code{average_quality} then sequence is filtered. Default is 0.
#' @param min_quality (numeric)  Similar as in average_quality, but depicts
#' the minimum quality for ALL nucleotides in given read. If one of nucleotides
#' has quality BELLOW \code{min_quality}, then the sequence is filtered.
#' Default is 20.
#' @param filter_n (boolean)  Whether to filter out reads containing N base.
#' @param batch_size (numeric) How many reads to analyze at a time? Needed for
#' filtering of large fastq files.
#' @param write_alignments_format (character vector) Whether
#' \code{amplicanPipeline} should write alignments results to separate files.
#' Alignments are also always saved as .rds object of
#' \code{\link{AlignmentsExperimentSet}} class.
#' Possible options are:
#' \itemize{
#'  \item{"fasta"}{ outputs alignments in fasta format where header indicates
#' experiment ID, read id and number of reads}
#'  \item{"txt"}{ simple format, read information followed by forward read and
#'  amplicon sequence followed by reverse read with its amplicon sequence
#'  eg.: \cr
#' \preformatted{
#' ID: ID_1 Count: 7
#' ACTGAAAAA--------
#' ACTG-----ACTGACTG
#'
#' ------G-ACTG
#' ACTGACTGACTG
#' }}
#' \item{"None"}{ Don't write any alignments to files.}
#' \item{c("fasta", "txt")}{ There are also possible combinations of
#' above formats, pass a vector to get alignments in multiple formats.}
#' }
#' @param scoring_matrix (matrix) Default is 'NUC44'. Pass desired matrix using
#' \code{\link{nucleotideSubstitutionMatrix}}.
#' @param gap_opening (numeric) The opening gap score.
#' @param gap_extension (numeric) The gap extension score.
#' @param fastqfiles (numeric) Normally you want to use both FASTQ files. But in
#' some special cases, you may want to use only the forward file, or only
#' the reverse file. Possible options:
#' \itemize{
#'  \item{0}{ Use both FASTQ files.}
#'  \item{0.5}{ Use both FASTQ files, but only for one of the reads (forward or
#'  reverse) is required to have primer perfectly matched to sequence - eg. use
#'  when reverse reads are trimmed of primers, but forward reads have forward
#'  primer in the sequence.}
#'  \item{1}{ Use only the forward FASTQ file.}
#'  \item{2}{ Use only the reverse FASTQ file.}
#' }
#' @param primer_mismatch (numeric) Decide how many mismatches are allowed
#' during primer matching of the reads, that groups reads by experiments.
#' When \code{primer_mismatch = 0} no mismatches are allowed, which can increase
#' number of unasssigned read.
#' @param donor_mismatch (numeric) How many events of length 1 (mismatches,
#' deletions and insertions of length 1) are allowed when aligning toward the
#' donor template. This parameter is only used when donor template is specified.
#' The higher the parameter the less strict will be algorithm accepting read as
#' HDR. Set to 0 if only perfect alignments to the donor template marked as HDR,
#' unadvised due to error rate of the sequencers.
#' @param donor_strict (logical) Applies more strict algorithm for HDR detection.
#' Only these reads that have all of the donor events will count as HDR. Tolerates `donor_mismatch`
#' level of noise, but no indels are allowed. Use this when your reads should span over the
#' whole window of the donor events. Might be more time consuming.
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @param event_filter (logical) Whether detection of offtarget reads,
#' should be enabled.
#' @param cut_buffer The number of bases by which extend expected cut sites
#' (specified as UPPER case letters in the amplicon) in 5' and 3' directions.
#' @param promiscuous_consensus (boolean) Whether rules of
#' \code{\link{amplicanConsensus}} should be \code{promiscuous}. When
#' promiscuous, we allow indels that have no confirmation on the other strand.
#' @param normalize (character vector) If column 'Control' in config table
#' has all FALSE/0 values then normalization is skipped. Otherwise,
#' normalization is strict, which means events that are
#' found in 'Control' TRUE group will be removed in 'Control' FALSE group.
#' This parameter by default uses columns 'guideRNA' and 'Group' to impose
#' additional restrictions on normalized events eg. only events created by the
#' same 'guideRNA' in the same 'Group' will be normalized.
#' @param min_freq (numeric) All events below this frequency are treated as
#' sequencing errors and rejected. This parameter is used during normalization
#' through \code{\link{amplicanNormalize}}.
#' @param continue (boolean) Default TRUE, decides whether to continue failed
#' ampliCan runs. In case of FALSE, all contents in `results` folder will
#' be removed.
#' @include amplicanAlign.R amplicanReport.R
#' @return (invisible) results_folder path
#' @export
#' @family analysis steps
#' @examples
#' # path to example config file
#' config <- system.file("extdata", "config.csv", package = "amplican")
#' # path to example fastq files
#' fastq_folder <- system.file("extdata", package = "amplican")
#' # output folder
#' results_folder <- tempdir()
#'
#' #full analysis, not knitting files automatically
#' amplicanPipeline(config, fastq_folder, results_folder, knit_reports = FALSE)
#'
# config <- system.file("extdata", "config.csv", package = "amplican")
# fastq_folder <- system.file("extdata", package = "amplican")
# results_folder <- tempdir()
# knit_reports = TRUE
# write_alignments_format = "txt"
# average_quality = 30
# min_quality = 0
# use_parallel = FALSE
# scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
#   match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA")
# gap_opening = 25
# gap_extension = 0
# fastqfiles = 0.5
# PRIMER_DIMER = 30
# event_filter = TRUE
# cut_buffer = 5
# primer_mismatch = 1
# promiscuous_consensus = TRUE
# normalize = c("guideRNA", "Group")
# donor_mismatch = 3
# min_freq = 0.01
amplicanPipeline <- amplicanPipe(0.01)


#' Wraps main package functionality into one function.
#'
#' amplicanPipelineIndexHopping is identical as amplicanPipeline except that
#' default \code{min_freq} threshold is set to 0.15. Setting this threshold
#' higher will decrease risks of inadequate normalization in cases of potential
#' Index Hopping, potentially decreasing precision of true editing rate calling.
#' Index Hopping can be mitigated with use  of unique dual indexing pooling
#' combinations. However, in cases when you might expect Index Hopping to occur
#' you should use this function instead of amplicanPipeline.
#'
#' \code{result_folder} and also knit prepared reports into 'reports' folder.
#' @inheritParams amplicanPipeline
#' @include amplicanAlign.R amplicanReport.R
#' @return (invisible) results_folder path
#' @export
#' @family analysis steps
#'
amplicanPipelineConservative <- amplicanPipe(0.15)
