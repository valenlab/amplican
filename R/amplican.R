#' amplican: fast and precise analysis of CRISPR experiments
#'
#' It has three main goals:
#'
#' \enumerate{
#' \item Alignment and analysis of the MiSeq or HiSeq data.
#'
#' \item Prepare automatic reports as .Rmd files that are flexible
#' and open for manipulation.
#' \item Provide specialized plots for deletions, insertions, mismatches.
#' }
#'
#' To learn more about amplican, start with the vignettes:
#' \code{browseVignettes(package = "amplican")}
#'
#' @docType package
#' @name amplican
#' @import methods BiocGenerics Biostrings data.table
"_PACKAGE"


#' Wraps main package functionality into one function.
#'
#' amplicanPipeline is convenient wrapper around all functionality of the
#' package with the most robust settings. It will generate all results in the
#' \code{result_folder} and also print prepared reports into 'reports' folder.
#' @param results_folder (string) Where do you want your results to be stored?
#' The package will create files in that folder so make sure you have writing
#' permissions.
#' @param knit_reports (boolean) whether function should "knit" all
#' reports automatically for you (it is time consuming, be patient), when false
#' reports will be prepared, but not knitted
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}
#' @param fastq_folder (string) Path to FASTQ files. If not specified,
#' FASTQ files should be in the same directory as config file.
#' @param total_processors (numeric) Set this to the number of processors you
#' want to use. Default is 1. Works only if you have
#' \code{\link[doParallel]{doParallel}} installed and accessible.
#' @param average_quality (numeric) The FASTQ file have a quality for each
#' nucleotide, depending on sequencing technology there exist many formats.
#' This package uses \code{\link[ShortRead]{readFastq}} to parse the reads.
#' If the average quality of the reads fall below value of
#' \code{average_quality} then sequence is filtered. Default is 0.
#' @param min_quality (numeric)  Similar as in average_quality, but this is
#' the minimum quality for ALL nucleotides in given read. If one of nucleotides
#' has quality BELLOW \code{min_quality}, then the sequence is filtered.
#' Default is 20.
#' @param write_alignments_format (character vector) Whether we should write
#' alignments results to separate files for each ID, reverse reads are reverse
#' complemented to minimize alignment bias, possible options are:
#' \itemize{
#'  \item{"fasta"}{ outputs alignments in fasta format where header indicates
#' experiment ID, read id and number of reads}
#'  \item{"txt"}{ simple format, read information followed by foward read and
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
#' above formats, pass a vector to get alginments in multiple formats.}
#' }
#' @param scoring_matrix (matrix) Default is 'NUC44'. Pass desired matrix using
#' \code{\link[Biostrings]{nucleotideSubstitutionMatrix}}
#' @param gap_opening (numeric) The opening gap score. Default is 50.
#' @param gap_extension (numeric) The gap extension score. Default is 30.
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
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @param cut_buffer The number of bases by which extend expected cut sites
#' (specified as UPPER case letters in the amplicon) in 5' and 3' directions.
#' @param normalize (character vector) If column 'Control' in config table
#' has all FALSE/0 values then normalization is skipped. Otherwise,
#' normalization is strict, which means events that are
#' found in 'Control' TRUE group will be removed in 'Control' FALSE group.
#' This parameter by default uses columns 'guideRNA' and 'Group' to impose
#' additional restrictions on normalized events eg. only events created by the
#' same 'guideRNA' will be normalized.
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
# total_processors = 1
# scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
#   match = 5, mismatch = -4, baseOnly = TRUE, type = "DNA")
# gap_opening = 50
# gap_extension = 0
# fastqfiles = 0
# PRIMER_DIMER = 30
# cut_buffer = 5
# normalize = c("guideRNA", "Group")

amplicanPipeline <- function(
  config, fastq_folder, results_folder, knit_reports = TRUE,
  write_alignments_format = "txt", average_quality = 30,
  min_quality = 0, total_processors = 1,
  scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
    match = 5, mismatch = -4, baseOnly = TRUE, type = "DNA"),
  gap_opening = 50, gap_extension = 0, fastqfiles = 0,
  PRIMER_DIMER = 30, cut_buffer = 5,
  normalize = c("guideRNA", "Group")) {

  message("Checking write access...")
  checkFileWriteAccess(results_folder)

  aln <- amplicanAlign(config = config,
                       fastq_folder = fastq_folder,
                       total_processors = total_processors,
                       average_quality = average_quality,
                       scoring_matrix = scoring_matrix,
                       gap_opening = gap_opening,
                       gap_extension = gap_extension,
                       min_quality = min_quality,
                       fastqfiles = fastqfiles)
  message("Saving alignments...")

  resultsFolder <- file.path(results_folder, "alignments")
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # save as .rds object
  saveRDS(aln, file.path(resultsFolder, "AlignmentsExperimentSet.rds"))
  # save as other formats
  if (!"None" %in% write_alignments_format) {
    for (frmt in write_alignments_format) {
      writeAlignments(aln, file.path(resultsFolder,
                                     paste0("alignments.", frmt)), frmt)
    }
  }
  message("Saving parameters...")
  logFileName <- file.path(results_folder, "RunParameters.txt")
  if (file.exists(logFileName)) {
    file.remove(logFileName)
  }
  logFileConn <- file(logFileName, open = "at")
  writeLines(c(paste("Config file:        ", config),
               paste("Processors used:    ", total_processors),
               paste("Average Quality:    ", average_quality),
               paste("Minimum Quality:    ", min_quality),
               paste("Write Alignments:   ", toString(write_alignments_format)),
               paste("Fastq files Mode:   ", fastqfiles),
               paste("Gap Opening:        ", gap_opening),
               paste("Gap Extension:      ", gap_extension),
               paste("PRIMER DIMER buffer:", PRIMER_DIMER),
               paste("Cut buffer:", cut_buffer),
               "Scoring Matrix:"), logFileConn)
  utils::write.csv(scoring_matrix, logFileConn, quote = FALSE, row.names = TRUE)
  close(logFileConn)

  message("Saving unassigned sequences...")
  data.table::fwrite(unassignedData(aln),
                     file.path(resultsFolder, "unassigned_reads.csv"))
  message("Saving barcode statistics...")
  data.table::fwrite(barcodeData(aln),
                     file.path(results_folder, "barcode_reads_filters.csv"))
  message("Translating alignments into events...")
  cfgT <- experimentData(aln)
  aln <- extractEvents(aln, total_processors = total_processors)
  message("Saving complete events - unfiltered...")
  data.table::fwrite(aln, file.path(resultsFolder, "raw_events.csv"))
  data.table::setDT(aln)
  seqnames <- read_id <- counts <- NULL
  # find PRIMER DIMERS
  PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

  # summarize how many PRIMER DIMER reads per ID
  onlyPD <- aln[PD, ]
  onlyPD <- unique(onlyPD, by = c("seqnames", "read_id"))
  summaryPD <- onlyPD[, list(counts  = sum(counts)), by = c("seqnames")]
  cfgT$PRIMER_DIMER <- 0
  cfgT$PRIMER_DIMER[match(summaryPD$seqnames, cfgT$ID)] <- summaryPD$counts
  cfgT$Reads_noPD <- cfgT$Reads - cfgT$PRIMER_DIMER

  # apply filter - remove all events that come from PD infected reads
  aln <- aln[!onlyPD, on = list(seqnames, read_id)]

  # filter events overlaping primers
  eOP <- findEOP(aln, cfgT)
  aln <- aln[!eOP, ]
  data.table::setDF(aln)

  # shift to relative (most left UPPER case is position 0)
  message("Shifting events as relative...")
  aln <- data.frame(map_to_relative(aln, cfgT), stringsAsFactors = FALSE)
  message("Saving shifted events - filtered...")
  data.table::fwrite(aln,
                     file.path(resultsFolder,
                               "events_filtered_shifted.csv"))
  # revert guides to 5'-3'
  cfgT$guideRNA[cfgT$Direction] <- revComp(cfgT$guideRNA[cfgT$Direction])
  # normalize
  message("Normalizing events...")
  aln <- amplicanNormalize(aln, cfgT, add = normalize)
  aln$overlaps <- amplicanOverlap(aln, cfgT, cut_buffer = cut_buffer)
  message("Saving normalized events...")
  data.table::fwrite(aln,
                     file.path(resultsFolder,
                               "events_filtered_shifted_normalized.csv"))
  # summarize
  cfgT <- amplicanSummarize(aln, cfgT)
  data.table::fwrite(
    cfgT[, c("ID", "Barcode", "Forward_Reads_File", "Reverse_Reads_File",
             "Group", "guideRNA", "Found_Guide", "Control", "Forward_Primer",
             "Reverse_Primer", "Direction", "Amplicon", "Reads", "PRIMER_DIMER",
             "Reads_noPD", "Reads_Cut", "Reads_Frameshifted",
             "Reads_Frameshifted_Overlapping")],
    file.path(results_folder, "config_summary.csv"))

  # reports
  reportsFolder <- file.path(results_folder, "reports")
  if (!dir.exists(reportsFolder)) {
    dir.create(reportsFolder)
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
