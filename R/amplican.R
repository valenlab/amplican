#' amplican: fast and precise analysis of CRISPR experiments
#'
#' It has three main goals:
#'
#' \itemize{
#' \item Alignment and analysis of the MiSeq or HiSeq data.
#' \item Prepare automatic reports as .Rmd files that are flexible and open for manipulation.
#' \item Provide specialized plots for deletions, insertions, mismatches.
#' }
#'
#' To learn more about amplican, start with the vignettes:
#' \code{browseVignettes(package = "amplican")}
#'
#'
#' @docType package
#' @name amplican
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("amplican", libpath)
}


#' Wraps main package functionality into one function.
#'
#' amplican_pipeline is convinient wrapper around all most popular settings of
#' \code{\link{amplicanAnalysis}} and \code{\link{amplicanReport}}. It will generate all
#' results in the \code{result_folder} and also print prepared reports into "reports" folder.
#'
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}
#' @param fastq_folder (string) Path to FASTQ files. If not specified, FASTQ files should be
#' in the same directory as config file.
#' @param results_folder (string) Where do you want your results to be stored. The
#' program will create files in that folder so make sure you have writing permissions.
#' @param total_processors (numeric) Set this to the number of processors you want to use.
#' Default is 1. Works only if you have "doParallel" installed and accessible.
#' @param min_quality (numeric)  Similar as in average_quality, but this is the minimum quality for
#' ALL nucleotides. If one of them has quality BELLOW this threshold, then the sequence is skipped.
#' Default is 0.
#' @param fastqfiles (numeric) Normally you want to use both FASTQ files. But in
#'                         some special cases, you may want to use only the
#'                         forward file, or only the reverse file.
#'                         0 - Use both FASTQ files
#'                         1 - Use only the forward FASTQ file
#'                         2 - Use only the reverse FASTQ file
#' @param PRIMER_DIMER (numeric) Value specyfying buffer for PRIMER DIMER detection. For a given read it will be
#' recognized as PRIMER DIMER when alignment will introduce gap of size bigger than:
#' length of amplicon - (lenghts of PRIMERS + PRIMER_DIMER value)
#' @param cut_buffer (numeric) Value specyfying a buffer for PAM, this will add from both sides to
#' a window defined from uppercase letters in the amplicon. Deletions overlapping this window will be considered a
#' valid cut (if confirmed by both forward and rewerse reads).
#' @return NULL
#' @export
#'
amplicanPipeline <- function(config, fastq_folder, results_folder,
                              total_processors = 1,
                              min_quality = 0,
                              fastqfiles = 2,
                              PRIMER_DIMER = 30,
                              cut_buffer = 5) {

  amplicanAnalysis(config, fastq_folder, results_folder,
                   total_processors = total_processors,
                   min_quality = min_quality,
                   fastqfiles = fastqfiles,
                   PRIMER_DIMER = PRIMER_DIMER,
                   cut_buffer = cut_buffer)

  reportsFolder <- paste0(results_folder, "/reports")
  if (!dir.exists(reportsFolder)) {
    dir.create(file.path(reportsFolder))
  }

  amplicanReport(results_folder, report_files = paste0(reportsFolder, c("/report_id",
                                                                        "/report_barcode",
                                                                        "/report_group",
                                                                        "/report_guide",
                                                                        "/report_amplicon",
                                                                        "/report_summary")))
}
