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
#' amplican_default is convinient wrapper around all default settings of
#' amplicanAnalysis, amplicanReport and amplicanSummary. It will generate all
#' results in the result_folder and also print prepared reports into reports folder.
#'
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}
#' @param fastq_folder (string) Path to FASTQ files. If not specified, FASTQ files should be
#' in the same directory as config file.
#' @param results_folder (string) Where do you want your results to be stored. The
#' program will create files in that folder so make sure you have writing permissions.
#' @return NULL
#' @export
#'
amplican_default <- function(config, fastq_folder, results_folder) {
  amplicanAnalysis(config, fastq_folder, results_folder)

  reportsFolder <- paste0(results_folder, "/reports")
  if (!dir.exists(reportsFolder)) {
    dir.create(file.path(reportsFolder))
  }

  amplicanSummary(results_folder, report_name = paste0(reportsFolder, "/summary_report"))
  amplicanReport(results_folder, report_files = paste0(reportsFolder, c("/report_id",
                                                                        "/report_barcode",
                                                                        "/report_group",
                                                                        "/report_guide",
                                                                        "/report_amplicon")))
}
