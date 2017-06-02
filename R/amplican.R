#' amplican: fast and precise analysis of CRISPR experiments
#'
#' It has three main goals:
#'
#' \enumerate{
#' \item Alignment and analysis of the MiSeq or HiSeq data.
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
NULL


#' Wraps main package functionality into one function.
#'
#' amplicanPipeline is convenient wrapper around all most popular settings of
#' \code{\link{amplicanAlign}} and \code{\link{amplicanReport}}.
#' It will generate all results in the \code{result_folder} and also print
#' prepared reports into 'reports' folder.
#'
#' @param make_reports (boolean) whether function should create and "knit" all
#' reports automatically for you (it is time consuming, but shows progress)
#' @inheritParams amplicanAlign
#' @importFrom rmarkdown render
#' @import knitr
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
#' amplicanPipeline(config, fastq_folder, results_folder, make_reports = FALSE)
#'
amplicanPipeline <- function(config,
                             fastq_folder,
                             make_reports = TRUE,
                             results_folder,
                             average_quality = 30,
                             min_quality = 20,
                             total_processors = 1,
                             fastqfiles = 0,
                             PRIMER_DIMER = 30,
                             cut_buffer = 5) {

    amplicanAlign(config,
                  fastq_folder,
                  results_folder,
                  average_quality = average_quality,
                  min_quality = min_quality,
                  total_processors = total_processors,
                  fastqfiles = fastqfiles,
                  PRIMER_DIMER = PRIMER_DIMER,
                  cut_buffer = cut_buffer)

    if (make_reports) {

      reportsFolder <- file.path(results_folder, "reports")
      if (!dir.exists(reportsFolder)) {
        dir.create(reportsFolder)
      }

      message(paste0("Making reports. Due to high quality ",
                     "figures it is time consuming. Use .Rmd templates for ",
                     "more control."))
      amplicanReport(results_folder,
                     report_files = file.path(reportsFolder,
                                              c("report_id",
                                                "report_barcode",
                                                "report_group",
                                                "report_guide",
                                                "report_amplicon",
                                                "index")))
    }

    invisible(results_folder)
}
