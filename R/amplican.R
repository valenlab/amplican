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
#' @param knit_files (boolean) whether function should "knit" all reports
#' automatically for you (it is time consuming, but shows progress)
#' @inheritParams amplicanAlign
#' @importFrom rmarkdown render
#' @import knitr
#' @include amplicanAlign.R amplicanReport.R
#' @return (string) results_folder path
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
#' amplicanPipeline(config, fastq_folder, results_folder, knit_files = FALSE)
#'
amplicanPipeline <- function(config,
                             fastq_folder,
                             results_folder,
                             knit_files = TRUE,
                             total_processors = 1,
                             min_quality = 0,
                             fastqfiles = 0,
                             PRIMER_DIMER = 30,
                             cut_buffer = 5) {

    amplicanAlign(config,
                  fastq_folder,
                  results_folder,
                  total_processors = total_processors,
                  min_quality = min_quality,
                  fastqfiles = fastqfiles,
                  PRIMER_DIMER = PRIMER_DIMER,
                  cut_buffer = cut_buffer)

    reportsFolder <- file.path(results_folder, "reports")
    if (!dir.exists(reportsFolder)) {
      dir.create(reportsFolder)
    }

    amplicanReport(results_folder,
                   report_files = file.path(reportsFolder,
                                            c("report_id",
                                              "report_barcode",
                                              "report_group",
                                              "report_guide",
                                              "report_amplicon",
                                              "report_summary")))

    if (knit_files) {
      rmarkdown::render(file.path(reportsFolder, "report_id.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
      rmarkdown::render(file.path(reportsFolder, "report_barcode.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
      rmarkdown::render(file.path(reportsFolder, "report_group.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
      rmarkdown::render(file.path(reportsFolder, "report_guide.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
      rmarkdown::render(file.path(reportsFolder, "report_amplicon.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
      rmarkdown::render(file.path(reportsFolder, "report_summary.Rmd"),
                        output_dir = file.path(reportsFolder, "html"))
    }

    invisible(reportsFolder)
}
