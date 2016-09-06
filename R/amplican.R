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
#' \code{browseVignettes(package = 'amplican')}
#'
#'
#' @docType package
#' @name amplican
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("amplican", libpath)
}


#' Wraps main package functionality into one function.
#'
#' amplican_pipeline is convinient wrapper around all most popular settings of
#' \code{\link{amplicanAlign}} and \code{\link{amplicanReport}}.
#' It will generate all results in the \code{result_folder} and also print
#' prepared reports into 'reports' folder.
#'
#' @param knit_files (bolean) whether we should "knit" all reports automatically for you
#' (it is time consuming)
#' @inheritParams amplicanAlign
#' @importFrom rmarkdown render
#' @import knitr
#' @return NULL All results are created in specified results_folder.
#' @export
#' @family analysis steps
#' @examples
#' config <- system.file("extdata", "config.csv", package = "amplican") #example config file
#' fastq_folder <- system.file("extdata", "", package = "amplican") #path to example fastq files
#' results_folder <- paste0(fastq_folder, "/results") #output folder
#'
#' #full analysis, not kniting files automatically
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

    reportsFolder <- paste0(results_folder, "/reports")
    if (!dir.exists(reportsFolder)) {
        dir.create(file.path(reportsFolder))
    }

    amplicanReport(results_folder, report_files = paste0(reportsFolder,
                                                         c("/report_id",
                                                           "/report_barcode",
                                                           "/report_group",
                                                           "/report_guide",
                                                           "/report_amplicon",
                                                           "/report_summary")))

    if (knit_files) {
        rmarkdown::render(paste0(reportsFolder, "/report_id.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
        rmarkdown::render(paste0(reportsFolder, "/report_barcode.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
        rmarkdown::render(paste0(reportsFolder, "/report_group.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
        rmarkdown::render(paste0(reportsFolder, "/report_guide.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
        rmarkdown::render(paste0(reportsFolder, "/report_amplicon.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
        rmarkdown::render(paste0(reportsFolder, "/report_summary.Rmd"),
                          output_dir = paste0(reportsFolder, "/html"))
    }
}
