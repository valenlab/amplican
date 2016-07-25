#' Prepare reports as .Rmd files.
#'
#' amplicanReport takes a configuration file, fastq reads and output directory
#' to prepare summary as editable .Rmd file. You can specify whether you want
#' to make summaries based on ID, Barcode, Group or even guideRNA and Amplicon.
#' @param results_folder (string) Folder containing results from the
#' \code{\link{amplicanAnalysis}} function, do not change names of the files.
#' @param levels (vector) Possible values are: 'id', 'barcode', 'group',
#' 'guide', 'amplicon', 'summary'. You can also input more than one value
#' eg. c('id', 'barcode') will create two separate reports for each level.
#' @param report_files (vector) You can supply your own names of the files.
#' For each of the levels there has to be one file name. Files are created
#' in current working directory by default.
#' @return NULL All results are written into current dir with specified names.
#' @export
#' @examples
#' results_folder <- system.file("extdata", "results", package = "amplican") #output folder
#' amplicanReport(results_folder, report_files = paste0(paste0(results_folder, "/reports"),
#'                                                      c("/report_id",
#'                                                        "/report_barcode",
#'                                                        "/report_group",
#'                                                        "/report_guide",
#'                                                        "/report_amplicon",
#'                                                        "/report_summary")))
#'
amplicanReport <- function(results_folder,
                           levels = c("id",
                                      "barcode",
                                      "group",
                                      "guide",
                                      "amplicon",
                                      "summary"),
                           report_files = c("report_id",
                                            "report_barcode",
                                            "report_group",
                                            "report_guide",
                                            "report_amplicon",
                                            "report_summary")) {

    existing_levels <- c("id", "barcode", "group",
                         "guide", "amplicon", "summary")
    invalid_levels <- !levels %in% existing_levels
    if (any(invalid_levels)) {
        stop(paste("Invalid levels:", levels[invalid_levels]))
    }
    if (length(levels) != length(report_files)) {
        stop("report_files must provide name for each of the levels")
    }

    for (i in 1:length(levels)) {
        report_name <- paste0(report_files[i], ".Rmd")
        isRmdReady <- file.create(report_name, showWarnings = TRUE)
        if (isRmdReady) {

            fileConn <- file(report_name)
            rmdContent <- switch(levels[i],
                                 id = make_id_rmd(results_folder),
                                 barcode = make_barcode_rmd(results_folder),
                                 group = make_group_rmd(results_folder),
                                 guide = make_guide_rmd(results_folder),
                                 amplicon = make_amplicon_rmd(results_folder),
                                 summary = make_summary_rmd(results_folder))
            writeLines(rmdContent, fileConn)
            close(fileConn)

        } else {
            stop(isRmdReady)
        }
    }
}
