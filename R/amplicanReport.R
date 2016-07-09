#' Prepare report as .Rmd file
#'
#' amplicanReport takes a configuration file, fastq reads and output directory to prepare
#' summary as editable .Rmd file. You can specify whether you want to make summaries based on
#' ID, Barcode, Experiment Type or even guideRNA.
#' @param results_folder (string) Folder containing results from the \code{\link{amplicanAnalysis}} function,
#' do not change names of the files.
#' @param levels (vector) Possible values are: "id", "barcode", "experiment_type", "guide". You can also input more than one
#' value eg. c("id", "barcode") will create two separate reports for each level.
#' @param report_files (vector) You can supply your own names of the files. For each of the levels there has to be one
#' file name. Files are created in current working directory by default.
#' @export
#'
amplicanReport <- function(results_folder,
                           levels = c("id", "barcode", "experiment_type", "guide"),
                           report_files = c("report_id", "report_barcode", "report_experiment", "report_guide")){

  existing_levels <- c("id", "barcode", "experiment_type", "guide")
  invalid_levels <- !levels %in% existing_levels
  if (any(invalid_levels)) {
    stop(paste("Invalid levels:", levels[invalid_levels]))
  }
  if (length(levels) != length(report_files)) {
    stop("report_files must provide name for each of the levels")
  }

  for (i in 1:length(levels)) {
    report_name <- paste0(report_files[i], ".Rmd")
    isRmdReady <- file.create(report_name, showWarnings = T)
    if (isRmdReady) {

      fileConn <- file(report_name)
      rmdContent <- switch(levels[i],
                           "id" = make_id_rmd(results_folder),
                           "barcode" = make_barcode_rmd(results_folder),
                           "experiment_type" = make_experiment_rmd(results_folder),
                           "guide" = make_guide_rmd(results_folder))
      writeLines(rmdContent, fileConn)
      close(fileConn)

    } else {
      stop(isRmdReady)
    }
  }
}
