#' Prepare reports as .Rmd files.
#'
#' amplicanReport takes a configuration file, fastq reads and output directory
#' to prepare summaries as an editable .Rmd file. You can specify whether you
#' want to make summaries based on ID, Barcode, Group or even guideRNA and
#' Amplicon. This function autmatically knits all reports after creation.
#' If you want to postpone kniting and edit reports, use .Rmd templates to
#' create your own version of reports instead of this function.
#' @param results_folder (string) Folder containing results from the
#' \code{\link{amplicanAlign}} function, do not change names of the files.
#' @param levels (vector) Possible values are: "id", "barcode", "group",
#' "guide", "amplicon", "summary". You can also input more than one value
#' eg. c("id", "barcode") will create two separate reports for each level.
#' @param report_files (vector) You can supply your own names of the files.
#' For each of the levels there has to be one file name. Files are created
#' in current working directory by default.
#' @param cut_buffer (numeric) Default 5. A number of bases that is used around
#' the specified cut site.
#' @param xlab_spacing (numeric) Default is 4. Spacing of the ticks on the x
#' axis of plots.
#' @param top (numeric) Default is 5. How many of the top most frequent
#' unassigned reads to report? It is only relevant when you used forward and
#' reverse reads. We align them to each other as we could not specify correct
#' amplicon.
#' @param knit_reports (boolean) Whether to knit reports automatically.
#' @include helpers_rmd.R
#' @return (string) Path to the folder with results.
#' @export
#' @family analysis steps
#' @examples
#' results_folder <- tempdir()
#' amplicanReport(results_folder, report_files = file.path(results_folder,
#'                                                         c("id_report",
#'                                                           "barcode_report",
#'                                                           "group_report",
#'                                                           "guide_report",
#'                                                           "amplicon_report",
#'                                                           "index")),
#'                knit_reports = FALSE)
amplicanReport <- function(results_folder,
                           levels = c("id",
                                      "barcode",
                                      "group",
                                      "guide",
                                      "amplicon",
                                      "summary"),
                           report_files = c("id_report",
                                            "barcode_report",
                                            "group_report",
                                            "guide_report",
                                            "amplicon_report",
                                            "index"),
                           cut_buffer = 5,
                           xlab_spacing = 4,
                           top = 5,
                           knit_reports = TRUE) {

  existing_levels <- c("id", "barcode", "group",
                       "guide", "amplicon", "summary")
  template_names <- c("id_report", "barcode_report", "group_report",
                      "guide_report", "amplicon_report", "index_report")
  invalid_levels <- !levels %in% existing_levels

  if (any(invalid_levels)) {
    stop(paste("Invalid levels:", levels[invalid_levels]))
  }
  if (length(levels) != length(report_files)) {
    stop("report_files must provide name for each of the levels")
  }
  template_names <- template_names[!invalid_levels]
  template_names <- template_names[match(levels, existing_levels)]

  links <- c("***\n", "# Other Reports\n", "***\n")
  if ("summary" %in% levels) {
    summary <- which(levels == "summary")
    lvl_no_sum <- levels[-summary]
    files_no_sum <- basename(report_files[-summary])
  }
  for (i in seq_along(lvl_no_sum)) {
    links <- c(links, paste0(i, ". [Report by ",
                             lvl_no_sum[i], "](",
                             file.path(".", files_no_sum[i]), ".html)"))
  }

  barcode_summary <- file.path(results_folder, 'barcode_reads_filters.csv')
  config_summary <- file.path(results_folder, 'config_summary.csv')
  alignments <- file.path(results_folder, 'alignments',
                          'events_filtered_shifted_normalized.csv')
  unassigned_folder<- file.path(results_folder,
                                'alignments',
                                'unassigned_sequences')

  id_am_params <- list(alignments = alignments,
                       config_summary = config_summary,
                       cut_buffer = cut_buffer,
                       xlab_spacing = xlab_spacing)
  gr_ga_params <- list(alignments = alignments,
                       config_summary = config_summary)

  for (i in seq_along(levels)) {

    report_name <- paste0(report_files[i], ".Rmd")
    rmarkdown::draft(file = report_name,
                     template = template_names[i],
                     package = "amplican",
                     edit = FALSE)

    rmdParamList <- switch(levels[i],
                           id = id_am_params,
                           barcode = list(alignments = alignments,
                                          config_summary = config_summary,
                                          unassigned_folder = unassigned_folder,
                                          top = top),
                           group = gr_ga_params,
                           guide = gr_ga_params,
                           amplicon = id_am_params,
                           summary = list(barcode_summary = barcode_summary,
                                          config_summary = config_summary,
                                          links = links))
    rmd_content <- readLines(report_name)
    for (k in seq_along(rmdParamList)) {
      # 11th line is params:
      new_param <- if (is.numeric(rmdParamList[[k]])) {
        paste0(": ", rmdParamList[[k]], collapse = "")
      } else {
        paste0(": '", rmdParamList[[k]], "'", collapse = "")
      }
      rmd_content[11 + k] <- gsub(":.*", new_param, rmd_content[11 + k])
    }
    cat(rmd_content, file = report_name, sep = "\n")

    if (knit_reports) {
      rmarkdown::render(report_name,
                        params = rmdParamList)
    }
  }

  invisible(results_folder)
}
