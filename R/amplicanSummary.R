#' Prepare summary report as .Rmd file.
#'
#' amplicanSummary takes alignments_folder as directory of results to prepare
#' summary of all filters during read grouping as editable .Rmd file.
#' @param alignments_folder (string) Folder containing results from the \code{\link{amplicanAnalysis}} function,
#' do not change names of the files.
#' @param report_name (string) Name of the summary report.
#' @return NULL
#' @export
#'
amplicanSummary <- function(alignments_folder, report_name = "summary_report"){
  report_name <- paste0(report_name, ".Rmd")

  isRmdReady <- file.create(report_name, showWarnings = T)
  if (isRmdReady) {

    fileConn <- file(report_name)
    writeLines(c(write_head("Summary Read Report"),
                 "# Explanation of variables\n",
                 "***\n",
                 "**read_count** - how many reads belongs to this barcode  ",
                 "**bad_base_quality** - how many reads had base quality worse than specified (default is 0)  ",
                 "**bad_average_quality** - how many reads had average base quality worse than specified (default is 0)  ",
                 "**bad_alphabet** - how many reads had alphabet with bases other than A, C, G, T  ",
                 "**filtered_read_count** - how many reads were left after filtering  ",
                 "**unique_reads** - how many reads (forward and reverse together) for this barcode is unique  ",
                 "**unassigned_reads/assigned_reads** - how many reads have been not assigned/assigned to any of the",
                 "experiments\n",
                 "***\n",
                 "# Summary Table\n",
                 "***\n",
                 "```{r echo = F}",
                 "library(knitr)",
                 paste0("summaryDF <- read.csv('", alignments_folder, "/barcode_reads_filters.csv')  "),
                 "kable(summaryDF)",
                 "```\n",
                 "Table 1. Reads distributed for each barcode\n",
                 "***\n",
                 "# Barplot\n",
                 "***\n",
                 "```{r fig.width=8, fig.height=12, echo = F}",
                 "library(ggplot2)",
                 "library(reshape2)",
                 "summaryDFmelt = melt(summaryDF, id.vars = c(\"barcode\"), measure.vars = c(\"bad_base_quality\",",
                 "                                                                           \"bad_average_quality\",",
                 "                                                                           \"bad_alphabet\",",
                 "                                                                           \"unassigned_reads\",",
                 "                                                                           \"assigned_reads\"))",
                 "ggplot(data = summaryDFmelt, aes(x = barcode, y = value, fill = variable)) +",
                 "  geom_bar(position=\"stack\", stat=\"identity\") +",
                 "  ylab(\"number of reads\") +",
                 "  theme(legend.position = \"top\",",
                 "        legend.direction = \"horizontal\",",
                 "        legend.title = element_blank()) +",
                 "  coord_flip()",
                 "```\n",
                 "Plot 1. Reads distribution for each barcode.\n",
                 "*** \n"), fileConn)
    close(fileConn)

  } else {
    stop(isRmdReady)
  }
}
