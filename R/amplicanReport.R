#' Prepare report as .Rmd file
#'
#' amplicanReport takes a configuration file, fastq reads and output directory to prepare
#' summary as editable .Rmd file. You can specify whether you want to make summaries based on
#' ID, Barcode or Experiment Type.
#' @param alignments_folder (string) Folder containing results from the \code{\link{amplicanAnalysis}} function,
#' do not change names of the files.
#' @param levels (vector) Possible values are: "ID", "Barcode", "Experiment_Type". You can also input more than one
#' value eg. c("ID", "Barcode") will create two separate reports for each level.
#' @param report_files (vector) You can supply your own names of the files. For each of the levels there has to be one
#' file name. Files are created in current working directory by default.
#' @export
#'
amplicanReport <- function(alignments_folder, levels = c("ID"), report_files = c("report_id")){
  #check whether alignments folder have alignments
  #check config file whether it is correct
  #check whether aggregate has values of "ID", "Barcode", "Experiment_Type"
  #for each of levels
  ##make .Rmd file with plots and stuff

  # message("Checking write access...")
  # checkFileWriteAccess(results_folder)
  # if (temp_folder != "") {checkFileWriteAccess(temp_folder)}
  #
  # message("Checking configuration file...")
  # configTable <- read.table(config, header=FALSE, sep="\t", strip.white=TRUE)
  # colnames(configTable) <- c("ID", "Barcode", "Forward_Reads_File", "Reverse_Reads_File", "Experiment_Type",
  #                            "Target_Primer", "Forward_Primer", "Reverse_Primer", "Strand", "Amplicon")
  # configTable$Forward_Reads_File <- paste(fastq_folder, configTable$Forward_Reads_File, sep = "/")
  # configTable$Reverse_Reads_File <- paste(fastq_folder, configTable$Reverse_Reads_File, sep = "/")
  # checkConfigFile(configTable, fastq_folder)
}
