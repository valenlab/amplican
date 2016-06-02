#' Prepare alignments.
#'
#' ampliCanAnalysis takes a configuration files, fastq reads and output directory to prepare
#' alignments and summary.
#' @param config (string) The path to your configuration file. For example:
#'                      /Home/johndoe/.../AmpliCan/res/Cas9_toy/run11.txt
#' @param fastq_folder (string) Path to FASTQ files. If not specified, FASTQ files should be
#' in the same directory as config file.
#' @param results_folder (string) Where do you want your results to be stored. The
#' program will create files in that folder so make sure you have writing permissions.
#' @param total_processors (int) Set this to the number of processors you want to use.
#' Default is 1. Works only if you have "doParallel" installed and accessible.
#' @param skip_bad_nucleotides (logical) Some sequences have faulty nucleotides labels
#' with N. If we find a sequence like that in either forwards or reverse, we skip that
#' aligment. Default is TRUE.
#' @param average_quality (int) The FASTQ file have a quality for each nucleotide,
#' being ! the lower and ~ the highest. In ASCII :
#'                              !"#$%&'()*+,-./0123456789:;<=>?@
#'                              ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
#'                              abcdefghijklmnopqrstuvwxyz{|}~
#'
#' This quality variable goes from 0 to 100, being 0 the lowest quality and 100 the highest.
#' You can write whatever number in between, and the program will find out the apropiate character
#' encoding. The filter works by converting each character to a number, and then finding the average.
#' If the average fall above this threshold then we take the sequence. Default is 0.
#' @param min_quality (int)  Similar as in average_quality, but this is the minimum quality for
#' ALL nucleotides. If one of them has quality BELLOW this threshold, then the sequence is skipped.
#' Default is 0.
#' @param write_alignments (int) How to write the aligments results into disk:
#'                         0 - Write nothing
#'                         1 - Write only the summary file
#'                         2 - Write also a verbose file with all the alignments
#'                             in the same .txt file (Default option)
#' @param scoring_matrix (string) For now the only option is "NUC44".
#' @param gap_opening (int) The opening gap score. Default is 50.
#' @param gap_extension (int) The gap extension score. Default is 0.
#' @param gap_ending (logical) If you want that the ending gap count for the
#' alignment score set this to TRUE. Default is FALSE.
#' @param far_indels (logical) If the ending/starting gap should be considered to be
#' an indel. Default is TRUE. If you want to filter these out from the plots,
#' there is an option to do so in the plot function. You don't need to do it here.
#' @param deletefq (logical) If you have fastq.gz files they will be uncompressed into
#' the same directory. Set this to true if you want to delete uncompressed files afterwards.
#' Default is FALSE.
#' @param temp_folder (string) Your FASTQ files can be compressed in a file. If this
#'                         happens, the program needs to uncompress them first.
#'                         In order to do so, we need a folder where they will
#'                         be placed. In here, you can specify the path where
#'                         you want this files to be placed. If you don't
#'                         specify a path, they will be placed in the same
#'                         folder where the original files are. In any case you
#'                         will need writing permissions in the folder in order
#'                         to uncompress everything.
#' @param fastqfiles (int) Normally you want to use both FASTQ files. But in
#'                         some special cases, you may want to use only the
#'                         forward file, or only the reverse file.
#'                         0 - Use both FASTQ files
#'                         1 - Use only the forward FASTQ file
#'                         2 - Use only the reverse FASTQ file
#' @include gotoh.R helpers_alignment.R helpers_filters.R helpers_warnings.R helpers_directory.R
#' @export
#'
amplicanAnalysis <- function(config,
                             fastq_folder,
                             results_folder,
                             total_processors = 1,
                             skip_bad_nucleotides = TRUE,
                             average_quality = 0,
                             min_quality = 0,
                             write_alignments = 2,
                             scoring_matrix = "NUC44",
                             gap_opening = 50,
                             gap_extension = 0,
                             gap_ending = FALSE,
                             far_indels = TRUE,
                             deletefq = FALSE,
                             temp_folder = "",
                             fastqfiles = 0){
  message("Checking write access...")
  checkFileWriteAccess(results_folder)
  if (temp_folder != "") {checkFileWriteAccess(temp_folder)}

  message("Checking configuration file...")
  configTable <- read.table(config, header=FALSE, sep="\t", strip.white=TRUE)
  colnames(configTable) <- c("ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Experiment_Type",
                             "Target_Primer","Forward_Primer","Reverse_Primer","Strand","Amplicon")
  configTable$Forward_Reads_File <- paste(fastq_folder, configTable$Forward_Reads_File, sep = "/")
  configTable$Reverse_Reads_File <- paste(fastq_folder, configTable$Reverse_Reads_File, sep = "/")
  checkConfigFile(configTable, fastq_folder)

  message("Preparing FASTQ files...")
  configTable <- unpackFastq(configTable, temp_folder)

  resultsFolder <- paste0(results_folder, "/alignments/")
  if(!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  unassignedFolder <- paste0(resultsFolder, "/unassigned_sequences")
  if (!dir.exists(unassignedFolder)) {
    dir.create(file.path(unassignedFolder))
  }

  # Parameters
  logFileName <- paste(results_folder, "/RunParameters.txt", sep = '')
  if (file.exists(logFileName)) {file.remove(logFileName)}
  logFileConn <- file(logFileName, open="at")
  writeLines(paste("Config file:           ", config), logFileConn)
  writeLines(paste("Total Processors:      ", total_processors), logFileConn)
  writeLines(paste("Skip Bad Nucleotides:  ", skip_bad_nucleotides), logFileConn)
  writeLines(paste("Average Quality:       ", average_quality), logFileConn)
  writeLines(paste("Minimum Quality:       ", min_quality), logFileConn)
  writeLines(paste("Write Alignments Mode: ", write_alignments), logFileConn)
  writeLines(paste("Scoring Matrix:        ", scoring_matrix), logFileConn)
  writeLines(paste("Gap Opening:           ", gap_opening), logFileConn)
  writeLines(paste("Gap Extension:         ", gap_extension), logFileConn)
  writeLines(paste("Gap Ending:            ", gap_ending), logFileConn)
  writeLines(paste("Far Indels:            ", far_indels), logFileConn)
  close(logFileConn)

  uBarcode <- unique(configTable$Barcode)

  # Several statistics about deletions, cuts and reads
  configTable$Sum_Target_Forward_Found <- 0
  configTable$Sum_Target_Reverse_Found <- 0
  configTable$Sum_Deletions_Forward <- 0
  configTable$Sum_Deletions_Reverse <- 0
  configTable$Sum_Insertions_Forward <- 0
  configTable$Sum_Insertions_Reverse <- 0
  configTable$Sum_Missmatches_Forward <- 0
  configTable$Sum_Missmatches_Reverse <- 0

  # Warnings
  configTable$Found_Target <- 0
  configTable$Found_Primers <- 0
  configTable$Found_AP <- 0

  if (requireNamespace("doParallel", quietly = TRUE) & total_processors > 1) {
    cl <- parallel::makeCluster(total_processors, outfile="")
    doParallel::registerDoParallel(cl)

    foreach::foreach(j=1:length(uBarcode), .export=c('getEventInfo', 'upperGroups',
                                                     'checkTarget', 'checkPrimers', 'checkPositions',
                                                     'goodBaseQuality', 'goodAvgQuality', 'alphabetQuality',
                                                     'gRCPP'),
                     .packages = c('Rcpp', 'R.utils', 'GenomicRanges', 'ShortRead', 'seqinr')) %dopar% {
      makeAlignment(configTable[configTable$Barcode == uBarcode[j],],
                    resultsFolder,
                    skip_bad_nucleotides,
                    average_quality,
                    min_quality,
                    write_alignments,
                    scoring_matrix,
                    gap_opening,
                    gap_extension,
                    gap_ending,
                    far_indels)
    }
    parallel::stopCluster(cl)
  } else {
      for(j in 1:length(uBarcode)){
        makeAlignment(configTable[configTable$Barcode == uBarcode[j],],
                      resultsFolder,
                      skip_bad_nucleotides,
                      average_quality,
                      min_quality,
                      write_alignments,
                      scoring_matrix,
                      gap_opening,
                      gap_extension,
                      gap_ending,
                      far_indels)
      }
  }

  message("Alignments done.")

  # Put all the logs and all the configs together
  totalLogs <- unifyFiles(resultsFolder, "SUBLOG", paste0(results_folder, "/alignmentLog.txt"), header = F)
  totalConfigs <- unifyFiles(resultsFolder, "configFile_results", paste0(results_folder, "/config_summary.csv"))
  totalBarcode <- unifyFiles(resultsFolder, "reads_filters.csv", paste0(results_folder, "/barcode_reads_filters.csv"))

  # If the user want to delete the uncompressed results, do it now.
  if (deletefq == TRUE) {
    message("Deleting temporary files...")
    deleteFiles(configTable)
  }
  message("Finished.")
}
