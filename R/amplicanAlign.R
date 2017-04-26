#' Align reads to amplicons.
#'
#' amplicanAlign takes a configuration files, fastq reads and output
#' directory to prepare alignments and summary. Finally, it returns a GRanges
#' object containing all mismatches, indels and insertions from our alignments.
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}
#' @param fastq_folder (string) Path to FASTQ files. If not specified,
#' FASTQ files should be in the same directory as config file.
#' @param results_folder (string) Where do you want your results to be stored.
#' The package will create files in that folder so make sure you have writing
#' permissions.
#' @param total_processors (numeric) Set this to the number of processors you
#' want to use. Default is 1. Works only if you have 'doParallel' installed
#' and accessible.
#' @param skip_bad_nucleotides (logical) Some sequences have faulty nucleotides
#' labels with N. If we find a sequence like that in either forwards or
#' reverse, we skip that alignment. Default is TRUE.
#' @param average_quality (numeric) The FASTQ file have a quality for each
#' nucleotide, being ! the lower and ~ the highest. In ASCII :
#'                              !'#$%&'()*+,-./0123456789:;<=>?@
#'                              ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
#'                              abcdefghijklmnopqrstuvwxyz{|}~
#'
#' This quality variable goes from 0 to 100, being 0 the lowest quality and
#' 100 the highest. You can write whatever number in between, and the program
#' will find out the appropriate character encoding. The filter works by
#' converting each character to a number, and then finding the average. If the
#' average fall above this threshold then we take the sequence. Default is 0.
#' @param min_quality (numeric)  Similar as in average_quality, but this is
#' the minimum quality for ALL nucleotides. If one of them has quality BELLOW
#' this threshold, then the sequence is skipped. Default is 20.
#' @param write_alignments (boolean) Whether we should write alignments results
#' to separate files for each ID
#' @param scoring_matrix (string) For now the only option is 'NUC44'.
#' @param gap_opening (numeric) The opening gap score. Default is 50.
#' @param gap_extension (numeric) The gap extension score. Default is 30.
#' @param fastqfiles (numeric) Normally you want to use both FASTQ files. But in
#'                         some special cases, you may want to use only the
#'                         forward file, or only the reverse file.
#'                         0 - Use both FASTQ files
#'                         0.5 - Use both FASTQ files,
#'                         but only for one of the reads (forward or reverse)
#'                         is required to have primer perfectly matched to
#'                         sequence - eg. use when reverse reads are
#'                         trimmed of primers, but forward reads have
#'                         forward primer in the sequence
#'                         1 - Use only the forward FASTQ file
#'                         2 - Use only the reverse FASTQ file
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than:
#' length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)
#' @param cut_buffer (numeric) Value specifying a buffer for PAM, this will add
#' from both sides to a window defined from uppercase letters in the amplicon.
#' Deletions overlapping this window will be considered a
#' valid cut (if confirmed by both forward and reverse reads).
#' @return (string) Path to results_folder.
#' @include helpers_alignment.R helpers_warnings.R helpers_directory.R
#' @import doParallel foreach GenomicRanges
#' @import BiocParallel Biostrings
#' @importFrom utils read.csv
#' @export
#' @family analysis steps
#' @examples
#' # path to example config file
#' config <- system.file("extdata", "config.csv", package = "amplican")
#' # path to example fastq files
#' fastq_folder <- system.file("extdata", package = "amplican")
#' # output folder
#' results_folder <- tempdir()
#' amplicanAlign(config, fastq_folder, results_folder)
#'
amplicanAlign <- function(config,
                          fastq_folder,
                          results_folder,
                          total_processors = 1,
                          skip_bad_nucleotides = TRUE,
                          average_quality = 30,
                          min_quality = 20,
                          write_alignments = TRUE,
                          scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
                            match = 5, mismatch = -4,
                            baseOnly = TRUE, type = "DNA"),
                          gap_opening = 50,
                          gap_extension = 0,
                          fastqfiles = 0,
                          PRIMER_DIMER = 30,
                          cut_buffer = 5) {

  message("Checking write access...")
  checkFileWriteAccess(results_folder)

  message("Checking configuration file...")
  configTable <- utils::read.csv(config, strip.white = TRUE, stringsAsFactors = FALSE)
  configTable[is.na(configTable)] <- ""
  colnames(configTable) <- c("ID",
                             "Barcode",
                             "Forward_Reads_File",
                             "Reverse_Reads_File",
                             "Group",
                             "guideRNA",
                             "Forward_Primer",
                             "Reverse_Primer",
                             "Direction",
                             "Amplicon")
  configTable$Barcode <- as.character(configTable$Barcode)
  configTable$Forward_Reads_File <- ifelse(configTable$Forward_Reads_File == "",
                                           "",
                                           file.path(fastq_folder, configTable$Forward_Reads_File))
  configTable$Reverse_Reads_File <- ifelse(configTable$Reverse_Reads_File == "",
                                           "",
                                           file.path(fastq_folder, configTable$Reverse_Reads_File))

  if (sum(configTable$Reverse_Reads_File == "") > 0) {
    message("Reverse_Reads_File has empty rows. Changing fastqfiles parameter to 1, operating only on forward reads.")
    fastqfiles <- 1
  }
  if (sum(configTable$Forward_Reads_File == "") > 0) {
    message("Forward_Reads_File has empty rows. Changing fastqfiles parameter to 2, operating only on reverse reads.")
    fastqfiles <- 2
  }
  checkConfigFile(configTable, fastq_folder)
  configTable$Reverse_PrimerRC <- revComp(configTable$Reverse_Primer)
  configTable <- checkPrimers(configTable, fastqfiles)

  revDir <- configTable[, "Direction"] == 1
  configTable[revDir, "guideRNA"] <- revComp(configTable[revDir, "guideRNA"])
  configTable$RguideRNA <- revComp(configTable$guideRNA)
  configTable$Found_Guide <- checkTarget(configTable)
  configTable$cutSites <- lapply(
    configTable$Amplicon, function(x) upperGroups(x) + cut_buffer)
  configTable$ampl_len <- nchar(configTable$Amplicon)

  cutSitesCheck <- sapply(configTable$cutSites, length) == 0
  if (any(cutSitesCheck)) {
    message("Warning: Config file row without upper case groups (PAM): ",
            toString(which(cutSitesCheck)))
    configTable$cutSites[cutSitesCheck] <- as.list(
      IRanges::tile(
        IRanges::IRanges(start = 1,
                         width = configTable$ampl_len[cutSitesCheck]),
        1))
  }

  resultsFolder <- file.path(results_folder, "alignments")
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  unassignedFolder <- file.path(resultsFolder, "unassigned_sequences")
  if (!dir.exists(unassignedFolder)) {
    dir.create(file.path(unassignedFolder))
  }

  # Parameters
  logFileName <- file.path(results_folder, "RunParameters.txt")
  if (file.exists(logFileName)) {
    file.remove(logFileName)
  }
  logFileConn <- file(logFileName, open = "at")
  writeLines(c(paste("Config file:           ", config),
               paste("Processors used:       ", total_processors),
               paste("Skip Bad Nucleotides:  ", skip_bad_nucleotides),
               paste("Average Quality:       ", average_quality),
               paste("Minimum Quality:       ", min_quality),
               paste("Write Alignments:      ", write_alignments),
               paste("Fastq files Mode:      ", fastqfiles),
               paste("Gap Opening:           ", gap_opening),
               paste("Gap Extension:         ", gap_extension),
               paste("PRIMER DIMER buffer:   ", PRIMER_DIMER),
               paste("Cut buffer:            ", cut_buffer),
               "Scoring Matrix:"), logFileConn)
  write.csv(scoring_matrix, logFileConn, quote = FALSE, row.names = TRUE)
  close(logFileConn)


  uBarcode <- unique(configTable$Barcode)
  configTable$ExperimentsCount <- table(configTable$Barcode)[configTable$Barcode]

  # Several statistics about deletions, cuts and reads
  configTable$Cut <- 0
  configTable$Frameshift <- 0
  configTable$PRIMER_DIMER <- 0
  configTable$Reads <- 0

  # Warnings
  configTable$Found_Guide <- 0
  configTable$Found_PAM <- !cutSitesCheck

  if (total_processors > 1) {
    cl <- parallel::makeCluster(total_processors, outfile = "")
    doParallel::registerDoParallel(cl)
    p <- DoparParam()

    configSplit <- split(configTable, f = configTable$Barcode)
    bplapply(configSplit, FUN = makeAlignment,
             resultsFolder,
             skip_bad_nucleotides,
             average_quality,
             min_quality,
             write_alignments,
             scoring_matrix,
             gap_opening,
             gap_extension,
             fastqfiles,
             PRIMER_DIMER,
             cut_buffer, BPPARAM=p)

    parallel::stopCluster(cl)

  } else {
    for (j in seq_along(uBarcode)) {
      makeAlignment(configTable[configTable$Barcode == uBarcode[j], ],
                    resultsFolder,
                    skip_bad_nucleotides,
                    average_quality,
                    min_quality,
                    write_alignments,
                    scoring_matrix,
                    gap_opening,
                    gap_extension,
                    fastqfiles,
                    PRIMER_DIMER,
                    cut_buffer)
    }
  }

  message("Alignments done. Creating results files...")

  # Put all the logs and all the configs together
  unifyFiles(resultsFolder, "configFile_results",
             file.path(results_folder, "config_summary.csv"))
  unifyFiles(resultsFolder, "alignment_ranges",
             file.path(results_folder, "alignments_events.csv"))
  unifyFiles(resultsFolder, "reads_filters.csv",
             file.path(results_folder, "barcode_reads_filters.csv"))
  message("Finished.")

  invisible(results_folder)
}
