#' This function check if the target RNA is in the amplicon.
#'
#'@param targetPrimer (string) A sequence of nucleotides in a string format representing the target
#'@param amplicon (string) A sequence of nucleotides in a string format representing the amplicon
#'@param ID (string) The ID from where this target and amplicon came.
#'@param barcode (String) The barcode from where this target and amplicon came.
#'@param logFileConn (connection) A file descriptor in R which is pointing to the log file.
#'@return (bool) , TRUE:  Everything went good.
#'
checkTarget <- function(targetPrimer, amplicon, ID, barcode, logFileConn){
  targetPositions <- grepl(targetPrimer, amplicon, ignore.case=TRUE)
  if (targetPositions) {
    warning("Target has not been found in the amplicon. Check the log file for more information.")
    writeLines("Couldn't find the target primer:", logFileConn)
    writeLines(targetPrimer, logFileConn)
    writeLines("In amplicon:", logFileConn)
    writeLines(amplicon, logFileConn)
    writeLines(paste0("For ID: ", ID, " and barcode: ", barcode), logFileConn)
    writeLines("\n", logFileConn)
  }
  return(targetPositions)
}

#' This function check if the forward anr reverse primer is in the amplicon.
#'
#' The forward primer should be the one coming from the configuration file, while the reverse primer should be
#' reverse and complement. This function DOES NOT reverse complement the input. So if you look up for something
#' that is suppose to containt a reverse complement it won't find it. However you can use it to find any pair of
#' primers regarless if they are reverse complemented or not.
#'
#' @param forwardPrimer (string) A sequence of nucleotides in a string format representing the forward primer
#' @param reversePrimerRC (string) A sequence of nucleotides in a string format representing the reverse primer. Usually
#'                         you want to reverse complement this BEFORE giving it to the function. The function WILL NOT
#'                         reverse complement it for you. (see: /libraries/tools.R , you have there the reverse
#'                         complement function)
#' @param amplicon (string) A sequence of nucleotides in a string format representing the amplicon
#' @param ID (string) The ID from where this target and amplicon came.
#' @param barcode (string) The barcode from where this target and amplicon came.
#' @param logFileConn (connection) The location on disk where you can find the config file with this ID and barcode.
#' @return (bool) TRUE when everything went good.
#'
checkPrimers <- function(forwardPrimer, reversePrimerRC, amplicon, ID, barcode, logFileConn){
  forwardPrimerPosition <- grepl(forwardPrimer, amplicon, ignore.case=TRUE)
  reversePrimerPosition <- grepl(reversePrimerRC, amplicon, ignore.case=TRUE)
  if (!(forwardPrimerPosition | reversePrimerPosition)) {
    warning("One of primer was not found in the amplicon. Check the log file for more information.")
    writeLines("Couldn't find the forward primer or reverse primer (reversed and complemented):", logFileConn)
    writeLines(toString(forwardPrimer), logFileConn)
    writeLines(toString(reversePrimerRC), logFileConn)
    writeLines("In amplicon:", logFileConn)
    writeLines(toString(amplicon), logFileConn)
    writeLines(paste0("For ID: ", ID, " and barcode: ", barcode), logFileConn)
    writeLines("\n", logFileConn)
  }
  return(forwardPrimerPosition | reversePrimerPosition)
}

#' This function checks if the given alignment positions are valid
#'
#' @param alignmentPositions (vector) of integers representing the aligning positions inside the amplicon
#' @param amplicon (string) A sequence of nucleotides in a string format representing the amplicon
#' @param ID (string) The ID for target and amplicon.
#' @param barcode (string) barcode belonging to target and amplicon.
#' @param logFileConn (connection) Path to config file.
#' @return (bool) TRUE when positions are greater than 0.
#'
checkPositions <- function(alignmentPositions, amplicon, ID, barcode, logFileConn){
  if (alignmentPositions[1] <= 0) {
    warning("Aligment position was not found in the amplicon. Find more information in the log file.")
    writeLines("Couldn't find alignment position for amplicon:", logFileConn)
    writeLines(toString(amplicon), logFileConn)
    writeLines("For ID and barcode:", logFileConn)
    writeLines(toString(paste(ID, barcode, sep=" ")), logFileConn)
    writeLines("\n", logFileConn)
  }
  return(alignmentPositions[1] > 0)
}

#' This function pre-process a config file and check that everything is in order
#'
#' Its takes care of the following:
#'   No IDs are duplicated
#'   Every combination of barcode, forward primer and reverse primer is unique.
#'   Each barcode has unique forward reads file and reverse read files.
#'   Check that the read files exist with read access.
#' @param configTable (data.frame) Config file.
#' @param fastq_folder (string) Path to fastq folder.
#' @return (Void) If anything goes wrong stops and prints error.
checkConfigFile <- function(configTable, fastq_folder){

  totalRows <- dim(configTable)[1]
  totalCols <- dim(configTable)[2]

  goodRows <- complete.cases(configTable)
  if (sum(goodRows) != totalRows) {
    stop(paste0("Config file has bad rows: ",
                paste(which(goodRows == F)),
                " due to NA/NULL values"))
  }

  if(length(unique(configTable[, "ID"])) != totalRows){
    stop(paste0("Config file has duplicates IDs in rows:",
                paste(which(duplicated(configTable[, "ID"])))))
  }

  barcode_primers_duple <- duplicated(configTable[c("Barcode", "Forward_Primer", "Reverse_Primer")])
  if(sum(barcode_primers_duple) != 0){
    stop(paste0("Config file has non unique combinations of barcode, forward primer and reverse primer.
                Duplicated rows: ", paste(which(barcode_primers_duple))))
  }

  barcode_files_duple <- duplicated(configTable[c("Barcode")])
  forward_reverse_files_duple <- duplicated(configTable[c("Forward_Reads_File", "Reverse_Reads_File")])
  fail_barcodes <- which(barcode_files_duple != forward_reverse_files_duple)
  if (length(fail_barcodes) > 0) {
    stop(paste0("Each of these rows are malfunctioned in the config file: ", paste(fail_barcodes),
                " For each barcode there can be only one set of paths for forward and reverse files."))
  }

  uniqueFilePaths <- unique(c(as.character(configTable$Forward_Reads_File),
                              as.character(configTable$Reverse_Reads_File)))
  access <- file.access(uniqueFilePaths, mode = 4) == -1
  if (sum(access) > 0) {
    stop(paste0("We either dont have read access or paths are incorrect. ",
                "Check specified paths in config for these files:\n",
                paste(uniqueFilePaths[access], sep = "\n")))
  }
}
