#' This function checks if the target RNA is in the amplicon.
#'
#' @param targetPrimer (string) A sequence of nucleotides in a string format
#' representing the target
#' @param amplicon (string) A sequence of nucleotides in a string format
#' representing the amplicon
#' @param ID (string) The ID from where this target and amplicon came.
#' @param barcode (string) The barcode from where this target and amplicon came.
#' @param logFileConn (connection) A file descriptor in R which is pointing to
#' the log file.
#' @return (boolean) TRUE when target is in amplicon
#'
checkTarget <- function(targetPrimer, amplicon, ID, barcode, logFileConn) {
    targetPositions <- grepl(targetPrimer, amplicon, ignore.case = TRUE)
    if (!targetPositions) {
        message("Warning: guideRNA has not been found in the amplicon.
                Check the log file for more information.")
        writeLines(paste0("Couldn't find the guideRNA in amplicon: ",
                          targetPrimer,
                          "\nFor ID: ",
                          ID,
                          " and barcode: ",
                          barcode,
                          "\n"),
                   logFileConn)
    }
    return(targetPositions)
}

#' This function checks if the forward and reverse primer are in the amplicon.
#'
#' The forward primer should be the one coming from the configuration file,
#' while the reverse primer should be reverse and complementary.
#' This function DOES NOT reverse complements the input.
#'
#' @param forwardPrimer (string) A sequence of nucleotides in a string format
#' representing the forward primer
#' @param reversePrimerRC (string) A sequence of nucleotides in a string format
#' representing the reverse primer. Usually you want to reverse complement this
#' BEFORE giving it to the function. The function WILL NOT reverse complement
#' it for you.
#' @param amplicon (string) A sequence of nucleotides in a string format
#' representing the amplicon
#' @param ID (string) The ID from where this target and amplicon came.
#' @param barcode (string) The barcode from where this target and amplicon came.
#' @param logFileConn (connection) The location on disk where you can find the
#' config file with this ID and barcode.
#' @return (boolean) TRUE when both primers are found in the amplicon
#'
checkPrimers <- function(forwardPrimer, reversePrimerRC, amplicon, ID,
                         barcode, logFileConn) {
    forwardPrimerPosition <- grepl(forwardPrimer, amplicon, ignore.case = TRUE)
    reversePrimerPosition <- grepl(reversePrimerRC, amplicon, ignore.case = TRUE)
    if (!(forwardPrimerPosition | reversePrimerPosition)) {
        message("Warning: One of primer was not found in the amplicon.
                Check the log file for more information.")
        writeLines(paste0("Couldn't find the forward primer: ",
                          toString(forwardPrimer),
                          "\nor reverse primer: ",
                          toString(reversePrimerRC),
                          "\nFor ID: ",
                          ID,
                          " and barcode: ",
                          barcode,
                          "\n"),
                   logFileConn)
    }
    return(forwardPrimerPosition | reversePrimerPosition)
}


#' This function pre-process a config file and checks that everything is in
#' order.
#'
#' Its takes care of the following:
#'   No IDs are duplicated.
#'   Every combination of barcode, forward primer and reverse primer is unique.
#'   Each barcode has unique forward reads file and reverse read files.
#'   Checks that the read files exist with read access.
#' @param configTable (data.frame) Config file.
#' @param fastq_folder (string) Path to fastq folder.
#' @return (void) If anything goes wrong stops and prints error.
#' @importFrom stats complete.cases
#'
checkConfigFile <- function(configTable, fastq_folder) {

    totalRows <- dim(configTable)[1]
    totalCols <- dim(configTable)[2]

    goodRows <- stats::complete.cases(configTable)
    if (sum(goodRows) != totalRows) {
        stop(paste0("Config file has bad rows: ",
                    paste(which(goodRows == FALSE) + 1),
                    " due to NA/NULL values"))
    }

    if (length(unique(configTable[, "ID"])) != totalRows) {
        stop(paste0("Config file has duplicates IDs in rows:",
                    paste(which(duplicated(configTable[, "ID"])) + 1)))
    }

    barcode_primers_duple <- duplicated(configTable[, c("Barcode",
                                                        "Forward_Primer",
                                                        "Reverse_Primer")])
    if (sum(barcode_primers_duple) != 0) {
        stop(paste0("Config file has non unique combinations of barcode,
                    forward primer and reverse primer. Duplicated rows: ",
                    toString(which(barcode_primers_duple) + 1)))
    }

    barcode_files_duple <- duplicated(configTable[c("Barcode")])
    forward_reverse_files_duple <- duplicated(configTable[c("Forward_Reads_File", "Reverse_Reads_File")])
    fail_barcodes <- which(barcode_files_duple != forward_reverse_files_duple) + 1
    if (length(fail_barcodes) > 0) {
        stop(paste0("Each of these rows are malfunctioned in the config file: ",
                    toString(fail_barcodes),
                    " For each barcode there can be only one set of paths for forward and reverse files."))
    }

    uniqueFilePaths <- unique(c(as.character(configTable$Forward_Reads_File),
                                as.character(configTable$Reverse_Reads_File)))
    access <- file.access(uniqueFilePaths, mode = 4) == -1
    if (sum(access) > 0) {
        stop(paste0("We either don't have read access or paths are incorrect. ",
                    "Check specified paths in config for these files:\n",
                    paste(uniqueFilePaths[access],
                          sep = "\n")))
    }
}
