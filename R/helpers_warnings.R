#' This function checks if the guideRNA is in the amplicon.
#'
#' @param configTable (data.frame) data frame of config file
#' @return (boolean vector) Prints warning when some guides can't be found.
#' @importFrom stringr str_detect
#'
checkTarget <- function(configTable) {
  targetPositions <- stringr::str_detect(tolower(configTable$Amplicon),
                                         tolower(configTable$guideRNA))
  if (any(!targetPositions)) {
    message(paste0(
      "Warning: guideRNA has not been found in the amplicon for line: ",
      toString(which(!targetPositions))
    ))
  }

  return(targetPositions)
}

#' This function checks if the forward and reverse primer are in the amplicon
#' and where they are located.
#'
#' @param configTable (data.frame) A data frame of config file.
#' @param fastqfiles (numeric) Which primers are important.
#' @return configTable (data.frame) A data frame of config file with additional
#' fields for start locations of the primers
#' @importFrom stringr str_locate
#'
checkPrimers <- function(configTable, fastqfiles) {

  configTable$forwardPrimerPosition <-
    stringr::str_locate(tolower(configTable$Amplicon),
                        tolower(configTable$Forward_Primer))[, 1]
  configTable[, c("reversePrimerPosition", "reversePrimerPosEnd")] <-
    stringr::str_locate(tolower(configTable$Amplicon),
                        tolower(configTable$Reverse_PrimerRC))

  fP <- if (fastqfiles != 2) which(is.na(configTable$forwardPrimerPosition)) else NULL
  rP <- if (fastqfiles != 1) which(is.na(configTable$reversePrimerPosition)) else NULL

  if (length(fP) > 0 | length(rP) > 0) {
    stop(paste0(
        "Error: One of primers was not found in the amplicon.",
        if (length(fP) > 0)
          paste0(" Could't locate forward primer in amplicon for row: ",
                 toString(fP))
        else "",
        if (length(rP) > 0)
          paste0(" Could't locate reverse primer in amplicon for row: ",
                 toString(rP))
        else ""
      )
    )
  }

  return(configTable)
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
#' @return (boolean) TRUE, If anything goes wrong stops and prints error.
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
    stop(paste0("Config file has non unique combinations of barcode, forward primer and reverse primer. Duplicated rows: ",
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
  invisible(TRUE)
}
