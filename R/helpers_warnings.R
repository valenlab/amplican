#' Checks if the guideRNA is in the amplicon.
#'
#' @keywords internal
#' @param configTable (data.frame) data frame of config file
#' @return (boolean vector) Prints warning when some guides can't be found.
#'
checkTarget <- function(configTable) {
  targetPositions <- stringr::str_detect(tolower(configTable$Amplicon),
                                         tolower(configTable$guideRNA))
  if (any(!targetPositions)) {
    message(
      "Warning: guideRNA has not been found in the amplicon for line: ",
      toString(which(!targetPositions) + 1)
    )
  }

  return(targetPositions)
}

#' Checks if the forward and reverse primer are in the amplicon
#' and where they are located.
#'
#' @keywords internal
#' @param configTable (data.frame) A data frame of config file.
#' @param fastqfiles (numeric) Which primers are important.
#' @return configTable (data.frame) A data frame of config file with additional
#' fields for start locations of the primers
#'
checkPrimers <- function(configTable, fastqfiles) {

  # str_locate is 1 based (min is 1)
  configTable[, c("fwdPrPos", "fwdPrPosEnd")] <-
    stringr::str_locate(tolower(configTable$Amplicon),
                        tolower(configTable$Forward_Primer))
  configTable[, c("rvePrPos", "rvePrPosEnd")] <-
    stringr::str_locate(tolower(configTable$Amplicon),
                        tolower(configTable$Reverse_PrimerRC))

  fP <- if (fastqfiles != 2) which(is.na(configTable$fwdPrPos)) else NULL
  rP <- if (fastqfiles != 1) which(is.na(configTable$rvePrPos)) else NULL

  if (length(fP) > 0 | length(rP) > 0) {
    stop("One of primers was not found in the amplicon.",
         if (length(fP) > 0)
           c(" Couldn't locate forward primer in amplicon for row: ",
             toString(fP + 1))
         else "",
         if (length(rP) > 0)
           c(" Couldn't locate reverse primer in amplicon for row: ",
             toString(rP + 1))
         else ""
    )
  }

  return(configTable)
}


#' Pre-process a config file and checks that everything is in
#' order.
#'
#' Its takes care of the following:
#'   No IDs are duplicated.
#'   Every combination of barcode, forward primer and reverse primer is unique.
#'   Each barcode has unique forward reads file and reverse read files.
#'   Checks that the read files exist with read access.
#' @keywords internal
#' @param configTable (data.frame) Config file.
#' @param fastq_folder (string) Path to fastq folder.
#' @return (boolean) TRUE, If anything goes wrong stops and prints error.
#'
checkConfigFile <- function(configTable, fastq_folder) {

  totalRows <- dim(configTable)[1]
  totalCols <- dim(configTable)[2]

  rp_num <- grepl("\\d", configTable$Forward_Primer)
  if (any(rp_num)) {
    stop("Config file has bad rows: ",
         toString(which(rp_num) + 1),
         " due to reverse primers containing numeric values.")
  }

  fp_num <- grepl("\\d", configTable$Reverse_Primer)
  if (any(fp_num)) {
    stop("Config file has bad rows: ",
         toString(which(fp_num) + 1),
         " due to forward primers containing numeric values.")
  }

  goodRows <- stats::complete.cases(configTable)
  if (sum(goodRows) != totalRows) {
    stop("Config file has bad rows: ",
         toString(which(goodRows == FALSE) + 1),
         " due to NA/NULL values")
  }

  if (length(unique(configTable$ID)) != totalRows) {
    stop("Config file has duplicates IDs in rows: ",
         toString(which(duplicated(configTable[, "ID"])) + 1))
  }

  barcode_primers_duple <- duplicated(configTable[, c("Barcode",
                                                      "Forward_Primer",
                                                      "Reverse_Primer")])
  if (sum(barcode_primers_duple) != 0) {
    stop("Config file has non unique combinations of barcode, forward ",
         "primer and reverse primer. Duplicated rows: ",
         toString(which(barcode_primers_duple) + 1))
  }

  barcode_files_duple <- duplicated(configTable$Barcode)
  forward_reverse_files_duple <- duplicated(
    configTable[c("Forward_Reads_File", "Reverse_Reads_File")])
  fail_barcodes <- which(barcode_files_duple != forward_reverse_files_duple) + 1
  if (length(fail_barcodes) > 0) {
    stop("Each of these rows are malfunctioned in the config file: ",
          toString(fail_barcodes),
          " For each barcode there can be only one set of paths ",
          "for forward and reverse files.")
  }

  uniqueFilePaths <- unique(
    c(configTable$Forward_Reads_File[configTable$Forward_Reads_File != ""],
      configTable$Reverse_Reads_File[configTable$Reverse_Reads_File != ""]))
  access <- file.access(uniqueFilePaths, mode = 4) == -1
  if (sum(access) > 0) {
    stop("We either don't have read access or paths are incorrect. ",
         "Check specified paths in config for these files:\n",
         toString(uniqueFilePaths[access]))
  }
  invisible(TRUE)
}
