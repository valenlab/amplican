#' Merge files of given type.
#'
#' This function takes all the files in a given folder that match a regex
#' pattern and merges them together.
#'
#' @param targetFolder (string) The folder where the files are suppose to be.
#' @param regex (string) The files to be combined, must comply with this regex.
#' @param finalFileName (String) The name of the final file.
#' @param header (boolean) If header is true, it will take the header of the
#'                     first file only and ignore the rest. Otherwise it
#'                     will repeat the header many times.
#' @param delete (boolean) If true, it will delete the original files
#' @param isrecursive (boolean) If TRUE will search recursively through regex
#' matching folders.
#' @return (void)
#'
unifyFiles <- function(targetFolder, regex, finalFileName, header = TRUE,
                       delete = TRUE, isrecursive = FALSE) {

  allFiles <- list.files(targetFolder, recursive = isrecursive)
  candidates <- grep(regex, allFiles, fixed = FALSE, ignore.case = FALSE)

  if (length(candidates) > 0) {
    if (file.exists(finalFileName)) {
      file.remove(finalFileName)
    }
    finalFileFD <- file(finalFileName, open = "w")

    for (i in seq_along(candidates)) {
      candidateName <- file.path(targetFolder, allFiles[candidates[i]])
      # Read the file into a string variable
      textFile <- readLines(candidateName, encoding = "UTF-8")
      if (header == TRUE && i != 1) {
        textFile <- textFile[2:length(textFile)]
      }
      writeLines(textFile, finalFileFD)

      if (delete == TRUE) {
        file.remove(candidateName)
      }
    }
    # Close the file descriptor for the final file
    close(finalFileFD)
  }
  return()
}

#' Remove forward and reverse fastq files.
#'
#' Delete files from configTables Forward_Reads_File and
#' Reverse_Reads_File (fastq files).
#'
#' @param configTable (data.frame) Contains Forward_Reads_File and
#' Reverse_Reads_File to be removed
#' @return (void) In case of fail, prints err.
#'
deleteFiles <- function(configTable) {
  for (i in seq_len(dim(configTable)[1])) {
    file.remove(configTable$Forward_Reads_File[i])
    file.remove(configTable$Reverse_Reads_File[i])
  }
}

#' This function checks if the given directory exist and can be written to.
#'
#' @param filePath (string) A string the path to the file.
#' @return (void) Stop if no access.
#'
checkFileWriteAccess <- function(filePath) {
  if (file.access(filePath[1], mode = 2) != 0) {
    stop(paste0("No write access to the path or it does not exists: ",
                filePath))
  }
}
