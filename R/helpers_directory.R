#' Checks if the given directory exist and can be written to.
#'
#' @param filePath (string) A string the path to the file.
#' @return (invisible) TRUE, Stop if no access.
#'
checkFileWriteAccess <- function(filePath) {
  if (file.access(filePath[1], mode = 2) != 0) {
    stop("No write access to the path or it does not exists: ",
         filePath)
  }
  invisible(TRUE)
}
