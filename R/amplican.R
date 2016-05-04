#' amplican: A package for something.
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#'
#' @section amplican functions:
#' The amplican functions ...
#'
#' @docType package
#' @name amplican
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("amplican", libpath)
}
