#' Pretty print forward and reverse reads aligned to each other.
#'
#' Usefull and needed for barcode reports.
#'
#' @param forward (character or vector of characters) Forward reads.
#' @param reverse (character or vector of characters) Will be reverse
#' complemented before alignment.
#' @return Vector with alignments ready to be printed.
#' @import Biostrings
#' @importFrom utils capture.output
#' @include helpers_general.R
#' @export
#'
amplican_print_reads <- function(forward, reverse) {

  alignForwardReverse <- Biostrings::pairwiseAlignment(forward,
                                                       revComp(reverse))
  wPA <- capture.output(
    Biostrings::writePairwiseAlignments(alignForwardReverse))
  wPA <- wPA[!grepl("#", wPA)]

  wPA
}


#' Helper function to calculate figure height
#' based on number of elements to plot.
#'
#' @param x (numeric) number of elements to fit onto height axis
#' @return (numeric)
#' @export
#'
plot_height <- function(x) {
  if (x < 10) {
    x + 1
  } else {
    ceiling(9 * log(x))
  }
}
