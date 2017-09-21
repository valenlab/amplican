#' Pretty print forward and reverse reads aligned to each other.
#'
#' Usefull and needed for barcode reports.
#'
#' @param forward (character or vector of characters) Forward reads.
#' @param reverse (character or vector of characters) Will be reverse
#' complemented before alignment.
#' @return Vector with alignments ready to be printed.
#' @include helpers_general.R
#' @export
#' @examples
#' # load example data
#' unassigned_file <- system.file('extdata', 'results',  'alignments',
#'                                'unassigned_reads.csv', package = 'amplican')
#' unassigned <- data.table::setDF(data.table::fread(unassigned_file))
#' # sort by frequency
#' unassigned <- unassigned[order(unassigned$BarcodeFrequency,
#'                                decreasing = TRUE), ]
#' # print alignment of most frequent unassigned reads
#' cat(amplican_print_reads(unassigned[1, 'Forward'],
#'                          unassigned[1, 'Reverse']),
#'           sep = "\n")
#'
amplican_print_reads <- function(forward, reverse) {

  alignForwardReverse <- Biostrings::pairwiseAlignment(forward,
                                                       revComp(reverse))
  wPA <- utils::capture.output(
    Biostrings::writePairwiseAlignments(alignForwardReverse))
  wPA <- wPA[!grepl("#", wPA)]

  wPA
}


#' Helper function to calculate figure height
#' based on number of elements to plot.
#'
#' @param x (numeric) number of elements to fit onto height axis
#' @return (numeric) In inches
#' @export
#' @examples
#' plot_height(20)
#'
plot_height <- function(x) {
    x <- grid::unit(x, "char") + grid::unit(x * 0.2, "mm") + grid::unit(2, "cm")
    ceiling(as.numeric(grid::convertUnit(x, "inches")))
}
