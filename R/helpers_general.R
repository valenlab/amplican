#' Reverse and complement given string or list of strings
#'
#' @param x (string or vector of strings)
#' @return (string or vector of strings) reverse complemented input
#' @importFrom Biostrings DNAStringSet reverseComplement
#'
revComp <- function(x) {
  return(
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  )
}
