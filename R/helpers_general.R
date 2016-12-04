#' Reverse complement string or vector of strings.
#'
#' @param toRevComp (vector of characters) Can be singular character.
#' @return (vector of characters)
#' @importFrom Biostrings DNAString DNAStringSet reverseComplement
#'
revC <- function(toRevComp) {

  as.character(Biostrings::reverseComplement(
    if (length(toRevComp) == 1) {
      Biostrings::DNAString(toRevComp)
    } else {
      Biostrings::DNAStringSet(toRevComp)
    }))
}
