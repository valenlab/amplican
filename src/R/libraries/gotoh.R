library(Rcpp) # Alignments
sourceCpp("gotoh3.cpp") # Alignments.

#' A Gotoh Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' gotoh_function()

gotoh <- function(pattern, subject, score_matrix, gap_opening, gap_extension, gap_ending, far_indels){

  return (nwRCPP(pattern, subject,score_matrix, gap_opening, gap_extension, gap_ending, far_indels))


}
  
  