#' This filters out sequences which have bad base quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param min (Int) This is the minimum quality that we accept for
#' every nucleotide. For example, if we have a sequence with nucleotides which have quality
#' 50-50-50-50-10, and we set the minimum to 30, the whole sequence will be a bad sequence.
#' The minimum is set to 0 by default.
#' @return (Logical) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead ShortReadQ
#'
goodBaseQuality <- function(reads, min = 0){
  apply(as(slot(reads, "quality"), "matrix"), 1, min, na.rm = T) >= min
}


#' This filters out sequences which have bad average quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param avg (Int) This is what the average score of the quality of
#' sequence should be (or greater). For example, if we have a sequence with nucleotides which have
#' quality 70-70-70, the average would be 70. If set the average to 70 or less the sequence will
#' pass. If we set the average to 71 the sequence will not pass. The average is set to 0 by default.
#' @return (Logical) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead ShortReadQ
#'
goodAvgQuality <- function(reads, avg = 0){
  apply(as(slot(reads, "quality"), "matrix"), 1, mean, na.rm = T) >= avg
}


#' This filters out sequences which have bad average quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @return (Logical) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead srFilter
#'
alphabetQuality <- function(reads){
  nucq <- polynFilter(nuc = c("A", "C", "T", "G"))
  !as.logical(nucq(reads))
}
