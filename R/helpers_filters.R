#' This filters out sequences which have bad base quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param min (numeric) This is the minimum quality that we accept for
#' every nucleotide. For example, if we have a sequence with nucleotides which
#' have quality 50-50-50-50-10, and we set the minimum to 30, the whole sequence
#' will be a bad sequence. The minimum is set to 0 by default.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead ShortReadQ
#' @importFrom matrixStats rowMins
#' @importFrom methods as slot
#'
goodBaseQuality <- function(reads, min = 30) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(matrixStats::rowMins(as(slot(reads, "quality"), "matrix"),
                              na.rm = TRUE) >= min)
}


#' This filters out sequences which have bad average quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param avg (numeric) This is what the average score of the quality of
#' sequence should be. For example, if we have a sequence with nucleotides which
#' have quality 70-70-70, the average would be 70. If set the average to 70 or
#' less the sequence will pass. If we set the average to 71 the sequence will
#' not pass. The average is set to 0 by default.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead ShortReadQ
#' @importFrom Matrix rowMeans
#' @importFrom methods as slot
#'
goodAvgQuality <- function(reads, avg = 30) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(Matrix::rowMeans(as(slot(reads, "quality"), "matrix"),
                          na.rm = TRUE) >= avg)
}


#' This filters out sequences which have nonstandard nucleotides.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#' @importFrom ShortRead ShortReadQ sread
#'
alphabetQuality <- function(reads) {
  if (is.logical(reads)) {
    return(reads)
  }

  nucq <- ShortRead::nFilter()
  return(as.logical(nucq(reads)))
}
