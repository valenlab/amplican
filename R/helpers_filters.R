#' Find Events Overlapping Primers.
#'
#' Very often alignments return deletions that are not real deletions, but
#' rather artifact of incomplete reads eg.: \cr
#' \preformatted{
#' ACTGAAAAA------- <- this "deletion" should be filtered
#' ACTG----ACTGACTG
#' }
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end.
#' @param cfgT (data.frame) Needs columns Forward_Primer, ReversePrimer and
#' Amplicon.
#' @return (logical vector) where TRUE indicates events that are overlapping
#' primers
#' @export
#' @family filters
#' @seealso \code{\link{findPD}}
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' findEOP(aln, cfgT)
#'
findEOP <- function(aln, cfgT) {
  mapID <- match(aln$seqnames, cfgT$ID)
  (aln$start < cfgT$fwdPrPosEnd[mapID]) | (aln$end > cfgT$rvePrPos[mapID])
}


#' Find PRIMER DIMER reads.
#'
#' Use to filter reads that are most likely PRIMER DIMERS.
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end.
#' @param cfgT (data.frame) Needs columns Forward_Primer, ReversePrimer and
#' Amplicon.
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @return (logical) Where TRUE indicates event classified as PRIMER DIMER
#' @export
#' @family filters
#' @seealso \code{\link{findEOP}}
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' findPD(aln, cfgT)
#'
findPD <- function(aln, cfgT, PRIMER_DIMER = 30) {

  PD_cutoff <- nchar(cfgT$Amplicon) -
    (nchar(cfgT$Forward_Primer) + nchar(cfgT$Reverse_Primer) + PRIMER_DIMER)
  PD_cutoff <- PD_cutoff[match(aln$seqnames, cfgT$ID)]

  aln$width > PD_cutoff
}


#' Filters out sequences which have bad base quality readings.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param min (numeric) This is the minimum quality that we accept for
#' every nucleotide. For example, if we have a sequence with nucleotides which
#' have quality 50-50-50-50-10, and we set the minimum to 30, the whole sequence
#' will be a bad sequence. The minimum is set to 0 by default.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#'
goodBaseQuality <- function(reads, min = 20) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(matrixStats::rowMins(methods::as(quality(reads), "matrix"),
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
#'
goodAvgQuality <- function(reads, avg = 30) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(Matrix::rowMeans(methods::as(quality(reads), "matrix"),
                          na.rm = TRUE) >= avg)
}


#' This filters out sequences which have nonstandard nucleotides.
#'
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#'
alphabetQuality <- function(reads) {
  if (is.logical(reads)) {
    return(reads)
  }

  nucq <- ShortRead::nFilter()
  return(as.logical(nucq(reads)))
}
