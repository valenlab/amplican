#' Filter Events Overlapping Primers and PRIMER DIMERS.
#'
#' Very often alignments return deletions that are not real deletions, but
#' rather artifact of incomplete reads eg.: \cr
#' \preformatted{
#' ACTGAAAAA------- <- this "deletion" should be filtered
#' ACTG----ACTGACTG
#' }
#' We call them Events Overlapping Primers and filter them together
#' with reads that are potentially PRIMER DIMERS.
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end.
#' @param cfgT (data.frame) Needs columns Forward_Primer, ReversePrimer and
#' Amplicon.
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @return (aln) Reduced by events classified as PRIMER DIMER or overlapping
#' primers.
#' @export
#' @family analysis steps
#' @seealso \code{\link{findPD}} and \code{\link{findEOP}}
#' @include helpers_filters.R
#'
amplicanFilter <- function(aln, cfgT, PRIMER_DIMER) {
  PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

  # PRIMER DIMER reads with unique ID and read_id
  onlyPD <- aln[PD, c("seqnames", "read_id")]
  onlyPD <- onlyPD[!duplicated(onlyPD), ]

  # apply filter
  PD <- which(aln$seqnames %in% onlyPD$seqnames)
  read_ids <- onlyPD$read_id[match(aln$seqnames[PD], onlyPD$seqnames)]
  PD <- PD[aln$read_id[PD] == read_ids]
  aln <- aln[-PD,]

  # filter events overlapping primers
  eOP <- findEOP(aln, cfgT)
  aln[!eOP, ]
}
