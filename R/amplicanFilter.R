#' Filter Events Overlapping Primers, PRIMER DIMERS and Low Alignment Score
#' Events.
#'
#' Very often alignments return deletions that are not real deletions, but
#' rather artifact of incomplete reads eg.: \cr
#' \preformatted{
#' ACTGAAAAA------- <- this "deletion" should be filtered
#' ACTG----ACTGACTG
#' }
#' We call them Events Overlapping Primers and filter them together
#' with reads that are potentially PRIMER DIMERS. This filter will also remove
#' all events coming from reads with low alignment score - potential
#' Off-targets.
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
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' amplicanFilter(aln, cfgT, 30)
#'
amplicanFilter <- function(aln, cfgT, PRIMER_DIMER) {
  seqnames <- NULL

  eOP <- findEOP(aln, cfgT)
  aln <- aln[!eOP, ]

  PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

  # PRIMER DIMER reads with unique ID and read_id
  onlyPD <- aln[PD, c("seqnames", "read_id")]
  onlyPD <- onlyPD[!duplicated(onlyPD), ]

  # apply filter
  PD <- which(aln$seqnames %in% onlyPD$seqnames)
  read_ids <- onlyPD$read_id[match(aln$seqnames[PD], onlyPD$seqnames)]
  PD <- PD[aln$read_id[PD] == read_ids]
  aln <- aln[-PD,]

  # alignment events filter
  for (i in seq_len(dim(cfgT)[1])) {
    aln_id <- aln[seqnames == cfgT$ID[i], ]
    aln_id <- aln_id[, list(events = (.N/length(unique(strand)))/max(end),
                            counts = max(counts)), by = "read_id"]
    threshS <- thresholdNEvents(aln_id$events)
    onlyBR <- aln_id[aln_id$events > threshS, ]
    onlyBR <- unique(onlyBR, by = "read_id")
    aln <- aln[!(aln$seqnames == cfgT$ID[i] &
                   aln$read_id %in% onlyBR$read_id), ]
  }

  aln
}
