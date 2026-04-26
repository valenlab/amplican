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

  eOP <- findEOP(aln, cfgT)
  aln <- aln[!eOP, ]

  PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

  # PRIMER DIMER reads with unique ID and read_id
  onlyPD <- unique(aln[PD, .(seqnames, read_id)])

  # Native data.table anti-join
  aln <- aln[!onlyPD, on = .(seqnames, read_id)]

  # alignment events filter
  # Using data.table fast binary search instead of O(N) full vector scan
  bad_reads_list <- lapply(seq_len(nrow(cfgT)), function(i) {
    aln_id <- aln[.(cfgT$ID[i]), on = "seqnames", nomatch = NULL]
    if (nrow(aln_id) == 0) return(NULL)
    onlyBR <- aln_id[findLQR(aln_id), ]
    onlyBR <- unique(onlyBR, by = "read_id")
    if (nrow(onlyBR) > 0) return(onlyBR[, c("seqnames", "read_id"), with = FALSE])
    return(NULL)
  })

  bad_reads <- data.table::rbindlist(bad_reads_list)
  if (nrow(bad_reads) > 0) {
    aln <- aln[!bad_reads, on = c("seqnames", "read_id")]
  }

  aln
}
