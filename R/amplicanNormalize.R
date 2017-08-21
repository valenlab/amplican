#' Remove events that can be found in Controls.
#'
#' This function can adjust events for small differences between
#' known annotations (amplicon sequences) and real DNA of the strain that
#' was sequenced. Events from the control are grouped by \code{add} and
#' their frequencies are calculated in respect to number of total reads in that
#' groups. In next step events from the control are filtered according to
#' \code{min_freq}, all events below are treated as sequencing errors and
#' rejected. Finally, all events that can be found in treatment group that find
#' their exact match in control group are removed. All events from control group
#' are returned back.
#' @param aln (data.frame) Contains events from alignments.
#' @param cfgT (data.frame) Config table with information about experiments.
#' @param add (character vector) Columns from cfgT that should be included
#' in event table for normalization matching. Defaults to c("guideRNA", "Group")
#' , which means that only those events created by the same guideRNA in the same
#' Group will be removed if found in Control.
#' @param skip (character vector) Specifies which column of aln to skip,
#' defaults to c('counts', 'score').
#' @param min_freq (numeric) All events from control group below this frequency
#' will be not included in filtering. Use this to filter out background noise
#' and sequencing errors.
#' @return (data.frame) Same as aln, but events are normalized. Events from
#' Control are not changed. Additionally columns from add are added to the
#' data.frame.
#' @family analysis steps
#' @export
#' @examples
#' aln <- data.frame(seqnames = 1:5, start = 1, end = 2, width = 2,
#'                   counts = 101:105)
#' cfgT <- data.frame(ID = 1:5, guideRNA = rep("ACTG", 5),
#'                    Reads_noPD = c(2, 2, 3, 3, 4),
#'                    Group = c("A", "A", "B", "B", "B"),
#'                    Control = c(TRUE, FALSE, TRUE, FALSE, FALSE))
#' # all events are same as in the control group, therefore are filtered out
#' # events from control groups stay
#' amplicanNormalize(aln, cfgT)
#' # events that are different from control group are preserved
#' aln[2, "start"] <- 3
#' amplicanNormalize(aln, cfgT)
#'
amplicanNormalize <- function(aln, cfgT,
                              add = c("guideRNA", "Group"),
                              skip = c("counts", "score"),
                              min_freq = 0.1){
  Reads_noPD <- frequency <- NULL
  if (!any(cfgT$Control)) {
    warning("Column 'Control' has no TRUE/1 values. Nothing to normalize.")
    return(aln)
  }
  data.table::setDT(aln)
  map <- match(aln$seqnames, cfgT$ID)
  for (column in add) {
    aln[[column]] <- cfgT[[column]][map]
  }
  cols <- names(aln)[!names(aln) %in% c(skip, "seqnames", "read_id")]

  ctr_indices <- cfgT[["Control"]][map]
  aln_ctr <- aln[ctr_indices, ]
  data.table::setkeyv(aln_ctr, cols)
  aln <- aln[!ctr_indices]
  data.table::setkeyv(aln, cols)

  data.table::setDT(cfgT)
  cfgT_total_reads <- cfgT[, list(Reads_noPD = sum(Reads_noPD)),
                          by = c(add, "Control")]
  cfgT_total_reads <- cfgT_total_reads[cfgT_total_reads$Control, ]

  aln_ctr_freq <- aln_ctr[, list(counts = sum(counts)), by = cols]
  aln_ctr_freq <-  merge(aln_ctr_freq, cfgT_total_reads, all = TRUE, by = add)
  aln_ctr_freq$frequency <- aln_ctr_freq$counts/aln_ctr_freq$Reads_noPD
  aln_ctr_freq <- aln_ctr_freq[frequency > min_freq, ]

  # dplyr used as data.table has issues handling too big dt
  aln <- dplyr::anti_join(aln, aln_ctr_freq, by = cols)
  aln <- data.table::rbindlist(list(aln, aln_ctr))
  aln <- aln[, colnames(aln)[!colnames(aln) %in% add], with=FALSE]
  data.table::setDF(aln)
  aln
}
