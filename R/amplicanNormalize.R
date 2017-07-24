#' Remove events that can be found in Controls.
#'
#' This function can adjust events for small differences between
#' known annotations (amplicon sequences) and real DNA of the strain that
#' was sequenced.
#' @param aln (data.frame) Contains events from alignments.
#' @param cfgT (data.frame) Config table with information about experiments.
#' @param add (character vector) Columns from cfgT that should be included
#' in event table for normaliztion matching. Defaults to c("guideRNA", "Group"),
#' which means that only those events created by the same guideRNA in the same
#' Group will be removed if found in Control.
#' @param skip (character vector) Specifies which column of aln to skip,
#' defaults to c('counts', 'score').
#' @return (data.frame) Same as aln, but events are normalized. Events from
#' Control are not changed. Additionaly columns from add are added to the
#' data.frame.
#' @family analysis steps
#' @export
#' @examples
#' aln <- data.frame(seqnames = 1:5,
#'                   start = 1,
#'                   end = 2,
#'                   width = 2,
#'                   counts = 101:105)
#' cfgT <- data.frame(ID = 1:5,
#'                    guideRNA = rep("ACTG", 5),
#'                    Group = c("A", "A", "B", "B", "B"),
#'                    Control = c(TRUE, FALSE, TRUE, FALSE, FALSE))
#'# all events are same as in the control group, therfore are filtered out
#'# events from control groups stay
#' amplicanNormalize(aln, cfgT)
#'# events that are different from control group are preserved
#' aln[2, "start"] <- 3
#' amplicanNormalize(aln, cfgT)
#'
amplicanNormalize <- function(aln, cfgT,
                              add = c("guideRNA", "Group"),
                              skip = c("counts", "score")){
  if (!any(cfgT$Control)) {
    warning("Column 'Control' has no TRUE/1 values. Nothing to normalize.")
    return(aln)
  }

  map <- match(aln$seqnames, cfgT$ID)
  for (column in add) {
    aln[[column]] <- cfgT[[column]][map]
  }
  aln[["Control"]] <- cfgT[["Control"]][map]
  aln_ctr <- aln[aln$Control, ]
  aln <- aln[!aln$Control, ]
  data.table::setDT(aln)
  data.table::setDT(aln_ctr)

  cols <- names(aln)[!names(aln) %in% c(skip, "Control", "seqnames", "read_id")]
  aln <- aln[!aln_ctr, on = cols]

  aln <- data.table::rbindlist(list(aln, aln_ctr))
  aln <- aln[, colnames(aln)[!colnames(aln) %in% c(add, "Control")], with=FALSE]
  data.table::setDF(aln)
  aln
}
