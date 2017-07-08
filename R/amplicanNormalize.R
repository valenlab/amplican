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

  skip <- c(skip, "Control", "seqnames", "read_id")
  skip_col <- which(names(aln) %in% skip)
  skip_col_ctr <- which(names(aln_ctr) %in% skip)
  # restrict number of potential events that repeat - for speed
  candid <- duplicated(rbind(aln[, -skip_col],
                             aln_ctr[, -skip_col_ctr]),
                       fromLast = TRUE)[seq_len(dim(aln)[1])]
  # check which of candidates are really replicates
  temp_aln_ctr <- aln_ctr[, -skip_col_ctr]
  next_row <- dim(temp_aln_ctr)[1] + 1
  to_norm <- apply(aln[candid, -skip_col], 1, function(x) {
    temp_aln_ctr[next_row, ] <- x
    duplicated(temp_aln_ctr, fromLast = FALSE)[next_row]
  })

  candid[candid] <- to_norm
  aln <- aln[!candid, ]

  rbind(aln, aln_ctr)
}
