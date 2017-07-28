#' Check which events overlap expected cut sites.
#'
#' To determine which deletions, insertions and mismatches (events) are probably
#' created by CRISPR we check whether they overlap expected cut sites. Expected
#' cut sites should be specified in UPPER CASE letters in the amplicon
#' sequences.
#' @param aln (data.frame) Contains relevant events in GRanges style.
#' @param cfgT (data.frame) Contains amplicon sequences.
#' @param cut_buffer (numeric) Number of bases that should expand 5' and 3' of
#' the specified expected cut sites.
#' @param relative (boolean) Default is TRUE, which means 'aln' events are
#' relative to thir cut sites, and most LEFT upper case letter specifies
#' position zero.
#' @return (bolean vector) Where TRUE means that given event overlaps cut site.
#' @export
#' @include helpers_general.R
#' @family analysis steps
#'
amplicanOverlap <- function(aln, cfgT, cut_buffer = 5, relative = TRUE) {

  cutSites <- lapply(
    cfgT$Amplicon, function(x) upperGroups(x) + cut_buffer)
  cutSitesCheck <- sapply(cutSites, length) == 0
  if (any(cutSitesCheck)) {
    message("Warning: Config file row without upper case groups (guideRNA): ",
            toString(which(cutSitesCheck)))
    cutSites[cutSitesCheck] <- as.list(IRanges::tile(
      IRanges::IRanges(start = 1,
                       width = cfgT$ampl_len[cutSitesCheck]),
      1))
  }
  if (relative) {
    cutSites <- lapply(cutSites, function(x) {
      IRanges::shift(x, -1 * IRanges::start(x)[1])
    })
  }

  alnIR <- IRanges::IRanges(aln$start, aln$end, aln$width)
  overlap <- vector(length = dim(aln)[1])
  # overlap assesment
  for(i in seq_along(cfgT$ID)) {
    map <- which(aln$seqnames == cfgT$ID[i])
    overlap[map] <- IRanges::overlapsAny(alnIR[map], cutSites[[i]])
  }
  overlap
}


#' Summarize counts of reads with frameshifts for each experiment
#'
#' @param aln (data.frame)
#' @param paired_end (boolean)
#' @return (data.frame)
#'
frameshifted_reads_by_ID <- function(aln, paired_end = TRUE) {
    if (!paired_end) {
      widthT <- stats::aggregate(width ~ seqnames + read_id + strand + counts,
                                 aln, sum)
      widthT <- widthT[widthT$width %% 3 != 0, ]
      widthT_final <- stats::aggregate(counts ~ seqnames, widthT, sum)
    } else {
      widthT <- stats::aggregate(width ~
                                   seqnames + read_id + strand + counts + score,
                                 aln, sum)
      widthT$uniqID <- paste(widthT$seqnames, widthT$read_id, sep = "_")
      widthT_fwd <- widthT[widthT$strand == "+",]
      widthT_rve <- widthT[widthT$strand == "-",]
      # confirmed by the reverse strand
      # trick with duplicates to figure out which fwd and rve are in agreement
      dup <- duplicated(rbind(widthT_fwd[, c("uniqID", "width")],
                              widthT_rve[, c("uniqID", "width")]),
                        fromLast = TRUE)[seq_len(dim(widthT_fwd)[1])]
      # otherwise higher alignment score decides
      for (i in which(!dup)) {
        j <- which(widthT_fwd[i, "uniqID"] == widthT_rve$uniqID)
        if (length(j) > 0 && widthT_fwd[i, "score"] < widthT_rve[j, "score"]) {
          widthT_fwd[i, ] <- widthT_rve[j, ]
          widthT_rve <- widthT_rve[-j, ]
        }
      }
      dup <- duplicated(rbind(widthT_rve[, c("uniqID", "width")],
                              widthT_fwd[, c("uniqID", "width")]),
                        fromLast = TRUE)[seq_len(dim(widthT_rve)[1])]
      widthT_fwd <- rbind(widthT_fwd, widthT_rve[!dup, ])
      widthT_fwd <- widthT_fwd[widthT_fwd$width %% 3 != 0, ]
      widthT_final <- stats::aggregate(counts ~ seqnames, widthT_fwd, sum)
    }
  widthT_final
}


#' Summarize how many reads has frameshift, frameshift overlapping cut site and
#' how many reads with deletions overlaping cut site.
#'
#' Adds columns to cfgT:
#' \itemize{
#' \item{ReadsCut}{ Count of reads with deletions overlapping expected
#' cut site.}
#' \item{Reads_Frameshifted}{ Count of reads with frameshift overlapping
#' expected cut site.}
#' \item{Reads_Frameshifted_Overlapping}{ Count of reads with frameshift
#' overlapping expected cut site.}
#' }
#' @param aln (data.frame) Contains events from the alignments.
#' @param cfgT (data.frame) Config file with the experiments details.
#' @param overlaps (character) Specifies which column of aln contains
#' output from \code{\link{amplicanOverlap}}. Defaults to "overlaps".
#' @return (data.frame) As cfgT, but with extra columns.
#' @export
#' @family analysis steps
#' @include helpers_general.R
#'
amplicanSummarize <- function(aln, cfgT, overlaps = "overlaps") {
  # rewrite to make use of data.table
  paired_end <- length(unique(aln$strand)) > 0

  aln <- aln[aln$type != "mismatch",]

  cfgT$Reads_Cut <- 0
  cfgT$Reads_Frameshifted <- 0
  cfgT$Reads_Frameshifted_Overlapping <- 0

  if (!paired_end) {
    # Frameshift
    widthT_final <- frameshifted_reads_by_ID(aln, FALSE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted[map] <- widthT_final$counts

    aln <- aln[aln[[overlaps]], ]
    # Frameshfit Overlaping
    widthT_final <-frameshifted_reads_by_ID(aln, FALSE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted_Overlapping[map] <- widthT_final$counts

    # Cut rate
    aln <- aln[aln$type == "deletion", ]
    aln_noD <- aln[!duplicated(aln[, c("seqnames", "read_id")]), ]
    widthT_final <- stats::aggregate(counts ~ seqnames, aln, sum)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Cut[map] <- widthT_final$counts
  } else {
    # Frameshift
    widthT_final <- frameshifted_reads_by_ID(aln, TRUE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted[map] <- widthT_final$counts

    aln <- aln[aln[[overlaps]], ]

    # Frameshfit Overlaping
    widthT_final <- frameshifted_reads_by_ID(aln, TRUE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted_Overlapping[map] <- widthT_final$counts

    # Cut rate
    aln <- aln[aln$type == "deletion", ]

    aln$uniqID <- paste(aln$seqnames, aln$read_id, sep = "_")
    aln_fwd <- aln[aln$strand == "+", ]
    aln_rve <- aln[aln$strand == "-", ]
    dup <- duplicated(rbind(aln_fwd[, c("uniqID", "start", "end")],
                            aln_rve[, c("uniqID", "start", "end")]),
                      fromLast = TRUE)[seq_len(dim(aln_fwd)[1])]
    # otherwise higher alignment score decides
    for (i in which(!dup)) {
      j <- which(aln_fwd[i, "uniqID"] == aln_rve$uniqID)
      if (length(j) > 0 && aln_fwd[i, "score"] < aln_rve[j, "score"]) {
        aln_fwd <- aln_fwd[-i, ]
        aln_fwd <- rbind(aln_fwd, aln_rve[j, ])
        aln_rve <- aln_rve[-j, ]
      }
    }
    # add to fwd, rve reads that are only in rve for when only
    # rve has events and fwd has none
    dup <- duplicated(rbind(aln_rve[, c("uniqID", "start", "end")],
                            aln_fwd[, c("uniqID", "start", "end")]),
                      fromLast = TRUE)[seq_len(dim(aln_rve)[1])]
    aln_fwd <- rbind(aln_fwd, aln_rve[!dup, ])
    aln_fwd <- aln_fwd[!duplicated(aln_fwd$uniqID), ]
    widthT_final <- stats::aggregate(counts ~ seqnames, aln_fwd, sum)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Cut[map] <- widthT_final$counts
  }
  cfgT
}
