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
#' when seqnames not presnt it means no frameshift in any of the reads for this
#' experiment. Normally counts for + and - are always the same.
#' @param aln (data.frame)
#' @param paired_end (boolean)
#' @return (data.frame)
#'
# test case:
# aln <- data.table::data.table(seqnames = rep(1:3, each = 2),
#                               read_id = rep(1:3, each = 2),
#                               strand = rep(c("+", "-"), times = 3),
#                               score = c(10, 10, 15, 10, 10, 15),
#                               counts = c(1:6),
#                               width = 1:6)
# aln <- rbind(aln, list(15, 1, "+", 10, 10, 10))
# aln <- rbind(aln, list(16, 1, "+", 10, 100, 10),
#                   list(16, 1, "-", 10, 100, 10))
# aln$type <- "deletion"
frameshifted_reads_by_ID <- function(aln, paired_end = TRUE) {
  seqnames <- read_id <- NULL

  aln <- aln[type != "mismatch"]
  aln[type == "deletion", width := width * -1L] # by reference

  if (!paired_end) {
    widthT <- aln[, list(width = sum(width)),
                  by = c("seqnames", "read_id", "counts")]
    widthT <- widthT[widthT$width %% 3 != 0, ]
    widthT_final <- widthT[, list(counts = sum(counts)), by = c("seqnames")]
  } else {
    widthT <- aln[, list(width = sum(width)),
                  by = c("seqnames", "read_id", "strand", "score", "counts")]
    widthT_fwd <- widthT[strand == "+",]
    widthT_rve <- widthT[strand == "-",]
    widthT_col <- names(widthT)
    # fwd that are not in rve - disregarding width and score
    fwd_not_rve <- widthT_fwd[!widthT_rve, on = list(seqnames, read_id)]
    rve_not_fwd <- widthT_rve[!widthT_fwd, on = list(seqnames, read_id)]
    # fwd and rve agree on width
    fwd_as_rve <- widthT_fwd[widthT_rve, .SD, on = list(width,
                                                        seqnames,
                                                        read_id),
                             nomatch = 0L, .SDcols = widthT_col]
    # drop out of fwd and rve rows that have equal width
    widthT_fwd2 <- widthT_fwd[!widthT_rve, on = list(width, seqnames, read_id)]
    widthT_rve2 <- widthT_rve[!widthT_fwd, on = list(width, seqnames, read_id)]

    # fwd that have score equal as rve and different widths
    # theoretical case - highly unlikely take forward
    fwd_eq_rve <- if (nrow(widthT_fwd2) > 0) {
      widthT_fwd2[widthT_rve2, .SD, on = list(score == score,
                                              seqnames, read_id),
                  nomatch = 0L, .SDcols = widthT_col]
    } else data.table::data.table()

    # fwd that have score higher than rve score and different widths
    fwd_big_rve <- if (nrow(widthT_fwd2) > 0) {
      widthT_fwd2[widthT_rve2, .SD, on = list(score > score,
                                              seqnames, read_id),
                  nomatch = 0L, .SDcols = widthT_col]
    } else data.table::data.table()

    # rve that have score higher than fwd score and different widths
    rve_big_fwd <- if (nrow(widthT_rve2) > 0) {
      widthT_rve2[widthT_fwd2, .SD, on = list(score > score,
                                              seqnames, read_id),
                  nomatch = 0L, .SDcols = widthT_col]
    } else data.table::data.table()

    # combine all reduced reads together
    widthT <- data.table::rbindlist(
      list(fwd_as_rve, fwd_big_rve, rve_big_fwd,
           fwd_not_rve, rve_not_fwd, fwd_eq_rve))

    widthT <- widthT[width %% 3 != 0]
    widthT_final <- widthT[, list(counts = sum(counts)), by = c("seqnames")]
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
  seqnames <- read_id <- counts <- start <- end <- score <- NULL
  # rewrite to make use of data.table
  paired_end <- length(unique(aln$strand)) > 1
  data.table::setDT(aln)
  aln <- aln[type != "mismatch"]

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
    aln_noD <- unique(aln, by = c("seqnames", "read_id"))
    widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Cut[map] <- widthT_final$counts
  } else {
    # Frameshift
    widthT_final <- frameshifted_reads_by_ID(aln, TRUE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted[map] <- widthT_final$counts

    aln <- aln[aln[, overlaps], ]

    # Frameshfit Overlaping
    widthT_final <- frameshifted_reads_by_ID(aln, TRUE)
    map <- match(widthT_final$seqnames, cfgT$ID)
    cfgT$Reads_Frameshifted_Overlapping[map] <- widthT_final$counts

    # Cut rate
    aln <- aln[aln$type == "deletion", ]
    aln_fwd <- aln[strand == "+",]
    aln_rve <- aln[strand == "-",]
    aln_col <- names(aln)
    # fwd that are not in rve - disregarding start, end and score
    fwd_not_rve <- aln_fwd[!aln_rve, on = list(seqnames, read_id)]
    rve_not_fwd <- aln_rve[!aln_fwd, on = list(seqnames, read_id)]
    # fwd and rve agree on start and end
    fwd_as_rve <- aln_fwd[aln_rve, .SD, on = list(start, end,
                                                  seqnames, read_id),
                          nomatch = 0L, .SDcols = aln_col]
    # drop out of fwd and rve rows that have equal width
    aln_fwd2 <- aln_fwd[!aln_rve, on = list(start, end, seqnames, read_id)]
    aln_rve2 <- aln_rve[!aln_fwd, on = list(start, end, seqnames, read_id)]

    # fwd that have score equal as rve and different widths
    # theoretical case - highly unlikely take forward
    fwd_eq_rve <- if (nrow(aln_fwd2) > 0) {
      aln_fwd2[aln_rve2, .SD, on = list(start == start, end == end,
                                        seqnames, read_id),
               nomatch = 0L, .SDcols = aln_col]
    } else data.table::data.table()

    # fwd that have score higher than rve score and different start and end
    fwd_big_rve <- if (nrow(aln_fwd2) > 0) {
      aln_fwd2[aln_rve2, .SD, on = list(score > score,
                                        seqnames, read_id),
               nomatch = 0L, .SDcols = aln_col]
    } else data.table::data.table()

    # rve that have score higher than fwd score and different start and end
    rve_big_fwd <- if (nrow(aln_rve2) > 0) {
      aln_rve2[aln_fwd2, .SD, on = list(score > score,
                                        seqnames, read_id),
               nomatch = 0L, .SDcols = aln_col]
    } else data.table::data.table()

    # combine all reduced reads together - all considered cuts
    aln <- data.table::rbindlist(
      list(fwd_as_rve, fwd_big_rve, rve_big_fwd,
           fwd_not_rve, rve_not_fwd, fwd_eq_rve))

    aln <- unique(aln, by = c("seqnames", "read_id"))
    aln <- aln[, list(counts = sum(counts)), by = c("seqnames")]
    map <- match(aln$seqnames, cfgT$ID)
    cfgT$Reads_Cut[map] <- aln$counts
  }

  cfgT
}
