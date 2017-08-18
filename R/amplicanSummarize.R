#' Extract consensus out of forward and reverse events.
#'
#' When forward and reverse reads are in agreement on the events (eg. deletion)
#' \code{amplicanConsensus} will mark forward event as TRUE indicating that he
#' represents consensus.
#' In cases where forward and reverse read agree only partially, for example,
#' they share the same start of the deletion, but they have different end
#' \code{amplicanConsensus} will pick the version of
#' read with higher alignment score, in situation where both of the reads
#' overlap expected cut site, otherwise both events will be rejected and marked
#' FALSE. When there are events only on one of the strands they will be
#' rejected.
#'
#' In situation where you have only forward or only reverse reads don't use this
#' function and assign all TRUE to all of your events.
#'
#' Consensus out of the forward + reverse reads is required for
#' \code{amplicanSummary}, and \code{amplicanConsensus} requires
#' \code{amplicanOverlap}.
#'
#' @param aln (data.frame) Contains relevant events in GRanges style.
#' @param overlaps (character) Specifies which metadata column of \code{aln}
#' indicates which events are overlapping expected cut site.
#' @return (bolean vector) Where TRUE means that given event represents
#' consensus out of forward and reverse reads.
#' @export
#' @include helpers_general.R
#' @family analysis steps
#'
amplicanConsensus <- function(aln, overlaps = "overlaps") {

  data.table::setDT(aln)
  consensus <- rep(FALSE, nrow(aln))

  aln$num <- seq_len(nrow(aln))
  aln_fwd <- aln[strand == "+"]
  aln_rve <- aln[strand == "-"]
  cols <- c("seqnames", "read_id", "start", "end", "replacement", "type")
  data.table::setkeyv(aln_fwd, cols)
  data.table::setkeyv(aln_rve, cols)

  # find those events that are confirmed by both fwd and rve
  f_both <- !is.na(aln_rve[aln_fwd, which = TRUE, mult = "first"])
  r_both <- !is.na(aln_fwd[aln_rve, which = TRUE, mult = "first"])

  consensus[aln_fwd$num[f_both]] <- TRUE # set fwd as true
  # filter these events from further calculations
  aln_fwd <- aln_fwd[!f_both & aln_fwd$`overlaps`]
  aln_rve <- aln_rve[!r_both & aln_rve$`overlaps`]
  # The last two columns should be the interval columns
  # find events that can are overlapping each other
  cols <- c("seqnames", "read_id", "strand", "score", "counts", "width",
            "type", "num", "originally","replacement", overlaps, "start", "end")
  data.table::setcolorder(aln_fwd, cols)
  data.table::setcolorder(aln_rve, cols)
  cols <- c("seqnames", "read_id", "type", "replacement", "start", "end")
  data.table::setkeyv(aln_fwd, cols)
  data.table::setkeyv(aln_rve, cols)
  oMatch <- data.table::foverlaps(aln_fwd, aln_rve,
                                  type = "any", which = TRUE,
                                  mult = "all", nomatch = 0)
  oScore <- aln_fwd$score[oMatch$xid] >= aln_rve$score[oMatch$yid]
  consensus[aln_fwd$num[oMatch$xid[oScore]]] <- TRUE
  consensus[aln_rve$num[oMatch$yid[!oScore]]] <- TRUE

  return(consensus)
}


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
#' relative to their cut sites, and most LEFT upper case letter specifies
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
  # overlap assessment
  for(i in seq_along(cfgT$ID)) {
    map <- which(aln$seqnames == cfgT$ID[i])
    overlap[map] <- IRanges::overlapsAny(alnIR[map], cutSites[[i]])
  }
  overlap
}


#' Summarize how many reads have frameshift and how many reads have deletions.
#'
#' Before using this function make sure events are filtered to represent
#' consensus with \code{amplicanConsensus}, if you use both forward and
#' reverse reads. If you want to calculate metrics over expected cut site,
#' filter events using \code{amplicanOverlap}.
#'
#' Adds columns to cfgT:
#' \itemize{
#' \item{ReadsCut}{ Count of reads with deletions overlapping expected
#' cut site.}
#' \item{Reads_Frameshifted}{ Count of reads with frameshift
#' overlapping expected cut site.}
#' }
#' @param aln (data.frame) Contains events from the alignments.
#' @param cfgT (data.frame) Config file with the experiments details.
#' @return (data.frame) As cfgT, but with extra columns.
#' @export
#' @family analysis steps
#' @include helpers_general.R
#'
amplicanSummarize <- function(aln, cfgT) {

  seqnames <- read_id <- counts <- start <- end <- score <- NULL
  data.table::setDT(aln)
  aln <- aln[type != "mismatch"] # mismatch has width of 1

  cfgT$Reads_Cut <- 0
  cfgT$Reads_Frameshifted <- 0

  # Frameshift
  aln[type == "deletion", width := width * -1L] # by reference
  widthT <- aln[, list(width = sum(width)),
                by = c("seqnames", "read_id", "counts")]
  widthT <- widthT[widthT$width %% 3 != 0, ]
  widthT_final <- widthT[, list(counts = sum(counts)), by = c("seqnames")]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Frameshifted[map] <- widthT_final$counts

  # Cut rate
  aln <- aln[aln$type == "deletion", ]
  aln_noD <- unique(aln, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Cut[map] <- widthT_final$counts

  cfgT
}
