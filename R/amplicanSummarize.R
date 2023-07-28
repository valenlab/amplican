getHits <- function(aln_fwd, aln_rve) {
  if (nrow(aln_fwd) == 0 | nrow(aln_rve) == 0) return(S4Vectors::Hits())
  suppressWarnings(GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(
      seqnames = paste0(aln_fwd$seqnames, "_", aln_fwd$read_id),
      ranges = IRanges::IRanges(start = aln_fwd$start,
                                end = aln_fwd$end),
      strand = "*"),
    GenomicRanges::GRanges(
      seqnames = paste0(aln_rve$seqnames, "_", aln_rve$read_id),
      IRanges::IRanges(start = aln_rve$start,
                       end = aln_rve$end),
      strand = "*"),
    type = "any", select = "all"))
}


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
#' @param cfgT (data.frame) Should be table containing at least
#' positions of primers in the amplicons and their identifiers
#' @param overlaps (character) Specifies which metadata column of \code{aln}
#' indicates which events are overlapping expected cut site.
#' @param promiscuous (boolean) Allows to relax consensus rules. When TRUE will
#' allow Indels that are not confirmed by the other strand (when both are used).
#' @return (bolean vector) Where TRUE means that given event represents
#' consensus out of forward and reverse reads.
#' @export
#' @include helpers_general.R
#' @family analysis steps
#' @examples
#' file_path <- system.file("test_data", "test_aln.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("test_data", "test_cfg.csv", package = "amplican"))
#' all(aln$consensus == amplicanConsensus(aln, cfgT))
#'
amplicanConsensus <- function(aln, cfgT, overlaps = "overlaps",
                              promiscuous = TRUE) {

  cols <- c("seqnames", "read_id", "start", "end")
  cols_all <- c("strand", "score", "counts", "width", "num", "originally",
                "replacement", "type", overlaps, cols)
  data.table::setDT(aln)

  aln <- aln[, which(colnames(aln) %in% cols_all), with = FALSE]
  if (dim(aln)[1] == 0) return(logical(0))
  consensus <- rep(FALSE, nrow(aln))
  aln$num <- seq_len(nrow(aln))

  # find EOP if any
  eop <- findEOP(aln, cfgT)
  eop_aln <- aln[eop & type == "deletion"]
  aln <- aln[!eop]

  eop_fwd <- eop_aln[strand == "+"]
  eop_rve <- eop_aln[strand == "-"]
  aln_fwd <- aln[strand == "+"]
  aln_rve <- aln[strand == "-"]

  data.table::setkeyv(aln_fwd, c(cols, "replacement", "type"))
  data.table::setkeyv(aln_rve, c(cols, "replacement", "type"))

  # find those events that are confirmed by both fwd and rve
  f_both <- !is.na(aln_rve[aln_fwd, which = TRUE, mult = "first"])
  r_both <- !is.na(aln_fwd[aln_rve, which = TRUE, mult = "first"])
  consensus[aln_fwd$num[f_both]] <- TRUE

  # filter these events from further calculations & leave only overlaps
  aln_fwd <- aln_fwd[!f_both & aln_fwd$`overlaps`]
  aln_rve <- aln_rve[!r_both & aln_rve$`overlaps`]

  # reads that have eop & overlaps true should take info from the other strand
  # unless the other strand is also broken
  b_rve <- eop_rve[eop_rve$overlaps, ]
  b_fwd <- eop_fwd[eop_fwd$overlaps, ]
  if (nrow(b_rve) + nrow(b_fwd) > 0) {
    # both strands are broken
    b_rve$seq_id <- paste0(b_rve$seqnames, "*_*", b_rve$read_id)
    b_fwd$seq_id <- paste0(b_fwd$seqnames, "*_*", b_fwd$read_id)
    b_rve <- b_rve[!b_rve$seq_id %in% b_fwd$seq_id, ]
    b_fwd <- b_fwd[!b_fwd$seq_id %in% b_rve$seq_id, ]
    # filter out events from those broken IDs
    aln_rve_seq_id <- paste0(aln_rve$seqnames, "*_*", aln_rve$read_id)
    aln_fwd_seq_id <- paste0(aln_fwd$seqnames, "*_*", aln_fwd$read_id)
    aln_rve <- aln_rve[!aln_rve_seq_id %in% b_rve$seq_id, ]
    aln_fwd <- aln_fwd[!aln_fwd_seq_id %in% b_fwd$seq_id, ]
  }

  # The last two columns should be the interval columns
  # find events that are overlapping each other

  for (useq in unique(aln_fwd$seqnames)) {
    oMatch <- getHits(aln_fwd[aln_fwd$seqnames == useq],
                      aln_rve[aln_rve$seqnames == useq])
    fi <- S4Vectors::from(oMatch)
    ri <- S4Vectors::to(oMatch)
    oScore <- aln_fwd[aln_fwd$seqnames == useq]$score[fi] >=
      aln_rve[aln_rve$seqnames == useq]$score[ri]

    oScore_fwd <- unique(fi[oScore])
    oScore_rve_not <- unique(ri[oScore])
    oScore_rve <- unique(ri[!oScore])
    oScore_fwd_not <- unique(fi[!oScore])
    consensus[aln_fwd[aln_fwd$seqnames == useq]$num[oScore_fwd]] <- TRUE
    consensus[aln_rve[aln_rve$seqnames == useq]$num[oScore_rve]] <- TRUE

    # filter scored events from further calculation
    if (length(c(oScore_fwd, oScore_fwd_not)) > 0) {
      aln_fwd <- rbind(aln_fwd[aln_fwd$seqnames == useq][-c(oScore_fwd, oScore_fwd_not), ],
                       aln_fwd[aln_fwd$seqnames != useq])
    }
    if (length(c(oScore_rve, oScore_rve_not)) > 0) {
      aln_rve <- rbind(aln_rve[aln_rve$seqnames == useq][-c(oScore_rve, oScore_rve_not), ],
                       aln_rve[aln_rve$seqnames != useq])
    }
  }



  if (!promiscuous) {
    # find events that overlap EOP from other strand and set them to true
    oMatch <- getHits(aln_fwd, eop_rve)
    consensus[aln_fwd$num[unique(S4Vectors::from(oMatch))]] <- TRUE
    oMatch <- getHits(aln_rve, eop_fwd)
    consensus[aln_rve$num[unique(S4Vectors::from(oMatch))]] <- TRUE
  } else { # not strict
    # all events that are left, don't overlap each other
    consensus[aln_fwd$num] <- TRUE
    consensus[aln_rve$num] <- TRUE
  }
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
#' @param relative (boolean) Sets whether events are relative to the position of
#' the target site.
#' @return (bolean vector) Where TRUE means that given event overlaps cut site.
#' @export
#' @include helpers_general.R
#' @family analysis steps
#' @examples
#' file_path <- system.file("test_data", "test_aln.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("test_data", "test_cfg.csv", package = "amplican"))
#' all(aln$overlaps == amplicanOverlap(aln, cfgT))
#'
amplicanOverlap <- function(aln, cfgT, cut_buffer = 5, relative = FALSE) {
  if (dim(aln)[1] == 0) return(logical(0))
  cutSites <- lapply(cfgT$ID, function(x) {
    upperGroups(get_seq(cfgT, x)) + cut_buffer})
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
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "events_filtered_shifted_normalized.csv",
#'                          package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' amplicanSummarize(aln, cfgT)
#'
amplicanSummarize <- function(aln, cfgT) {
  seqnames <- read_id <- counts <- start <- end <- score <- NULL
  data.table::setDT(aln)
  cfgT$HDR <- cfgT$Reads_Del <- cfgT$Reads_In <-
    cfgT$Reads_Edited <- cfgT$Reads_Frameshifted <- 0

  # HDR
  aln_noD <- unique(aln[aln$readType, ], by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$HDR[map] <- widthT_final$counts

  # Reads that had deletion
  alnD <- aln[aln$type == "deletion", ]
  aln_noD <- unique(alnD, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Del[map] <- widthT_final$counts

  # Reads that had insertion
  alnI <- aln[aln$type == "insertion", ]
  aln_noD <- unique(alnI, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_In[map] <- widthT_final$counts

  # Reads that had editing: deletions or insertions (or HDR)
  aln_E <- aln[aln$type %in% c("insertion", "deletion") | aln$readType, ]
  aln_noD <- unique(aln_E, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Edited[map] <- widthT_final$counts

  # Frameshift
  aln <- aln[type != "mismatch"] # mismatch has width of 1
  aln[type == "deletion", width := width * -1L] # by reference
  widthT <- aln[, list(width = sum(width)),
                by = c("seqnames", "read_id", "counts")]
  widthT <- widthT[widthT$width %% 3 != 0, ]
  widthT_final <- widthT[, list(counts = sum(counts)), by = c("seqnames")]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Frameshifted[map] <- widthT_final$counts
  cfgT
}
