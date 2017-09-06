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
amplicanConsensus <- function(aln, cfgT, overlaps = "overlaps") {

  cols <- c("seqnames", "read_id", "start", "end")
  cols_all <- c("strand", "score", "counts", "width", "num", "originally",
                "replacement", "type", overlaps, cols)
  data.table::setDT(aln)

  aln <- aln[, which(colnames(aln) %in% cols_all), with = FALSE]
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

  consensus[aln_fwd$num[f_both]] <- TRUE # set fwd as true
  # filter these events from further calculations & leave only overlaps
  aln_fwd <- aln_fwd[!f_both & aln_fwd$`overlaps`]
  aln_rve <- aln_rve[!r_both & aln_rve$`overlaps`]
  # The last two columns should be the interval columns
  # find events that are overlapping each other
  data.table::setcolorder(aln_fwd, cols_all)
  data.table::setcolorder(aln_rve, cols_all)
  data.table::setkeyv(aln_fwd, cols)
  data.table::setkeyv(aln_rve, cols)
  oMatch <- data.table::foverlaps(aln_fwd, aln_rve,
                                  type = "any", which = TRUE,
                                  mult = "all", nomatch = 0)
  oScore <- aln_fwd$score[oMatch$xid] >= aln_rve$score[oMatch$yid]
  oScore_fwd <- unique(oMatch$xid[oScore])
  oScore_rve <- unique(oMatch$yid[!oScore])
  consensus[aln_fwd$num[oScore_fwd]] <- TRUE
  consensus[aln_rve$num[oScore_rve]] <- TRUE
  # filter scored events from further calculation
  aln_fwd <- aln_fwd[-oScore_fwd]
  aln_rve <- aln_rve[-oScore_rve]

  # find events that overlap EOP from other strand and set them to true
  data.table::setcolorder(eop_fwd, cols_all)
  data.table::setcolorder(eop_rve, cols_all)
  data.table::setkeyv(eop_fwd, cols)
  data.table::setkeyv(eop_rve, cols)
  oMatch <- data.table::foverlaps(aln_fwd, eop_rve,
                                  type = "any", which = TRUE,
                                  mult = "all", nomatch = 0)
  consensus[aln_fwd$num[unique(oMatch$xid)]] <- TRUE
  data.table::setkeyv(aln_rve, cols)
  oMatch <- data.table::foverlaps(aln_rve, eop_fwd,
                                  type = "any", which = TRUE,
                                  mult = "all", nomatch = 0)
  consensus[aln_rve$num[unique(oMatch$xid)]] <- TRUE

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
amplicanOverlap <- function(aln, cfgT, cut_buffer = 5) {
  cutSites <- lapply(cfgT$ID, function(x) {
    upperGroups(get_amplicon(cfgT, x)) + cut_buffer})
  cutSitesCheck <- sapply(cutSites, length) == 0
  if (any(cutSitesCheck)) {
    message("Warning: Config file row without upper case groups (guideRNA): ",
            toString(which(cutSitesCheck)))
    cutSites[cutSitesCheck] <- as.list(IRanges::tile(
      IRanges::IRanges(start = 1,
                       width = cfgT$ampl_len[cutSitesCheck]),
      1))
  }
  if (any(aln$start <= 0 | aln$end <= 0)) {
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
  aln <- aln[type != "mismatch"] # mismatch has width of 1

  cfgT$Reads_Del <- 0
  cfgT$Reads_In <- 0
  cfgT$Reads_Indel <- 0
  cfgT$Reads_Frameshifted <- 0

  # Frameshift
  aln[type == "deletion", width := width * -1L] # by reference
  widthT <- aln[, list(width = sum(width)),
                by = c("seqnames", "read_id", "counts")]
  widthT <- widthT[widthT$width %% 3 != 0, ]
  widthT_final <- widthT[, list(counts = sum(counts)), by = c("seqnames")]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Frameshifted[map] <- widthT_final$counts

  # Reads that had deletions or insertions
  aln_noD <- unique(aln, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Indel[map] <- widthT_final$counts

  # Reads that had deletion
  alnD <- aln[aln$type == "deletion", ]
  aln_noD <- unique(alnD, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_Del[map] <- widthT_final$counts

  # Reads that had insertion
  aln <- aln[aln$type == "insertion", ]
  aln_noD <- unique(aln, by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$Reads_In[map] <- widthT_final$counts

  cfgT
}
