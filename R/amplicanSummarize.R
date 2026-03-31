getHits <- function(aln_fwd, aln_rve) {
  if (nrow(aln_fwd) == 0 | nrow(aln_rve) == 0) return(S4Vectors::Hits())

  fwd <- aln_fwd[, .(seqnames, read_id, start, end)]
  fwd[, fwd_idx := .I]
  rve <- aln_rve[, .(seqnames, read_id, start, end)]
  rve[, rve_idx := .I]

  data.table::setkey(rve, seqnames, read_id, start, end)
  hits <- data.table::foverlaps(fwd, rve, type = "any", nomatch = NULL)

  if (nrow(hits) == 0) return(S4Vectors::Hits())
  S4Vectors::Hits(from = hits$fwd_idx, to = hits$rve_idx,
                  nLnode = nrow(aln_fwd), nRnode = nrow(aln_rve),
                  sort.by.query = FALSE)
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
  if (nrow(aln) == 0) return(logical(0))
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
    # both strands are broken — remove reads that appear in the other strand
    b_rve <- b_rve[!b_fwd, on = .(seqnames, read_id)]
    b_fwd <- b_fwd[!b_rve, on = .(seqnames, read_id)]
    # filter out events from those broken IDs
    aln_rve <- aln_rve[!b_rve, on = .(seqnames, read_id)]
    aln_fwd <- aln_fwd[!b_fwd, on = .(seqnames, read_id)]
  }

  # Single vectorized getHits run avoiding nested seqnames iteration loop matching
  fwd_to_remove <- integer()
  rve_to_remove <- integer()

  oMatch <- getHits(aln_fwd, aln_rve)

  if (length(oMatch) > 0) {
    fi <- S4Vectors::from(oMatch)
    ri <- S4Vectors::to(oMatch)
    oScore <- aln_fwd$score[fi] >= aln_rve$score[ri]

    oScore_fwd <- unique(fi[oScore])
    oScore_rve_not <- unique(ri[oScore])
    oScore_rve <- unique(ri[!oScore])
    oScore_fwd_not <- unique(fi[!oScore])

    consensus[aln_fwd$num[oScore_fwd]] <- TRUE
    consensus[aln_rve$num[oScore_rve]] <- TRUE

    # filter scored events from further calculation
    fwd_to_remove <- unique(c(oScore_fwd, oScore_fwd_not))
    rve_to_remove <- unique(c(oScore_rve, oScore_rve_not))
  }

  if (length(fwd_to_remove) > 0) aln_fwd <- aln_fwd[-fwd_to_remove, ]
  if (length(rve_to_remove) > 0) aln_rve <- aln_rve[-rve_to_remove, ]



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
  if (nrow(aln) == 0) return(logical(0))
  cutSites <- lapply(seq_along(cfgT$ID), function(i) {
    upperGroups(get_seq(cfgT, cfgT$ID[i], row = i)) + cut_buffer})
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
  overlap <- logical(nrow(aln))
  # map each event to its experiment row index once
  map <- match(aln$seqnames, cfgT$ID)
  idx_by_exp <- split(seq_len(nrow(aln)), map)
  for (exp_i in names(idx_by_exp)) {
    rows <- idx_by_exp[[exp_i]]
    overlap[rows] <- IRanges::overlapsAny(alnIR[rows], cutSites[[as.integer(exp_i)]])
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
#' \describe{
#' \item{HDR}{Count of reads identified as Homology Directed Repair events.}
#' \item{Reads_Del}{Count of reads containing at least one deletion.}
#' \item{Reads_In}{Count of reads containing at least one insertion.}
#' \item{Reads_Edited}{Count of reads with any edit (insertion, deletion, or HDR).}
#' \item{Reads_Frameshifted}{Count of reads with a frameshift (net indel length is not a multiple of 3).}
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
  seqnames <- read_id <- counts <- type <- readType <- width <- NULL
  has_HDR <- has_Del <- has_In <- has_Edit <- is_FS <- NULL
  i.HDR <- i.Reads_Del <- i.Reads_In <- i.Reads_Edited <- i.Reads_Frameshifted <- NULL
  
  data.table::setDT(aln)
  data.table::setDT(cfgT)
  
  # Vectorize net_width globally (avoids expensive ifelse inside grouping)
  aln[, net_width := 0L]
  aln[type == "deletion", net_width := -width]
  aln[type == "insertion", net_width := width]

  # Step 1: Summarize what happened in each read_id in one go
  read_summary <- aln[, .(
    has_HDR = any(readType == TRUE),
    has_Del = any(type == "deletion"),
    has_In  = any(type == "insertion"),
    is_FS   = sum(net_width) %% 3 != 0,
    read_counts = max(counts) # Assuming counts are identical for the same read_id
  ), by = .(seqnames, read_id)]

  aln[, net_width := NULL]  # clean up
  
  read_summary[, has_Edit := has_Del | has_In | has_HDR]
  
  # Step 2: Aggregate up to the seqnames (experiment) level
  exp_summary <- read_summary[, .(
    HDR = sum(read_counts * has_HDR),
    Reads_Del = sum(read_counts * has_Del),
    Reads_In = sum(read_counts * has_In),
    Reads_Edited = sum(read_counts * has_Edit),
    Reads_Frameshifted = sum(read_counts * is_FS)
  ), by = seqnames]
  
  # Step 3: Update cfgT by reference
  cfgT[exp_summary, `:=`(
    HDR = i.HDR,
    Reads_Del = i.Reads_Del,
    Reads_In = i.Reads_In,
    Reads_Edited = i.Reads_Edited,
    Reads_Frameshifted = i.Reads_Frameshifted
  ), on = .(ID = seqnames)]
  
  # Fill NAs with 0 for experiments that had no events
  cols <- c("HDR", "Reads_Del", "Reads_In", "Reads_Edited", "Reads_Frameshifted")
  for (j in cols) data.table::set(cfgT, which(is.na(cfgT[[j]])), j, 0)
  
  return(data.table::setDF(cfgT))
}
