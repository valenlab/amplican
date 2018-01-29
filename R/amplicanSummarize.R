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

  if (!promiscuous) {
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


#' Realign reads to the donor template, report those read ids that have better
#' score against the donor template and are perfect alignments in regions of the
#' donor template differences toward amplicon sequence.
#'
#' @param reads (PairwiseAlignmentsSingleSubject)
#' @param donor (character)
#' @param hdr_events (IRanges) events from alignment amplican to donor without
#' shifting and flipping
#' @param scoring_matrix Use the same as in amplicanPipeline.
#' @param gapOpening (numeric)
#' @param gapExtension (numeric)
#' @param type (character) Same as in pairwiseAlignment.
#' @return (numeric) Read ids from reads that are HDR.
#'
hdr_read_ids <- function(reads, donor, hdr_events, scoring_matrix,
                         gapOpening = 25, gapExtension = 0, type = "overlap") {
  alignD <- Biostrings::pairwiseAlignment(
    DNAStringSet(gsub("-", "", pattern(reads))),
    DNAStringSet(toupper(donor)),
    substitutionMatrix = scoring_matrix,
    type = type, gapOpening = gapOpening, gapExtension = gapExtension)
  better_scores <- score(alignD) > score(reads)
  if (sum(better_scores) == 0) return(c())
  # comparison <- compareStrings(pattern(alignD[better_scores]),
  #                              subject(alignD[better_scores]))
  # comparison <- IRanges::RleList(strsplit(comparison, split = ""))
  # comparison <- IRanges::IRangesList(comparison %in% c("?", "+", "-"))
  # comparison <- IRanges::shift(
  #   comparison, start(subject(alignD[better_scores])) - 1)
  #
  # pure_hdr <- !as.logical(sum(IRanges::overlapsAny(
  #   comparison,
  #   IRanges::IRangesList(rep(list(hdr_events), length(comparison))),
  #   type = "any")))

  seq_len(length(reads))[better_scores]#[pure_hdr]
}

#' Find events that are correct HDR reads.
#'
#' Before using this function make sure events are filtered to represent
#' consensus with \code{amplicanConsensus}, if you use both forward and
#' reverse reads. If you want to calculate metrics over expected cut site,
#' filter events using \code{amplicanOverlap}.
#'
#' Takes donor sequences, aligns them against amplicons to determine expected
#' HDR events. Next searches for those events in aln data.frame and returns
#' T/F vector corresponding to every input event.
#'
#' @param aln (data.frame) Contains events from the alignments. Events have to
#' be filtered, shifted and normalized already.
#' @param aeSet (AlignmentExperimentSet) Contains whole alignments. Load with
#' readRDS().
#' @param cfgT (data.frame) Config file with the experiments details.
#' @param donors (character vector) Vector of donor templates, where sequences
#' correspond to the experiments from the cfgt$ID column.
#' @param scoring_matrix Use the same as in amplicanPipeline.
#' @param gapOpening (numeric)
#' @param gapExtension (numeric)
#' @param type (character) Same as in pairwiseAlignment.
#' @return (data.frame) As cfgT, but with extra columns.
#' @export
#' @family analysis steps
#' @include helpers_general.R
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "events_filtered_shifted_normalized.csv",
#'                          package = "amplican")
#' aln <- data.table::fread(file_path)
#' aeSet <- readRDS(system.file("extdata", "results", "alignments",
#'                              "AlignmentsExperimentSet.rds",
#'                              package = "amplican"))
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' !all(amplicanHDR(aln, aeSet, cfgT, cfgT$Amplicon))
#' # HDR will be zero as donors are exact amplicons and there is no events that
#' # comply
#'
#cfgT <- this_r
#donors <- donor$donor
amplicanHDR <- function(
  aln, aeSet, cfgT, donors,
  scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
    match = 5, mismatch = -4, baseOnly = TRUE, type = "DNA"),
  gapOpening = 25, gapExtension = 0, type = "overlap") {

  colnames(cfgT) <- c("ID", "Barcode", "Forward_Reads_File",
                      "Reverse_Reads_File", "Group", "guideRNA", "Found_Guide",
                      "Control", "Forward_Primer", "Reverse_Primer", "Direction",
                      "Amplicon", "Donor",
                      if (dim(cfgT)[2] > 13) colnames(cfgT)[14:dim(cfgT)[2]])
  seqnames <- read_id <- counts <- start <- end <- score <- NULL

  align <- Biostrings::pairwiseAlignment(
    DNAStringSet(toupper(donors)),
    DNAStringSet(toupper(cfgT$Amplicon)),
    substitutionMatrix = scoring_matrix,
    type = "overlap", gapOpening = 25, gapExtension = 0)
  pat <- pattern(align)
  subj <-  subject(align)
  scores <- score(align)

  hdr <- c()
  for (i in seq_len(nrow(cfgT))) {
    # extract events we want to find to quantify read as fully HDR
    this_id <- amplican::getEvents(pat[i], subj[i], scores = scores[i],
                                   ID = cfgT$ID[i], strand_info = "+",
                                   ampl_start = start(subj)[i])
    names(this_id) <- NULL
    hdr <- c(hdr, this_id)
  }
  hdr <- unlist(GenomicRanges::GRangesList(hdr), use.names = FALSE)
  names(hdr) <- NULL
  hdr <- as.data.frame(hdr)
  data.table::setDF(cfgT)
  hdr_nf <- flipRanges(hdr, cfgT)
  hdr_nf <- data.frame(amplicanMap(hdr_nf, cfgT), stringsAsFactors = FALSE)
  data.table::setDT(hdr_nf)

  # align to the donor
  HDR <- rep(FALSE, nrow(aln))

  for (i in seq_len(nrow(cfgT))) {

    whichI <- which(names(aeSet) == cfgT$ID[i])
    this_hdr <- hdr[hdr$seqnames == cfgT$ID[i], ]
    this_hdr <- IRanges::IRanges(this_hdr$start, this_hdr$end)

    hdr_ids_fwd <- hdr_read_ids(
      fwdReads(aeSet[whichI])[[1]], donors[i], this_hdr,
      scoring_matrix = scoring_matrix,
      gapOpening = gapOpening,
      gapExtension = gapExtension,
      type = type)
    HDR[aln$seqnames == cfgT$ID[i] &
          aln$read_id %in% hdr_ids_fwd & aln$strand == "-"] <- TRUE

    hdr_ids_rve <- hdr_read_ids(
      rveReads(aeSet[whichI])[[1]], donors[i], this_hdr,
      scoring_matrix = scoring_matrix,
      gapOpening = gapOpening,
      gapExtension = gapExtension,
      type = type)
    HDR[aln$seqnames == cfgT$ID[i] &
          aln$read_id %in% hdr_ids_rve & aln$strand == "-"] <- TRUE
  }
  HDR
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


#' Summarize how many reads comes from HDR.
#'
#' Before using this function make sure events are filtered to represent
#' consensus with \code{amplicanConsensus}, if you use both forward and
#' reverse reads. If you want to calculate metrics over expected cut site,
#' filter events using \code{amplicanOverlap}.
#'
#' Adds columns to cfgT:
#' \itemize{
#' \item{HDR}{ Count of reads with HDR events}
#' }
#' @param aln (data.frame) Contains events from the alignments.
#' @param cfgT (data.frame) Config file with the experiments details.
#' @param HDR (character) Name of the column from aln that contains HDR
#' information.
#' @return (data.frame) As cfgT, but with extra columns.
#' @export
#' @family analysis steps
#' @include helpers_general.R
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "events_filtered_shifted_normalized.csv",
#'                          package = "amplican")
#' aln <- data.table::fread(file_path)
#' aeSet <- readRDS(system.file("extdata", "results", "alignments",
#'                              "AlignmentsExperimentSet.rds",
#'                              package = "amplican"))
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' # check which events come from HDR, in the example amplicons are donor
#' # templates so we expect no events
#' aln$HDR <- amplicanHDR(aln, aeSet, cfgT, cfgT$Amplicon)
#' # finally calculate HDR counts per experiment
#' amplicanSummarizeHDR(aln, cfgT)
#'
amplicanSummarizeHDR <- function(aln, cfgT, HDR = "HDR") {

  seqnames <- read_id <- counts <- NULL
  cfgT$HDR <- 0
  aln_noD <- unique(aln[aln$`HDR`, ], by = c("seqnames", "read_id"))
  widthT_final <- aln_noD[, list(counts = sum(counts)), by = seqnames]
  map <- match(widthT_final$seqnames, cfgT$ID)
  cfgT$HDR[map] <- widthT_final$counts

  cfgT
}
