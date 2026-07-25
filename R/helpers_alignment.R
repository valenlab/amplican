#' Generate all combinations along string exchanging m characters at a time with
#' dictionary letters.
#'
#' Generate all combinations along string \code{seq} swapping \code{m}
#' characters at a time with letters defined in dictionary \code{letters}.
#' Allows, for instance, to create a list of possible primers with two
#' mismatches.
#'
#' @param seq (character) input character to permutate
#' @param m (integer) number of elements to permutate at each step
#' @param letters (character vector) dictionary source for combinations of
#' elements
#' @return (character vector) all unique combinations of permutated string
#' @export
#' @examples
#' comb_along("AC")
#' comb_along("AAA", 1)
#' comb_along("AAA")
#' comb_along("AAA", 3)
#' comb_along("AAAAAAAAAA")
#'
comb_along <- function(seq, m = 2, letters = c("A", "C", "T", "G")) {
  seq <- as.list(strsplit(seq, "")[[1]])
  indices <- utils::combn(seq_along(seq), m)
  letters <- list(letters)

  seq <- apply(indices, 2, function(x) {
    seq[x] <- letters
    do.call(paste0, expand.grid(seq))
  })

  unique(as.vector(seq))
}

#' Find full or partial primer start positions with mismatches.
#'
#' This function finds a primer that can be partially truncated at the 5' end.
#' It handles both full matches inside the read and partial matches at the
#' beginning of the read.
#'
#' @param reads A character vector of DNA sequences.
#' @param primer A single character string for the primer sequence.
#' @param m The maximum number of allowed mismatches.
#' @param min_overlap The minimum number of base pairs the primer must overlap
#'   with the read. A good value is often around half the primer length.
#' @return A numeric vector of the EARLIEST start position for a valid match,
#'   with NA for no match.
locate_pr_start <- function(reads, primer,
                            m = 3,
                            min_overlap = ceiling(nchar(primer) / 2)) {
  primer_len <- nchar(primer)
  if (min_overlap > primer_len) {
    stop("min_overlap cannot be greater than the primer length.")
  }

  # 1. Align the single full primer against all reads using overlap alignment.
  # This is the most robust way to handle 5' truncations and mismatches.
  pwa <- pwalign::pairwiseAlignment(
    pattern = Biostrings::DNAStringSet(reads),
    subject = Biostrings::DNAString(primer),
    type = "overlap",
    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
      match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA"
    ),
    gapOpening = 1, # Disallow indels
    gapExtension = 1
  )

  valid_indices <- which(pwalign::nmatch(pwa) >= min_overlap & pwalign::nedit(pwa) <= m)
  results <- rep(NA_real_, length(reads))
  if (length(valid_indices) == 0) {
    return(results)
  }
  results[valid_indices] <- start(pwalign::pattern(pwa))[valid_indices]
  return(results)
}

#' Determine which reads conform to HDR using a donor template (permissive).
#'
#' Aligns each read to the donor and accepts it as HDR when two conditions hold:
#' \enumerate{
#'   \item The read aligns to the donor at least as well as to the amplicon
#'         (\code{score(read vs donor) >= score(read vs amplicon)}).
#'   \item The number of events in the read that \emph{overlap the donor-event
#'         positions} (i.e. the positions that differ between donor and amplicon)
#'         and have width 1 does not exceed \code{donor_mismatch} after subtracting
#'         the donor events themselves.
#' }
#'
#' \strong{Important}: \code{donor_mismatch} counts only events whose coordinates
#' overlap the donor-vs-amplicon event window. Mismatches or indels in the read
#' that fall \emph{outside} that window are invisible to this threshold and never
#' disqualify a read.
#'
#' \strong{Score tie behaviour}: reads that align equally well to the donor and
#' the amplicon (equal scores) are treated as candidates (\code{>=}) and proceed
#' to the event-overlap check.
#'
#' Use \code{\link{is_hdr_strict}} when you require every donor-specific event to
#' be present verbatim in the read.
#'
#' @param reads (\code{\link[Biostrings]{DNAStringSet}}) Aligned reads.
#' @param scores (numeric) Alignment scores of \code{reads} against the amplicon
#'   (e.g. from \code{pwalign::pairwiseAlignment}).
#' @param amplicon (character) Amplicon sequence (single string).
#' @param donor (character) Donor template sequence (single string).
#' @param type (character) Alignment type passed to
#'   \code{\link[pwalign]{pairwiseAlignment}}.  Default \code{"overlap"}.
#' @param scoring_matrix Substitution matrix (e.g. from
#'   \code{\link[pwalign]{nucleotideSubstitutionMatrix}}).
#' @param gap_opening (numeric) Gap-opening penalty. Default 25.
#' @param gap_extension (numeric) Gap-extension penalty. Default 0.
#' @param donor_mismatch (numeric) Maximum number of width-1 events (single-base
#'   mismatches, single-base deletions, single-base insertions) that are allowed
#'   to overlap the donor-event positions.  Only events within the donor-vs-amplicon
#'   event coordinate window are counted; events elsewhere in the read are ignored.
#'   Set to 0 to require the donor region to match perfectly (note: sequencing error
#'   rate makes 0 inadvisable in practice). Default 3.
#' @keywords internal
#' @return (logical vector) TRUE for each read classified as HDR.
#' @seealso \code{\link{is_hdr_strict}} for an event-presence-based alternative.
is_hdr <- function(reads, scores, amplicon, donor, type = "overlap",
                   scoring_matrix, gap_opening = 25, gap_extension = 0,
                   donor_mismatch = 3) {
  align <- pwalign::pairwiseAlignment(
    DNAStringSet(toupper(donor)), DNAStringSet(toupper(amplicon)),
    substitutionMatrix = scoring_matrix, type = type,
    gapOpening = gap_opening, gapExtension = gap_extension
  )
  pat <- pattern(align)
  subj <- subject(align)
  # extract events we want to find to quantify read as fully HDR
  hdr_events <- amplican::getEvents(pat, subj,
    scores = score(align),
    ID = "HDR", strand_info = "+",
    ampl_start = start(subj)
  )
  names(hdr_events) <- NULL
  hdr_events <- IRanges::ranges(hdr_events)

  # now align reads to donor
  alignD <- pwalign::pairwiseAlignment(reads,
    DNAStringSet(toupper(donor)),
    type = type, substitutionMatrix = scoring_matrix,
    gapOpening = gap_opening, gapExtension = gap_extension
  )
  better_scores <- score(alignD) >= scores
  is_hdr <- rep(FALSE, length(reads))
  if (sum(better_scores) == 0) {
    return(is_hdr)
  }
  comparison <- pwalign::compareStrings(
    pattern(alignD[better_scores]),
    subject(alignD[better_scores])
  )
  comparison <- IRanges::RleList(strsplit(comparison, split = ""))

  mm <- IRanges::IRangesList(comparison == "?") # need to tile
  reads_ids <- seq_along(comparison)
  names(mm) <- reads_ids
  mm <- unlist(mm, use.names = TRUE)
  mm_l <- mm[IRanges::width(mm) > 1]
  mm <- mm[IRanges::width(mm) == 1]
  mm_ln <- names(mm_l)
  mm_l <- IRanges::tile(mm_l, width = 1L)
  names(mm_l) <- mm_ln
  mm_l <- unlist(mm_l, use.names = TRUE)
  mm <- c(mm, mm_l)

  comparison <- IRanges::IRangesList(comparison %in% c("+", "-"))
  names(comparison) <- seq_along(comparison)
  comparison <- unlist(comparison, use.names = TRUE)
  comparison <- c(comparison, mm)

  shft <- start(subject(alignD[better_scores]))[as.numeric(names(comparison))]
  comparison <- IRanges::shift(comparison, shft - 1)
  overlaps_hdr <- IRanges::overlapsAny(comparison, hdr_events, type = "any")
  all_e_not_overlap <- sapply(split(!overlaps_hdr, names(comparison)), all)
  if (length(all_e_not_overlap) > 0) {
    all_e_not_overlap <- as.integer(names(all_e_not_overlap)[all_e_not_overlap])
  } else {
    all_e_not_overlap <- integer(0)
  }

  # tolerate some noise level
  # no events + events not overlapping hdr + allow n events of length 1 in those
  # overlapping
  overlaps_e <- comparison[overlaps_hdr]
  overlaps_e <- IRanges::IRangesList(split(overlaps_e, names(overlaps_e)))
  overlaps_e_w1 <- sum(IRanges::width(overlaps_e) == 1)
  overlaps_e_w1[overlaps_e_w1 > donor_mismatch] <- donor_mismatch
  overlaps_e <- sum(IRanges::width(overlaps_e))
  overlaps_e <- overlaps_e - overlaps_e_w1

  ok_hdr <- c(
    reads_ids[!reads_ids %in% as.integer(names(comparison))],
    all_e_not_overlap,
    as.integer(names(overlaps_e[overlaps_e <= 0]))
  )

  is_hdr[better_scores][ok_hdr] <- TRUE
  is_hdr
}


#' Determine which reads conform to HDR using the donor (strict, event-presence).
#'
#' This is the strict counterpart to \code{\link{is_hdr}}.  A read is marked HDR
#' if and only if \emph{every} event that distinguishes the donor from the amplicon
#' is present verbatim (same \code{start}, \code{end}, \code{width},
#' \code{originally}, \code{replacement}, and \code{type}) in that read's
#' consensus events.
#'
#' \strong{Key behaviours}:
#' \itemize{
#'   \item Only rows where \code{consensus == TRUE} are used to match donor events.
#'         However, the \code{readType} flag is written back to \emph{all} rows
#'         sharing the same \code{read_id} (including non-consensus rows).
#'   \item When \code{donor_mismatch = Inf} (default), additional events in the
#'         read (noise mismatches, extra indels) do \strong{not} disqualify a
#'         read — only the \emph{absence} of a required donor event does.
#'   \item When \code{donor_mismatch} is finite (e.g. 0), extra consensus events
#'         \emph{within the amplicon UPPERCASE window} (expanded by
#'         \code{cut_buffer}) are counted.  If more than \code{donor_mismatch}
#'         extra events fall in that window, the read is rejected.  Events
#'         outside the window (e.g. near primers) are ignored.
#'   \item A different substitution at the same position as a donor event (e.g.
#'         the donor has G->A but the read has G->C) is rejected because the
#'         \code{replacement} column differs in the inner-join merge.
#'   \item If the donor and amplicon are identical (no events), or if no read
#'         events match any donor event, the function returns \code{aln} unchanged.
#' }
#'
#' @param aln (data.table) Consensus-filtered, shifted, and normalised events
#'   table (must contain columns \code{seqnames}, \code{read_id}, \code{consensus},
#'   \code{readType}, \code{start}, \code{end}, \code{width}, \code{originally},
#'   \code{replacement}, \code{type}).
#' @param cfgT (data.table or data.frame) Config table with at least columns
#'   \code{ID}, \code{Amplicon}, \code{Donor}, and \code{Direction}.
#' @param scoring_matrix Substitution matrix passed to
#'   \code{\link[pwalign]{pairwiseAlignment}} for the donor-vs-amplicon alignment.
#' @param gap_opening (numeric) Gap-opening penalty. Default 25.
#' @param gap_extension (numeric) Gap-extension penalty. Default 0.
#' @param donor_mismatch (numeric) Maximum number of extra consensus events
#'   (beyond the required donor events) allowed within the amplicon UPPERCASE
#'   window.  Set to \code{Inf} (default) to allow unlimited noise, or to
#'   \code{0} to require no extra events in the window.
#' @param cut_buffer (numeric) Number of bases to expand the UPPERCASE window
#'   on each side.  Same semantics as in \code{\link{amplicanOverlap}}.
#'   Default 5.
#' @export
#' @return (data.table) Same as \code{aln} on entry, but \code{readType} is set
#'   to \code{TRUE} for every row whose \code{read_id} contains all donor events
#'   and at most \code{donor_mismatch} additional events in the window.
#' @seealso \code{\link{is_hdr}} for the permissive, score-based alternative.
#'
is_hdr_strict <- function(aln, cfgT, scoring_matrix,
                          gap_opening = 25,
                          gap_extension = 0,
                          donor_mismatch = Inf,
                          cut_buffer = 5) {
  setDT(aln)

  for (i in seq_len(dim(cfgT)[1])) {
    amplicon <- get_seq(cfgT, cfgT$ID[i])
    donor <- get_seq(cfgT, cfgT$ID[i], "Donor")
    aln_id <- !is.na(aln$seqnames) & aln$seqnames == cfgT$ID[i]

    if (!any(aln_id) | donor == "") next()

    # donor vs amplicon
    d_a_aln <- pwalign::pairwiseAlignment(
      DNAStringSet(toupper(donor)),
      DNAStringSet(toupper(amplicon)),
      substitutionMatrix = scoring_matrix, type = "overlap",
      gapOpening = gap_opening, gapExtension = gap_extension
    )
    pat <- pattern(d_a_aln)
    subj <- subject(d_a_aln)
    # extract events we want to find to quantify read as fully HDR
    hdr_events <- amplican::getEvents(pat, subj,
      scores = score(d_a_aln),
      ID = cfgT$ID[i], strand_info = "+",
      ampl_start = start(subj)
    )
    if (length(hdr_events) == 0) next()
    hdr_events <- amplicanMap(hdr_events, cfgT)

    # this is strict algorithm
    # we take only consensus events
    events <- aln[aln_id & aln$consensus, ]
    if (nrow(events) == 0) next()

    hits <- data.table::merge.data.table(as.data.table(events),
      as.data.table(hdr_events),
      all.x = F, all.y = F,
      by = c(
        "start", "end", "width",
        "originally", "replacement", "type"
      )
    )
    if (nrow(hits) == 0) next()
    hits <- as.data.table(hits)
    hits <- hits[, .(n = .N), by = "read_id.x"]
    hits <- hits$read_id.x[hits$n == length(hdr_events)] # make sure all events are represented

    # enforce donor_mismatch: count extra events within the UPPERCASE window
    if (is.finite(donor_mismatch) && length(hits) > 0) {
      ug <- upperGroups(amplicon)
      if (length(ug) > 0) {
        # shift to relative coords (matching amplicanMap) and expand
        ug <- IRanges::shift(ug, -1L * IRanges::start(ug)[1])
        ug <- ug + cut_buffer
        # count events in window per hit read
        hit_events <- events[read_id %in% hits]
        hit_ranges <- IRanges::IRanges(
          start = hit_events$start, end = hit_events$end
        )
        in_window <- IRanges::overlapsAny(hit_ranges, ug)
        windowed <- hit_events[in_window, .(n = .N), by = "read_id"]
        # subtract donor events that fall in the window
        hdr_dt <- as.data.table(hdr_events)
        hdr_ranges <- IRanges::IRanges(
          start = hdr_dt$start, end = hdr_dt$end
        )
        n_hdr_in_window <- sum(IRanges::overlapsAny(hdr_ranges, ug))
        windowed[, extra := n - n_hdr_in_window]
        bad <- windowed$read_id[windowed$extra > donor_mismatch]
        hits <- hits[!hits %in% bad]
      }
    }
    aln[aln_id, readType := read_id %in% hits]
  }
  return(aln)
}


#' Make alignments helper.
#'
#' Aligning reads to the amplicons for each ID in this barcode, constructing
#' amplicanAlignment. Assume that all IDs here belong to the same barcode.
#' @keywords internal
#' @param cfgT config file as data table
#' @inheritParams amplicanAlign
#' @include helpers_general.R helpers_filters.R AlignmentsExperimentSet-class.R
#' @return amplicanAlignment object for this barcode experiments
#'
makeAlignment <- function(cfgT,
                          average_quality,
                          min_quality,
                          filter_n,
                          batch_size,
                          scoring_matrix,
                          gap_opening,
                          gap_extension,
                          fastqfiles,
                          primer_mismatch,
                          donor_mismatch,
                          donor_strict,
                          temp_folder = NULL,
                          sample = 0,
                          seed = 0) {
  barcode <- cfgT$Barcode[1]

  if (!is.null(temp_folder)) {
    temp_file <- file.path(temp_folder, paste0(barcode, "_aln.rds"))
    if (file.exists(temp_file)) {
      message("Skipping alignments for ", barcode, " (already exists)")
      return(temp_file)
    }
  }

  message("Aligning reads for ", barcode)

  fwdA <- vector("list", length(cfgT$ID))
  names(fwdA) <- cfgT$ID
  rveA <- countsA <- fwdAType <- rveAType <- fwdA # pre-allocate alignment lists

  # Read Reads for this Barcode
  if (fastqfiles != 2) {
    if (sample > 0) {
      fwdStream <- ShortRead::FastqSampler(cfgT$Forward_Reads_File[1], n = sample)
    } else {
      fwdStream <- ShortRead::FastqStreamer(cfgT$Forward_Reads_File[1], n = batch_size)
    }
    on.exit(close(fwdStream), add = TRUE)
  }
  if (fastqfiles != 1) {
    if (sample > 0) {
      rveStream <- ShortRead::FastqSampler(cfgT$Reverse_Reads_File[1], n = sample)
    } else {
      rveStream <- ShortRead::FastqStreamer(cfgT$Reverse_Reads_File[1], n = batch_size)
    }
    on.exit(close(rveStream), add = TRUE)
  }

  unqT_list <- list()
  chunk_count <- 0L
  bad_base_quality <- 0
  bad_average_quality <- 0
  bad_alphabet <- 0
  read_count <- 0
  filtered_read_count <- 0

  sampled_already <- FALSE

  repeat {
    if (sample > 0 && sampled_already) break
    if (sample > 0) set.seed(seed)
    fwdT <- if (fastqfiles != 2) ShortRead::yield(fwdStream) else NULL
    if (sample > 0) set.seed(seed)
    rveT <- if (fastqfiles != 1) ShortRead::yield(rveStream) else NULL
    sampled_already <- TRUE

    if (fastqfiles == 1) {
      if (length(fwdT) == 0) break
    } else if (fastqfiles == 2) {
      if (length(rveT) == 0) break
    } else {
      if (length(fwdT) == 0 && length(rveT) == 0) break
    }

    if (fastqfiles == 1) {
      rveT <- rep(TRUE, length(fwdT))
    }
    if (fastqfiles == 2) {
      fwdT <- rep(TRUE, length(rveT))
    }

    read_count <- read_count + length(fwdT)

    # Filter Reads
    goodq <- goodBaseQuality(fwdT, min = min_quality, batch_size = batch_size) &
      goodBaseQuality(rveT, min = min_quality, batch_size = batch_size)
    avrq <- goodAvgQuality(fwdT, avg = average_quality, batch_size = batch_size) &
      goodAvgQuality(rveT, avg = average_quality, batch_size = batch_size)
    nucq <- if (filter_n) {
      alphabetQuality(fwdT, batch_size = batch_size) &
        alphabetQuality(rveT, batch_size = batch_size)
    } else {
      rep(TRUE, length(avrq))
    }
    goodReads <- goodq & avrq & nucq

    bad_base_quality <- bad_base_quality + sum(!goodq)
    bad_average_quality <- bad_average_quality + sum(!avrq)
    bad_alphabet <- bad_alphabet + sum(!nucq)
    filtered_read_count <- filtered_read_count + sum(goodReads)

    if (sum(goodReads) > 0) {
      fwdT_good <- fwdT[goodReads]
      rveT_good <- rveT[goodReads]

      chunk_unqT <- data.table::data.table(
        Forward = if (fastqfiles == 2) "" else as.character(ShortRead::sread(fwdT_good)),
        Reverse = if (fastqfiles == 1) "" else as.character(ShortRead::sread(rveT_good))
      )
      chunk_unqT <- chunk_unqT[, .(Total = .N), by = .(Forward, Reverse)]
      chunk_count <- chunk_count + 1L
      unqT_list[[chunk_count]] <- chunk_unqT
    }
  }

  barcodeTable <- data.frame(
    Barcode = barcode,
    experiment_count = length(unique(cfgT$ID)),
    read_count = read_count,
    bad_base_quality = bad_base_quality,
    bad_average_quality = bad_average_quality,
    bad_alphabet = bad_alphabet,
    filtered_read_count = filtered_read_count,
    stringsAsFactors = FALSE
  )

  if (chunk_count > 0) {
    unqT <- data.table::rbindlist(unqT_list[seq_len(chunk_count)])
    unqT <- unqT[, .(Total = sum(Total)), by = .(Forward, Reverse)]
  } else {
    unqT <- data.table::data.table(Forward = character(), Reverse = character(), Total = integer())
  }

  if (nrow(unqT) == 0) {
    barcodeTable$unique_reads <- 0
    barcodeTable$unassigned_reads <- 0
    barcodeTable$assigned_reads <- 0
    aes <- methods::new("AlignmentsExperimentSet",
      fwdReads = fwdA,
      rveReads = rveA,
      fwdReadsType = fwdAType,
      rveReadsType = rveAType,
      readCounts = countsA,
      unassignedData = NULL,
      experimentData = cfgT,
      barcodeData = barcodeTable
    )
    if (!is.null(temp_folder)) {
      temp_file <- file.path(temp_folder, paste0(barcode, "_aln.rds"))
      temp_file_writing <- file.path(temp_folder, paste0(barcode, "_aln.rds.temp"))
      saveRDS(aes, temp_file_writing)
      file.rename(temp_file_writing, temp_file)
      return(temp_file)
    }
    return(aes)
  }

  data.table::set(unqT, j = "BarcodeFrequency", value = unqT$Total / sum(unqT$Total))
  data.table::setorder(unqT, Forward, Reverse)
  data.table::set(unqT, j = "Asigned", value = FALSE)
  data.table::set(unqT, j = "Forward", value = toupper(as.character(unqT$Forward)))
  data.table::set(unqT, j = "Reverse", value = toupper(as.character(unqT$Reverse)))
  barcodeTable$unique_reads <- nrow(unqT)

  fwd_primer_cache <- list()
  rve_primer_cache <- list()

  # for each experiment
  n_cfg <- nrow(cfgT)
  for (i in seq_len(n_cfg)) {
    # Primers and amplicon info
    fwdPrimer <- toupper(cfgT$Forward_Primer[i])
    rvePrimer <- toupper(cfgT$Reverse_Primer[i])
    amplicon <- toupper(cfgT$Amplicon[i])
    donor <- toupper(cfgT$Donor[i])

    # Search for the forward, reverse and targets
    if (fastqfiles == 2 | fwdPrimer == "") {
      data.table::set(unqT, j = "fwdPrInReadPos", value = NA)
      data.table::set(unqT, j = "forwardFound", value = FALSE)
    } else {
      if (!fwdPrimer %in% names(fwd_primer_cache)) {
        fwd_primer_cache[[fwdPrimer]] <- locate_pr_start(unqT$Forward, fwdPrimer, primer_mismatch)
      }
      data.table::set(unqT, j = "fwdPrInReadPos", value = fwd_primer_cache[[fwdPrimer]])
      data.table::set(unqT, j = "forwardFound", value = is.finite(unqT$fwdPrInReadPos))
    }

    if (fastqfiles == 1 | rvePrimer == "") {
      data.table::set(unqT, j = "rvePrInReadPos", value = NA)
      data.table::set(unqT, j = "reverseFound", value = FALSE)
    } else {
      if (!rvePrimer %in% names(rve_primer_cache)) {
        rve_primer_cache[[rvePrimer]] <- locate_pr_start(unqT$Reverse, rvePrimer, primer_mismatch)
      }
      data.table::set(unqT, j = "rvePrInReadPos", value = rve_primer_cache[[rvePrimer]])
      data.table::set(unqT, j = "reverseFound", value = is.finite(unqT$rvePrInReadPos))
    }

    primersFound <-
      if (fastqfiles == 0.5) {
        unqT$forwardFound | unqT$reverseFound
      } else if (fastqfiles == 1) {
        unqT$forwardFound
      } else if (fastqfiles == 2) {
        unqT$reverseFound
      } else {
        unqT$forwardFound & unqT$reverseFound
      }
    data.table::set(unqT, j = "Asigned", value = unqT$Asigned | primersFound)
    IDunqT <- unqT[primersFound, ]
    # most abundant fwd + rve combination from the top
    data.table::setorder(IDunqT, -Total)

    if (nrow(IDunqT) > 0) {
      # when some reads match to correct primers
      # do alignments with removed dna bases before primers to allow
      # reads and amplicons start with primer
      # when calling events into GRanges shift_ampl is used to adapt for
      # subtractions happening here
      if (fastqfiles != 2) {
        rF <- Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[["Forward"]]),
          start = IDunqT$fwdPrInReadPos
        )
        fwdA[[cfgT$ID[i]]] <-
          pwalign::pairwiseAlignment(
            rF,
            Biostrings::subseq(amplicon,
              start = cfgT$fwdPrPos[i],
              end = cfgT$rvePrPosEnd[i]
            ),
            type = "overlap", substitutionMatrix = scoring_matrix,
            gapOpening = gap_opening, gapExtension = gap_extension
          )

        if (donor != "") {
          fwdAType[[cfgT$ID[i]]] <- if (donor_strict) {
            rep(FALSE, length(rF))
          } else {
            is_hdr(
              rF, score(fwdA[[cfgT$ID[i]]]),
              amplicon, donor,
              type = "overlap", scoring_matrix = scoring_matrix,
              gap_opening = gap_opening, gap_extension = gap_extension,
              donor_mismatch = donor_mismatch
            )
          }
        }
      }

      if (fastqfiles != 1) {
        rR <- Biostrings::reverseComplement(
          Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[["Reverse"]]),
            start = IDunqT$rvePrInReadPos
          )
        )
        rveA[[cfgT$ID[i]]] <- pwalign::pairwiseAlignment(
          rR,
          Biostrings::subseq(amplicon,
            start = cfgT$fwdPrPos[i],
            end = cfgT$rvePrPosEnd[i]
          ),
          type = "overlap", substitutionMatrix = scoring_matrix,
          gapOpening = gap_opening, gapExtension = gap_extension
        )

        if (donor != "") {
          rveAType[[cfgT$ID[i]]] <- if (donor_strict) {
            rep(FALSE, length(rR))
          } else {
            is_hdr(
              rR, score(rveA[[cfgT$ID[i]]]),
              amplicon, donor,
              type = "overlap", scoring_matrix = scoring_matrix,
              gap_opening = gap_opening, gap_extension = gap_extension,
              donor_mismatch = donor_mismatch
            )
          }
        }
      }
      countsA[[cfgT$ID[i]]] <- IDunqT$Total
    }
    cfgT$Reads[i] <- sum(IDunqT$Total)
  }

  barcodeTable$unassigned_reads <- sum(!unqT$Asigned)
  barcodeTable$assigned_reads <- sum(unqT$Asigned)
  unassignedTable <- unqT[!unqT$Asigned, ]

  if (nrow(unassignedTable) > 0) {
    data.table::set(unassignedTable, j = "Barcode", value = barcode)
    data.table::setorder(unassignedTable, -Total)
    data.table::setDF(unassignedTable)
  } else {
  }

  aes <- methods::new("AlignmentsExperimentSet",
    fwdReads = fwdA,
    rveReads = rveA,
    fwdReadsType = fwdAType,
    rveReadsType = rveAType,
    readCounts = countsA,
    unassignedData = unassignedTable,
    experimentData = cfgT,
    barcodeData = barcodeTable
  )
  if (!is.null(temp_folder)) {
    temp_file <- file.path(temp_folder, paste0(barcode, "_aln.rds"))
    temp_file_writing <- file.path(temp_folder, paste0(barcode, "_aln.rds.temp"))
    saveRDS(aes, temp_file_writing)
    file.rename(temp_file_writing, temp_file)
    return(temp_file)
  }
  return(aes)
}
