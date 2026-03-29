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
                            min_overlap = ceiling(nchar(primer)/2)) {
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
      match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA"),
    gapOpening = 1, # Disallow indels
    gapExtension = 1)

  valid_indices <- which(pwalign::nmatch(pwa) >= min_overlap & pwalign::nedit(pwa) <= m)
  results <- rep(NA_real_, length(reads))
  if (length(valid_indices) == 0) {
    return(results)
  }
  results[valid_indices] <- start(pwalign::pattern(pwa))[valid_indices]
  return(results)
}

is_hdr <- function(reads, scores, amplicon, donor, type = "overlap",
                   scoring_matrix, gap_opening = 25, gap_extension = 0,
                   donor_mismatch = 3) {

  align <- pwalign::pairwiseAlignment(
    DNAStringSet(toupper(donor)), DNAStringSet(toupper(amplicon)),
    substitutionMatrix = scoring_matrix, type = type,
    gapOpening = gap_opening, gapExtension = gap_extension)
  pat <- pattern(align)
  subj <-  subject(align)
  # extract events we want to find to quantify read as fully HDR
  hdr_events <- amplican::getEvents(pat, subj, scores = score(align),
                                 ID = "HDR", strand_info = "+",
                                 ampl_start = start(subj))
  names(hdr_events) <- NULL
  hdr_events <- IRanges::ranges(hdr_events)

  # Fast exact match bypass
  reads_dna <- DNAStringSet(reads)
  donor_dna <- DNAStringSet(toupper(donor))
  is_exact_donor <- reads_dna == donor_dna

  is_hdr <- rep(FALSE, length(reads))
  if (all(is_exact_donor)) {
    is_hdr[] <- TRUE
    return(is_hdr)
  }

  reads_to_align <- reads[!is_exact_donor]

  # now align reads to donor
  alignD <- pwalign::pairwiseAlignment(reads_to_align,
    donor_dna,
    type = type, substitutionMatrix = scoring_matrix,
    gapOpening = gap_opening, gapExtension = gap_extension)
  better_scores <- score(alignD) >= scores[!is_exact_donor]

  if (sum(better_scores) == 0) {
    is_hdr[is_exact_donor] <- TRUE
    return(is_hdr)
  }
  comparison <- pwalign::compareStrings(pattern(alignD[better_scores]),
                                        subject(alignD[better_scores]))
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
  comparison <- IRanges::shift(comparison,  shft - 1)
  overlaps_hdr <- IRanges::overlapsAny(comparison, hdr_events, type = "any")
  all_e_not_overlap <- sapply(split(!overlaps_hdr, names(comparison)), all)
  if (length(all_e_not_overlap) > 0) {
    all_e_not_overlap <- as.integer(names(all_e_not_overlap)[all_e_not_overlap])
  } else {
    all_e_not_overlap <- NULL
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
    as.integer(names(overlaps_e[overlaps_e <= 0])))

  is_hdr[is_exact_donor] <- TRUE
  is_hdr[!is_exact_donor][better_scores][ok_hdr] <- TRUE
  is_hdr
}


#' Figure out which reads conform to the HDR using the donor.
#'
#' This is strict detection as compared to `is_hdr` which was designed to be
#' less specific and allow for all kinds of donors. This method requires that
#' you have exactly the same events (mismatches, insertions, deletions) as the difference
#' between amplicon and donor sequences. It ignores everything else, so other mismatches and small
#' indels etc. as noise are allowed here for valid HDR.
#'
#' @param aln (data.table) This are events that contain already consensus column,
#' they are also shifted and normalized.
#' @param cfgT (data.table) Config data.table with columns for amplicon and donor.
#' @param scoring_matrix (scoring matrix)
#' @param gap_opening (integer)
#' @param gap_extension (integer)
#' @export
#' @return (aln) same as aln on entry, but readType is updated to TRUE when read is recognized as HDR
#'
is_hdr_strict <- function(aln, cfgT, scoring_matrix,
                          gap_opening = 25,
                          gap_extension = 0) {
  . <- NULL

  for (i in seq_len(dim(cfgT)[1])) {
    donor <- get_seq(cfgT, cfgT$ID[i], "Donor")
    if (donor == "") next()
    amplicon <- get_seq(cfgT, cfgT$ID[i])

    # donor vs amplicon
    d_a_aln <- pwalign::pairwiseAlignment(
      DNAStringSet(toupper(donor)),
      DNAStringSet(toupper(amplicon)),
      substitutionMatrix = scoring_matrix, type = "overlap",
      gapOpening = gap_opening, gapExtension = gap_extension)
    pat <- pattern(d_a_aln)
    subj <-  subject(d_a_aln)
    # extract events we want to find to quantify read as fully HDR
    hdr_events <- amplican::getEvents(pat, subj, scores = score(d_a_aln),
                                      ID = cfgT$ID[i], strand_info = "+",
                                      ampl_start = start(subj))
    if (length(hdr_events) == 0) next()
    hdr_events <- amplicanMap(hdr_events, cfgT)

    # this is strict algorithm, using binary lookups instead of sequential vectors
    events <- aln[.(cfgT$ID[i]), on = "seqnames", nomatch = NULL]
    if (nrow(events) == 0) next()

    events <- events[consensus == TRUE]
    if (nrow(events) == 0) next()

    hits <- data.table::merge.data.table(as.data.table(events),
                                         as.data.table(hdr_events),
                                         all.x = FALSE, all.y = FALSE,
                                         by = c("start", "end", "width",
                                                "originally", "replacement", "type"))
    if (nrow(hits) == 0) next()
    hits <- as.data.table(hits)
    hits <- hits[, .(n = .N), by = "read_id.x"]
    valid_reads <- hits$read_id.x[hits$n == length(hdr_events)]
    if(length(valid_reads) == 0) next()

    aln[.(cfgT$ID[i]), readType := read_id %in% valid_reads, on = "seqnames"]
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
                          temp_folder = NULL) {

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
    fwdStream <- ShortRead::FastqStreamer(cfgT$Forward_Reads_File[1], n = batch_size)
    on.exit(close(fwdStream), add = TRUE)
  }
  if (fastqfiles != 1) {
    rveStream <- ShortRead::FastqStreamer(cfgT$Reverse_Reads_File[1], n = batch_size)
    on.exit(close(rveStream), add = TRUE)
  }

  unqT_list <- list()
  bad_base_quality <- 0
  bad_average_quality <- 0
  bad_alphabet <- 0
  read_count <- 0
  filtered_read_count <- 0

  repeat {
    fwdT <- if (fastqfiles != 2) ShortRead::yield(fwdStream) else NULL
    rveT <- if (fastqfiles != 1) ShortRead::yield(rveStream) else NULL

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

      chunk_unqT <- data.frame(
        Forward = if (fastqfiles == 2) "" else as.character(ShortRead::sread(fwdT_good)),
        Reverse = if (fastqfiles == 1) "" else as.character(ShortRead::sread(rveT_good)),
        stringsAsFactors = FALSE
      )
      chunk_unqT$Total <- paste0(chunk_unqT$Forward, chunk_unqT$Reverse)
      chunk_unqT <- stats::aggregate(Total ~ Forward + Reverse, chunk_unqT, length)
      unqT_list[[length(unqT_list) + 1]] <- chunk_unqT
    }
  }

  barcodeTable <- data.frame(Barcode = barcode,
                             experiment_count = length(unique(cfgT$ID)),
                             read_count = read_count,
                             bad_base_quality = bad_base_quality,
                             bad_average_quality = bad_average_quality,
                             bad_alphabet = bad_alphabet,
                             filtered_read_count = filtered_read_count,
                             stringsAsFactors = FALSE)

  if (length(unqT_list) > 0) {
    unqT <- data.table::rbindlist(unqT_list)
    unqT <- unqT[, .(Total = sum(Total)), by = .(Forward, Reverse)]
    data.table::setDF(unqT)
  } else {
    unqT <- data.frame(Forward=character(), Reverse=character(), Total=integer(), stringsAsFactors=FALSE)
  }

  if (dim(unqT)[1] == 0) {
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
                        barcodeData = barcodeTable)
    if (!is.null(temp_folder)) {
      temp_file <- file.path(temp_folder, paste0(barcode, "_aln.rds"))
      temp_file_writing <- file.path(temp_folder, paste0(barcode, "_aln.rds.temp"))
      saveRDS(aes, temp_file_writing)
      file.rename(temp_file_writing, temp_file)
      return(temp_file)
    }
    return(aes)
  }

  unqT$BarcodeFrequency <- unqT$Total / sum(unqT$Total)
  unqT <- unqT[order(unqT$Forward, unqT$Reverse), ]
  unqT$Asigned <- FALSE
  unqT[c("Forward", "Reverse")] <- lapply(unqT[c("Forward", "Reverse")],
                                          function(x) toupper(as.character(x)))
  barcodeTable$unique_reads <- nrow(unqT)

  fwd_primer_cache <- list()
  rve_primer_cache <- list()

  # for each experiment
  for (i in seq_len(dim(cfgT)[1])) {

    # Primers and amplicon info
    fwdPrimer <- toupper(cfgT$Forward_Primer[i])
    rvePrimer <- toupper(cfgT$Reverse_Primer[i])
    amplicon <- toupper(cfgT$Amplicon[i])
    donor <- toupper(cfgT$Donor[i])

    # Search for the forward, reverse and targets
    if (fastqfiles == 2 | fwdPrimer == "") {
      unqT$fwdPrInReadPos <- NA
      unqT$forwardFound <- FALSE
    } else {
      if (!fwdPrimer %in% names(fwd_primer_cache)) {
        fwd_primer_cache[[fwdPrimer]] <- locate_pr_start(unqT$Forward, fwdPrimer, primer_mismatch)
      }
      unqT$fwdPrInReadPos <- fwd_primer_cache[[fwdPrimer]]
      unqT$forwardFound <- is.finite(unqT$fwdPrInReadPos)
    }

    if (fastqfiles == 1 | rvePrimer == "") {
      unqT$rvePrInReadPos <- NA
      unqT$reverseFound <- FALSE
    } else {
      if (!rvePrimer %in% names(rve_primer_cache)) {
        rve_primer_cache[[rvePrimer]] <- locate_pr_start(unqT$Reverse, rvePrimer, primer_mismatch)
      }
      unqT$rvePrInReadPos <- rve_primer_cache[[rvePrimer]]
      unqT$reverseFound <- is.finite(unqT$rvePrInReadPos)
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
    unqT$Asigned <- unqT$Asigned | primersFound
    IDunqT <- unqT[primersFound, ]
    # most abundant fwd + rve combination from the top
    IDunqT <- IDunqT[order(-IDunqT$Total), ]

    if (dim(IDunqT)[1] > 0) {
      # when some reads match to correct primers
      # do alignments with removed dna bases before primers to allow
      # reads and amplicons start with primer
      # when calling events into GRanges shift_ampl is used to adapt for
      # subtractions happening here
      if (fastqfiles != 2) {
        rF <- Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Forward"]),
                                 start = IDunqT$fwdPrInReadPos)
        fwdA[[cfgT$ID[i]]] <-
          pwalign::pairwiseAlignment(
            rF,
            Biostrings::subseq(amplicon,
                               start = cfgT$fwdPrPos[i],
                               end = cfgT$rvePrPosEnd[i]),
            type = "overlap", substitutionMatrix = scoring_matrix,
            gapOpening = gap_opening, gapExtension = gap_extension)

        if (donor != "") {
          fwdAType[[cfgT$ID[i]]] <- if (donor_strict) {
            rep(FALSE, length(rF))
          } else {
            is_hdr(
              rF, score(fwdA[[cfgT$ID[i]]]),
              amplicon, donor,
              type = "overlap", scoring_matrix =  scoring_matrix,
              gap_opening = gap_opening, gap_extension = gap_extension,
              donor_mismatch = donor_mismatch)
          }
        }
      }

      if (fastqfiles != 1) {
        rR <- Biostrings::reverseComplement(
          Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Reverse"]),
                             start = IDunqT$rvePrInReadPos))
        rveA[[cfgT$ID[i]]] <- pwalign::pairwiseAlignment(
          rR,
          Biostrings::subseq(amplicon,
                             start = cfgT$fwdPrPos[i],
                             end = cfgT$rvePrPosEnd[i]),
          type = "overlap", substitutionMatrix = scoring_matrix,
          gapOpening = gap_opening, gapExtension = gap_extension)

        if (donor != "") {
          rveAType[[cfgT$ID[i]]] <- if (donor_strict) {
            rep(FALSE, length(rR))
          } else {
            is_hdr(
              rR, score(rveA[[cfgT$ID[i]]]),
              amplicon, donor,
              type = "overlap", scoring_matrix =  scoring_matrix,
              gap_opening = gap_opening, gap_extension = gap_extension,
              donor_mismatch = donor_mismatch)
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

  if (dim(unassignedTable)[1] > 0) {
    unassignedTable$Barcode <- barcode
    unassignedTable <- unassignedTable[order(-unassignedTable$Total),]
  } else {
    unassignedTable <- NULL
  }

  aes <- methods::new("AlignmentsExperimentSet",
               fwdReads = fwdA,
               rveReads = rveA,
               fwdReadsType = fwdAType,
               rveReadsType = rveAType,
               readCounts = countsA,
               unassignedData = unassignedTable,
               experimentData = cfgT,
               barcodeData = barcodeTable)
  if (!is.null(temp_folder)) {
    temp_file <- file.path(temp_folder, paste0(barcode, "_aln.rds"))
    temp_file_writing <- file.path(temp_folder, paste0(barcode, "_aln.rds.temp"))
    saveRDS(aes, temp_file_writing)
    file.rename(temp_file_writing, temp_file)
    return(temp_file)
  }
  return(aes)
}
