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


locate_pr_start <- function(reads, primer, m = 0) {
  primer <- comb_along(primer, m)
  primer <- sapply(primer, function(pr) {
    stringr::str_locate(reads, pr)[, 1]
  }, simplify = TRUE, USE.NAMES = FALSE)
  reads <- matrixStats::rowMaxs(primer, na.rm = TRUE)
  reads[!is.finite(reads)] <- NA
  reads
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
                          scoring_matrix,
                          gap_opening,
                          gap_extension,
                          fastqfiles,
                          primer_mismatch) {

  barcode <- cfgT$Barcode[1]
  message("Aligning reads for ", barcode)

  fwdA <- vector("list", length(cfgT$ID))
  names(fwdA) <- cfgT$ID
  rveA <- countsA <- fwdAType <- rveAType <- fwdA # pre-allocate alignment lists

  # Read Reads for this Barcode
  fwdT <- if (fastqfiles == 2) NULL else ShortRead::readFastq(
    cfgT$Forward_Reads_File[1])
  rveT <- if (fastqfiles == 1) NULL else ShortRead::readFastq(
    cfgT$Reverse_Reads_File[1])
  if (fastqfiles == 1) {
    rveT <- rep(TRUE, length(fwdT))
  }
  if (fastqfiles == 2) {
    fwdT <- rep(TRUE, length(rveT))
  }

  # Filter Reads
  goodq <- goodBaseQuality(fwdT, min = min_quality) &
    goodBaseQuality(rveT, min = min_quality)
  avrq <- goodAvgQuality(fwdT, avg = average_quality) &
    goodAvgQuality(rveT, avg = average_quality)
  nucq <- alphabetQuality(fwdT) & alphabetQuality(rveT)
  goodReads <- goodq & avrq & nucq

  barcodeTable <- data.frame(Barcode = barcode,
                             experiment_count = length(unique(cfgT$ID)),
                             read_count = length(goodReads),
                             bad_base_quality = sum(!goodq),
                             bad_average_quality = sum(!avrq),
                             bad_alphabet = sum(!nucq),
                             filtered_read_count = sum(goodReads),
                             stringsAsFactors = FALSE)

  fwdT <- fwdT[goodReads]
  rveT <- rveT[goodReads]

  # Unique reads
  unqT <- data.frame(
    if (fastqfiles == 2) "" else as.character(ShortRead::sread(fwdT)),
    if (fastqfiles == 1) "" else as.character(ShortRead::sread(rveT)))
  colnames(unqT) <- c("Forward", "Reverse")
  unqT$Total <- paste0(unqT$Forward, unqT$Reverse)
  unqT <- stats::aggregate(Total ~ Forward + Reverse, unqT, length)
  unqT$BarcodeFrequency <- unqT$Total / sum(unqT$Total)
  unqT <- unqT[order(unqT$Forward, unqT$Reverse), ]
  unqT$Asigned <- FALSE
  unqT[c("Forward", "Reverse")] <- lapply(unqT[c("Forward", "Reverse")],
                                          function(x) toupper(as.character(x)))
  barcodeTable$unique_reads <- nrow(unqT)

  # for each experiment
  for (i in seq_len(dim(cfgT)[1])) {

    # Primers and amplicon info
    fwdPrimer <- toupper(cfgT$Forward_Primer[i])
    rvePrimer <- toupper(cfgT$Reverse_Primer[i])
    amplicon <- toupper(cfgT$Amplicon[i])
    donor <- toupper(cfgT$Donor[i])

    # Search for the forward, reverse and targets
    unqT$fwdPrInReadPos <- locate_pr_start(
      unqT$Forward, fwdPrimer, primer_mismatch)
    unqT$forwardFound <-
      if (fastqfiles == 2 | fwdPrimer == "") {
        FALSE
      } else {
        is.finite(unqT$fwdPrInReadPos)
      }
    unqT$rvePrInReadPos <- locate_pr_start(
      unqT$Reverse, rvePrimer, primer_mismatch)
    unqT$reverseFound <-
      if (fastqfiles == 1) {
        FALSE
      } else {
        is.finite(unqT$rvePrInReadPos)
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
        fwdA[[cfgT$ID[i]]] <-
          Biostrings::pairwiseAlignment(
            Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Forward"]),
                               start = IDunqT$fwdPrInReadPos),
            Biostrings::subseq(amplicon,
                               start = cfgT$fwdPrPos[i],
                               end = cfgT$rvePrPosEnd[i]),
            type = "overlap",
            substitutionMatrix =  scoring_matrix,
            gapOpening = gap_opening,
            gapExtension = gap_extension)
        if (donor != "") {
          donorA <-
            Biostrings::pairwiseAlignment(
              Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Forward"]),
                                 start = IDunqT$fwdPrInReadPos),
              Biostrings::DNAString(donor),
              type = "overlap",
              substitutionMatrix =  scoring_matrix,
              gapOpening = gap_opening,
              gapExtension = gap_extension)
          fwdAType[[cfgT$ID[i]]] <- score(donorA) > score(fwdA[[cfgT$ID[i]]])
        }
      }

      if (fastqfiles != 1) {
        rveA[[cfgT$ID[i]]] <- Biostrings::pairwiseAlignment(
          Biostrings::reverseComplement(
            Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Reverse"]),
                               start = IDunqT$rvePrInReadPos)),
          Biostrings::subseq(amplicon,
                             start = cfgT$fwdPrPos[i],
                             end = cfgT$rvePrPosEnd[i]),
          type = "overlap",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension)

        if (donor != "") {
          donorA <- Biostrings::pairwiseAlignment(
            Biostrings::reverseComplement(
              Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Reverse"]),
                                 start = IDunqT$rvePrInReadPos)),
            Biostrings::DNAString(donor),
            type = "overlap",
            substitutionMatrix =  scoring_matrix,
            gapOpening = gap_opening,
            gapExtension = gap_extension)
          rveAType[[cfgT$ID[i]]] <- score(donorA) > score(rveA[[cfgT$ID[i]]])
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

  methods::new("AlignmentsExperimentSet",
               fwdReads = fwdA,
               rveReads = rveA,
               fwdReadsType = fwdAType,
               rveReadsType = rveAType,
               readCounts = countsA,
               unassignedData = unassignedTable,
               experimentData = cfgT,
               barcodeData = barcodeTable)
}
