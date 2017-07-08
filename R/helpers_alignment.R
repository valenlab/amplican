#' Make alignments helper.
#'
#' Aligning reads to the amplicons for each ID in this barcode, constructing
#' amplicanAlignment. Assume that all IDs here belong to the same barcode.
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
                          fastqfiles) {

  barcode <- cfgT$Barcode[1]
  message(paste0("Aligning reads for ", barcode))

  fwdA <- vector("list", length(cfgT$ID))
  names(fwdA) <- cfgT$ID
  rveA <- fwdA # pre-allocate alignment lists
  countsA <- fwdA

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

  barcodeTable <- data.frame(barcode = barcode,
                             experiment_count = length(unique(cfgT$ID)),
                             read_count = length(goodReads),
                             bad_base_quality = sum(!goodq),
                             bad_average_quality = sum(!avrq),
                             bad_alphabet = sum(!nucq),
                             filtered_read_count = sum(goodReads))

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

    # Search for the forward, reverse and targets
    unqT$fwdPrInReadPos <- stringr::str_locate(unqT$Forward, fwdPrimer)[,1]
    unqT$forwardFound <-
      if (fastqfiles == 2 | fwdPrimer == "") {
        FALSE
      } else {
        !is.na(unqT$fwdPrInReadPos)
      }
    unqT$rvePrInReadPos <- stringr::str_locate(unqT$Reverse, rvePrimer)[,1]
    unqT$reverseFound <-
      if (fastqfiles == 1) {
        FALSE
      } else {
        !is.na(unqT$rvePrInReadPos)
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
      # substractions happening here
      if (fastqfiles != 2) {
        fwdA[[cfgT$ID[i]]] <-
          Biostrings::pairwiseAlignment(
            Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Forward"]),
                               start = IDunqT$fwdPrInReadPos),
            Biostrings::subseq(amplicon,
                               start = cfgT$fwdPrPos[i]),
            type = "global",
            substitutionMatrix =  scoring_matrix,
            gapOpening = gap_opening,
            gapExtension = gap_extension)
      }

      if (fastqfiles != 1) {
        rveA[[cfgT$ID[i]]] <- Biostrings::pairwiseAlignment(
          Biostrings::reverseComplement(
            Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Reverse"]),
                               start = IDunqT$rvePrInReadPos)),
          Biostrings::subseq(amplicon,
                             start = 1,
                             end = cfgT$rvePrPosEnd[i]),
          type = "global",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension)
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
    unassignedTable <- data.frame()
  }

  methods::new("AlignmentsExperimentSet",
               fwdReads = fwdA,
               rveReads = rveA,
               readCounts = countsA,
               unassignedData = unassignedTable,
               experimentData = cfgT,
               barcodeData = barcodeTable)
}
