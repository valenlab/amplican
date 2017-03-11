#' Helper to construct GRanges with additional metadata columns.
#'
#' @param x (IRanges)
#' @param ID (string)
#' @param type (string)
#' @param strand_info (string) Either '+', '-'
#' @param originally (string) Base pairs on the amplicon.
#' @param replacement (string) Base pairs on the read.
#' @return (GRanges) Object with meta-data
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#'
defGR <- function(x,
                  ID,
                  strand_info = "+",
                  type = "deletion",
                  originally = "",
                  replacement = "") {

  x <- GRanges(
    ranges = x,
    strand = Rle(rep(strand_info, length(x))),
    seqnames = Rle(rep(ID, length(x)))
  )
  x$originally = originally
  x$replacement = replacement
  x$type = type
  return(x)
}

#' This function takes alignments and gives back the events coordinates.
#'
#' @param align (PairwiseAlignmentsSingleSubject)
#' @param ID (string)
#' @param ampl_shift (numeric) Shift events additionaly by this value.
#' PairwiseAlignmentsSingleSubject returns truncated alignments.
#' @param strand_info (string) Either '+', '-' or default '*'
#' @return (List of GRanges) Object with meta-data for insertion, deletion, mismatch
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#'
getEventInfo <- function(align, ID, ampl_shift, strand_info = "+") {

  del <- deletion(align)[[1]]
  del <- shift(del, c(0, cumsum(width(del))[-length(del)]))
  ins <- insertion(align)[[1]]
  mm <- mismatchSummary(summary(align))$subject

  if (length(ins) > 0 & dim(mm)[1] > 0) {
    ss <- sapply(mm$SubjectPosition, function(x) {
      sum(width(ins)[x > start(ins)])
    })
    mm$SubjectPosition <- mm$SubjectPosition + ss
  }
  ins <- shift(ins, c(0, cumsum(width(ins))[-length(ins)]))

  finalGR <- GRanges()

  if (length(del) > 0) {
    finalGR <- defGR(del, ID, strand_info)
  }

  if (length(ins) > 0) {
    ins_seq <-
      substr(rep(as.character(pattern(align)), length(ins)),
                 start = start(ins), stop = end(ins))
    finalGR <- c(finalGR, defGR(ins, ID, strand_info, "insertion", "", ins_seq))
  }

  if (dim(mm)[1] > 0) {
    finalGR <- c(finalGR,
                   defGR(IRanges(mm$SubjectPosition, width = 1),
                         ID,
                         strand_info,
                         "mismatch",
                         mm$Subject,
                         mm$Pattern))
  }

  if (strand_info == "+") {
    finalGR <- shift(finalGR, ampl_shift - 1)
  } else {
    subj_bas <- stringr::str_count(as.character(subject(align)), "[ATCG]")
    finalGR <- shift(finalGR, ampl_shift - subj_bas)
  }

  return(finalGR)
}


#' Detect uppercases as ranges object.
#'
#' For a given string, detect how many groups of uppercases is inside, where
#' are they, and how long they are.
#'
#' For example:
#'    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
#'
#' Has 4 groups of uppercases of length 7, 4, 1 and 3.
#' @param candidate (string) A string with the nucleotide sequence.
#' @import IRanges
#' @return (Ranges) A Ranges object with uppercases groups for given candidate
#' string
#' @importFrom stringr str_detect
#'
upperGroups <- function(candidate) {
  return(IRanges::reduce(IRanges(
    start = which(stringr::str_detect(strsplit(candidate, "")[[1]], "[[:upper:]]")),
    width = 1
  )))
}


#' Make alignments helper.
#'
#' Main functionality of the package, aligning reads to the amplicon.
#' @param configTable config file as data table
#' @param resultsFolder (string) path to resultsFolder
#' @inheritParams amplicanAlign
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @importFrom ShortRead readFastq sread
#' @include helpers_general.R helpers_filters.R
#' @return resultsFolder as invisible
#'
makeAlignment <- function(configTable,
                          resultsFolder,
                          skip_bad_nucleotides,
                          average_quality,
                          min_quality,
                          write_alignments,
                          scoring_matrix,
                          gap_opening,
                          gap_extension,
                          fastqfiles,
                          PRIMER_DIMER,
                          cut_buffer) {

  div <- c("deletion", "insertion")
  barcode <- configTable$Barcode[1]
  alignmentRangesBar <- GRanges()
  GenomeInfoDb::seqlevels(alignmentRangesBar) <- unique(configTable$ID)

  message(paste0("Aligning reads for ", barcode))

  # Read Reads for this Barcode
  forwardsTable <- if (fastqfiles == 2) NULL else ShortRead::readFastq(configTable$Forward_Reads_File[1])
  reversesTable <- if (fastqfiles == 1) NULL else ShortRead::readFastq(configTable$Reverse_Reads_File[1])
  if (fastqfiles == 1) {
    rewersesTable <- rep(TRUE, length(forwardsTable))
  }
  if (fastqfiles == 2) {
    forwardsTable <- rep(TRUE, length(reversesTable))
  }

  # Filter Reads
  goodq <- goodBaseQuality(forwardsTable, min = min_quality) & goodBaseQuality(reversesTable, min = min_quality)
  avrq <- goodAvgQuality(forwardsTable, avg = average_quality) & goodAvgQuality(reversesTable, avg = average_quality)
  nucq <- alphabetQuality(forwardsTable) & alphabetQuality(reversesTable)
  goodReads <- goodq & avrq & nucq

  barcodeTable <- data.frame(barcode = barcode,
                             experiment_count = configTable$ExperimentsCount[1],
                             read_count = length(goodReads),
                             bad_base_quality = sum(!goodq),
                             bad_average_quality = sum(!avrq),
                             bad_alphabet = sum(!nucq),
                             filtered_read_count = sum(goodReads))

  forwardsTable <- forwardsTable[goodReads]
  reversesTable <- reversesTable[goodReads]

  # Unique reads
  uniqueTable <- data.frame(if (fastqfiles == 2) "" else as.character(sread(forwardsTable)),
                            if (fastqfiles == 1) "" else as.character(sread(reversesTable)))
  colnames(uniqueTable) <- c("Forward", "Reverse")
  uniqueTable$Total <- paste0(uniqueTable$Forward, uniqueTable$Reverse)
  uniqueTable <- stats::aggregate(Total ~ Reverse + Forward, uniqueTable, length)
  uniqueTable$BarcodeFrequency <- uniqueTable$Total / sum(uniqueTable$Total)
  uniqueTable <- uniqueTable[order(uniqueTable$Forward, uniqueTable$Reverse), ]
  uniqueTable$Asigned <- FALSE
  uniqueTable$PRIMER_DIMER <- FALSE
  uniqueTable[c("Reverse", "Forward")] <- lapply(uniqueTable[c("Reverse", "Forward")],
                                                 function(x) toupper(as.character(x)))
  barcodeTable$unique_reads <- nrow(uniqueTable)

  # for each experiment
  for (i in seq_len(dim(configTable)[1])) {

    alignmentRanges <- GRanges()
    currentID <- configTable[i, "ID"]
    # Primers and amplicon info
    forwardPrimer <- toupper(configTable[i, "Forward_Primer"])
    reversePrimer <- toupper(configTable[i, "Reverse_Primer"])
    guideRNA <- toupper(configTable[i, "guideRNA"])
    rguideRNA <- toupper(configTable[i, "RguideRNA"])
    amplicon <- configTable[i, "Amplicon"]

    # Search for the forward, reverse and targets
    uniqueTable$forPrInReadPos <- stringr::str_locate(uniqueTable$Forward, forwardPrimer)[,1]
    uniqueTable$forwardFound <-
      if (fastqfiles == 2 | forwardPrimer == "") {
        FALSE
      } else {
        !is.na(uniqueTable$forPrInReadPos)
      }
    uniqueTable$revPrInReadPos <- stringr::str_locate(uniqueTable$Reverse, reversePrimer)[,1]
    uniqueTable$reverseFound <-
      if (fastqfiles == 1) {
        FALSE
      } else {
        !is.na(uniqueTable$revPrInReadPos)
      }
    uniqueTable$guideFoundForward <- stringr::str_detect(uniqueTable$Forward, guideRNA)
    uniqueTable$guideFoundReverse <- stringr::str_detect(uniqueTable$Reverse, rguideRNA)
    primersFound <-
      if (fastqfiles == 0.5) {
        uniqueTable$forwardFound | uniqueTable$reverseFound
      } else if (fastqfiles == 1) {
        uniqueTable$forwardFound
      } else if (fastqfiles == 2) {
        uniqueTable$reverseFound
      } else {
        uniqueTable$forwardFound & uniqueTable$reverseFound
      }
    uniqueTable$Asigned <- uniqueTable$Asigned | primersFound
    IDuniqueTable <- uniqueTable[primersFound, ]

    if (dim(IDuniqueTable)[1] > 0) {

      if (fastqfiles != 2) {
        alignForward <- Biostrings::pairwiseAlignment(
          Biostrings::subseq(Biostrings::DNAStringSet(IDuniqueTable[, "Forward"]),
                             start = IDuniqueTable$forPrInReadPos),
          Biostrings::subseq(toupper(amplicon),
                             start = configTable$forwardPrimerPosition[i]),
          type = "global",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension
        )
      }
      if (fastqfiles != 1) {
        alignReverse <- Biostrings::pairwiseAlignment(
          Biostrings::reverseComplement(Biostrings::subseq(Biostrings::DNAStringSet(IDuniqueTable[, "Reverse"])),
                                        start = IDuniqueTable$revPrInReadPos),
          Biostrings::subseq(toupper(amplicon),
                             start = 1,
                             end = configTable$reversePrimerPosEnd[i]),
          type = "global",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension
        )
      }
      # Write the alignments
      if (write_alignments >= 1) {
        algn_file <-
          file.path(resultsFolder, paste0(currentID, "_", barcode, ".txt"))
        if (file.exists(algn_file)) {
          file.remove(algn_file)
        }
        algn_file_con <- file(algn_file, open = "at")
        writeLines(as.vector(rbind(
          paste("ID:", currentID,
                "Count:", format(IDuniqueTable$Total)),
          if (fastqfiles != 2) rbind(as.character(pattern(alignForward)),
                                     as.character(subject(alignForward))),
          if (fastqfiles < 1) "",
          if (fastqfiles != 1) rbind(
            as.character(pattern(alignReverse)),
            as.character(subject(alignReverse))),
          ""
        )),  algn_file_con)
        close(algn_file_con)
      }

      for (r in seq_len(dim(IDuniqueTable)[1])) {

        forwardData <- getEventInfo(alignForward[r], currentID, configTable$forwardPrimerPosition[i], "+")
        reverseData <- getEventInfo(alignReverse[r], currentID, configTable$reversePrimerPosEnd[i], "-")

        # Filter PRIMER DIMERS and sum how many
        PD_cutoff <- nchar(amplicon) -
          (nchar(forwardPrimer) + nchar(reversePrimer) + PRIMER_DIMER)
        isPD <- any(c(width(forwardData),
                      width(reverseData)) > PD_cutoff)
        IDuniqueTable[r, "PRIMER_DIMER"] <- isPD
        configTable$PRIMER_DIMER[i] <- configTable$PRIMER_DIMER[i] +
          isPD * IDuniqueTable$Total[r]
        if (isPD) {
          next
        }

        # Filters
        forwardDataFiltered <-
          forwardData[!(forwardData$type %in% div &
                          end(forwardData) >= configTable$reversePrimerPosition[i])]
        reverseDataFiltered <-
          reverseData[!(reverseData$type %in% div &
                          end(reverseData) >= configTable$reversePrimerPosition[i])]
        forwardDataFiltered <-
          forwardData[!(forwardData$type %in% div &
                          start(forwardData) <= configTable$forwardPrimerPositionEnd[i])]
        reverseDataFiltered <-
          reverseData[!(reverseData$type %in% div &
                          start(reverseData) <= configTable$forwardPrimerPositionEnd[i])]

        # Frameshift table
        frameshift <- FALSE
        if (fastqfiles == 0 | fastqfiles == 0.5) {
          frameshift <-
            sum(width(forwardDataFiltered[forwardDataFiltered$type == "insertion" |
                                            forwardDataFiltered$type == "deletion"])) %%
            3 !=  0 &
            sum(width(reverseDataFiltered[reverseDataFiltered$type == "insertion" |
                                            reverseDataFiltered$type == "deletion"])) %%
            3 != 0
        } else if (fastqfiles == 1) {
          frameshift <-
            sum(width(forwardDataFiltered[forwardDataFiltered$type == "insertion" |
                                            forwardDataFiltered$type == "deletion"])) %%
            3 != 0
        } else {
          frameshift <-
            sum(width(reverseData[reverseDataFiltered$type == "insertion" |
                                    reverseDataFiltered$type == "deletion"])) %%
            3 != 0
        }
        configTable$Frameshift[i] <- configTable$Frameshift[i] +
          frameshift * IDuniqueTable$Total[r]

        # definitions
        if (length(forwardData) > 0) {
          forwardData$cut <- FALSE
          forwardData$count <- IDuniqueTable$Total[r]
          forwardData$frequency <- 0  #prepare field
          forwardData$read_id <- r
        }

        if (length(reverseData) > 0) {
          reverseData$cut <- FALSE
          reverseData$count <- IDuniqueTable$Total[r]
          reverseData$frequency <- 0  #prepare field
          reverseData$read_id <- r
        }

        # cut assessment
        overlapFd <- subsetByOverlaps(ranges(forwardDataFiltered[forwardDataFiltered$type == "deletion"]),
                                      configTable$cutSites[[i]])
        overlapRe <- subsetByOverlaps(ranges(reverseDataFiltered[reverseDataFiltered$type == "deletion"]),
                                      configTable$cutSites[[i]])
        if (fastqfiles == 0 | fastqfiles == 0.5) {
          # forward and reverse have to agree on deletion
          overlapFd <- overlapFd[!is.na(match(overlapFd, overlapRe))]
          overlapRe <- overlapRe[!is.na(match(overlapRe, overlapFd))]
          forwardData$cut[!is.na(match(ranges(forwardData), overlapFd))] <- TRUE
          reverseData$cut[!is.na(match(ranges(reverseData), overlapRe))] <- TRUE
        } else if (fastqfiles == 1) {
          forwardData$cut[!is.na(match(ranges(forwardData), overlapFd))] <- TRUE
        } else {
          reverseData$cut[!is.na(match(ranges(reverseData), overlapRe))] <- TRUE
        }

        configTable$Cut[i] <- configTable$Cut[i] +
          if (fastqfiles == 2) {
            any(reverseData$cut) * IDuniqueTable$Total[r]
          } else {
            any(forwardData$cut) * IDuniqueTable$Total[r]
          }

        if (length(forwardData) > 0 | length(reverseData) > 0) {
          alignmentRanges <- c(alignmentRanges,
                               forwardData,
                               reverseData)
        }
      }

      configTable$Reads[i] <- sum(IDuniqueTable$Total)
      # fill frequency
      if (length(alignmentRanges) > 0) {
        reads_count_filtered <- configTable$Reads[i] - configTable$PRIMER_DIMER[i]
        alignmentRanges$frequency <- alignmentRanges$count / reads_count_filtered
      }
    }

    alignmentRangesBar <- c(alignmentRangesBar, alignmentRanges)
  }

  if (length(alignmentRangesBar) > 0) {
    utils::write.csv(as.data.frame(alignmentRangesBar),
                     file.path(resultsFolder,
                               paste0(barcode, "_alignment_ranges.csv")),
                     row.names = FALSE)
  }

  barcodeTable$unassigned_reads <- sum(!uniqueTable$Asigned)
  barcodeTable$assigned_reads <- sum(uniqueTable$Asigned)
  utils::write.csv(uniqueTable[!uniqueTable$Asigned, ],
                   file = file.path(resultsFolder,
                                    "unassigned_sequences",
                                    paste0(barcode, "_unassigned_reads.csv")),
                   quote = FALSE, row.names = FALSE)

  utils::write.csv(configTable[c("ID",
                                 "Barcode",
                                 "Forward_Reads_File",
                                 "Reverse_Reads_File",
                                 "Group","guideRNA",
                                 "Forward_Primer",
                                 "Reverse_Primer",
                                 "Direction",
                                 "Amplicon",
                                 "ExperimentsCount",
                                 "Cut",
                                 "Frameshift",
                                 "PRIMER_DIMER",
                                 "Reads",
                                 "Found_Guide",
                                 "Found_PAM")],
                   file.path(resultsFolder,
                             paste0(barcode, "_configFile_results.csv")),
                   row.names = FALSE)
  utils::write.csv(barcodeTable,
                   file.path(resultsFolder,
                             paste0(barcode, "_reads_filters.csv")),
                   row.names = FALSE)

  invisible(resultsFolder)
}
