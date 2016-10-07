#' This function takes a lite string and gives back the events coordinates.
#'
#' For the given example:
#' 0         1         2         3         4         5
#' 012345678901234567890123456789012345678901234567890123456789
#' NNNNNNNNNNN-----NNNNNNN-------NNNNNNNN--------NNNNNNNN------    Subject
#' ------NNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNN---NNNNNNNNNN    Pattern
#'
#' The RCPP library returns the information regarding the start and ends of the
#' events in the following way:
#' (Notice that there is a +1 shift in the index)
#'
#' #returned: @4@12,16*24,30*39,46*55,60*!@3@1,6*31,33*48,50*!@0@
#' #example: @0@!@1@152,207*!@1@25,t,A*
#' @param liteString (string)
#' @param ID (string)
#' @param strand (string) Either '+', '-' or default '*'
#' @param read_alignment (string) Sequence of aligned read to amplicon.
#' @return (GRanges) Object with meta-data for insertion, deletion, mismatch
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#'
getEventInfo <- function(liteString, ID, strand = "*", read_alignment) {

  if (is.na(liteString) | nchar(liteString) < 12) {
    return(GRanges())
  }

  alignmentData <- unlist(strsplit(liteString, "!", fixed = TRUE))
  insertionsData <- alignmentData[1]
  deletionsData <- alignmentData[2]
  mismatchesData <- alignmentData[3]

  insertionsR <- GRanges()
  insCount <- as.numeric(unlist(strsplit(insertionsData,
                                         "@",
                                         fixed = TRUE))[2])
  if (insCount > 0) {
    insertions <- unlist(strsplit(insertionsData,
                                  "@",
                                  fixed = TRUE))[3]
    insertions <- unlist(strsplit(insertions, "*", fixed = TRUE))
    insertionsPairs <- as.numeric(unlist(strsplit(insertions, ",",
                                                  fixed = TRUE)))
    insertionsR <- GRanges(
      ranges = IRanges(start = insertionsPairs[c(TRUE, FALSE)],
                       end = insertionsPairs[c(FALSE, TRUE)]),
      strand = S4Vectors::Rle(rep(strand, insCount)),
      seqnames = Rle(rep(ID, insCount))
    )
    insertionsR$mm_originally = ""
    insertionsR$mm_replacement = ""
    insertionsR$type = "insertion"
  }

  deletionsR <- GRanges()
  delCount <- as.numeric(unlist(strsplit(deletionsData, "@",
                                         fixed = TRUE))[2])
  if (delCount > 0) {
    deletions <- unlist(strsplit(deletionsData, "@", fixed = TRUE))[3]
    deletions <- unlist(strsplit(deletions, "*", fixed = TRUE))
    deletionsPairs <- as.numeric(unlist(strsplit(deletions, ",",
                                                 fixed = TRUE)))
    deletionsR <- GRanges(
      ranges = IRanges(start = deletionsPairs[c(TRUE, FALSE)],
                       end = deletionsPairs[c(FALSE, TRUE)]),
      strand = Rle(rep(strand, delCount)),
      seqnames = Rle(rep(ID, delCount))
    )
    deletionsR$mm_originally = ""
    deletionsR$mm_replacement = ""
    deletionsR$type = "deletion"
  }

  mismatchesR <- GRanges()
  mmCount <- as.numeric(unlist(strsplit(mismatchesData, "@",
                                        fixed = TRUE))[2])
  if (mmCount > 0) {
    mismatches <- unlist(strsplit(mismatchesData, "@", fixed = TRUE))[3]
    mismatches <- unlist(strsplit(mismatches, "*", fixed = TRUE))
    mismatches <- unlist(strsplit(mismatches, ",", fixed = TRUE))
    mismatchesR <- GRanges(
      ranges = IRanges(start = as.numeric(mismatches[c(TRUE, FALSE, FALSE)]),
                       width = 1),
      strand = Rle(rep(strand, mmCount)),
      seqnames = Rle(rep(ID, mmCount))
    )
    mismatchesR$mm_originally = mismatches[c(FALSE, TRUE, FALSE)]
    mismatchesR$mm_replacement = mismatches[c(FALSE, FALSE, TRUE)]
    mismatchesR$type = "mismatch"
  }

  finalEvents <- c(insertionsR, deletionsR, mismatchesR)
  return(
    if (length(finalEvents) > 0) {
      finalEvents
    } else {
      GRanges(ranges = IRanges(start = 0, end = 0),
              seqnames = ID,
              strand = strand,
              mm_originally = "",
              mm_replacement = "",
              type = "no_events")
    }
  )
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
#'
upperGroups <- function(candidate) {
  return(IRanges::reduce(
    IRanges(start = c(which(grepl("[[:upper:]]",
                                  strsplit(candidate, "")[[1]]))), width = 1)))
}


#' Make alignments helper.
#'
#' Main functionality of the package, aligning reads to the amplicon using
#' gotoh implementation.
#' @param configTable config file as data table
#' @param resultsFolder (string) path to resultsFolder
#' @inheritParams amplicanAlign
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @importFrom ShortRead readFastq sread
#' @return resultsFolder as invisible
#'
makeAlignment <- function(configTable,
                          resultsFolder,
                          skip_bad_nucleotides = TRUE,
                          average_quality = 0,
                          min_quality = 0,
                          write_alignments = 1,
                          scoring_matrix = "NUC44",
                          gap_opening = 50,
                          gap_extension = 0,
                          gap_ending = FALSE,
                          far_indels = TRUE,
                          fastqfiles = 0,
                          PRIMER_DIMER = 30,
                          cut_buffer = 5) {

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
                                                 as.character)

  barcodeTable$unique_reads <- nrow(uniqueTable)

  revDir <- configTable[, "Direction"] == 1
  configTable[revDir, "guideRNA"] <-
    as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAStringSet(
          configTable[revDir, "guideRNA"])))
  configTable$RReverse_Primer <-
    as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAStringSet(
          configTable$Reverse_Primer)))
  configTable$RguideRNA <-
    as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAStringSet(
          configTable$guideRNA)))

  for (i in seq_len(dim(configTable)[1])) {

    if (write_alignments >= 1) algn_file_con <- c()
    if (write_alignments >= 2) algn_det_file_con <- c()
    sublog_file_con <- ""
    # for each experiment
    alignmentRanges <- GRanges()
    currentID <- configTable[i, "ID"]
    # Primers and amplicon info
    forwardPrimer <- configTable[i, "Forward_Primer"]
    reversePrimer <- configTable[i, "Reverse_Primer"]
    guideRNA <- configTable[i, "guideRNA"]
    amplicon <- configTable[i, "Amplicon"]

    # Folder for each ID
    currentIDFolderName <- file.path(resultsFolder, paste0(currentID, "_", barcode))
    if (!dir.exists(currentIDFolderName)) {
      dir.create(file.path(currentIDFolderName), showWarnings = FALSE)
    }

    warning_checkTarget <- checkTarget(guideRNA, amplicon, currentID, barcode)
    configTable$Found_Guide[i] <- is.null(warning_checkTarget)
    sublog_file_con <- c(sublog_file_con, warning_checkTarget)
    sublog_file_con <- c(sublog_file_con,
                         checkPrimers(forwardPrimer,
                                      configTable[i, "RReverse_Primer"],
                                      amplicon, currentID, barcode))

    cutSites <- upperGroups(amplicon) + cut_buffer
    if (length(cutSites) == 0) {
      message("Warning: Aligment position was not found in the amplicon.
              Find more information in the log file.")
      sublog_file_con <- c(sublog_file_con,
                           paste0("Couldn't find upper case groups (PAM) position in amplicon.",
                                  "\nFor ID: ",
                                  currentID,
                                  " and barcode: ",
                                  barcode,
                                  "\n"))
      configTable$Found_PAM[i] <- 0
      cutSites <- IRanges(start = 1, width = nchar(amplicon))
    }

    # Search for the forward, reverse and targets
    uniqueTable$forwardFound <-
      if (fastqfiles == 2) {
        FALSE
      } else {
        grepl(forwardPrimer, uniqueTable$Forward, ignore.case = TRUE)
      }
    uniqueTable$reverseFound <-
      if (fastqfiles == 1) {
        FALSE
      } else {
        grepl(reversePrimer, uniqueTable$Reverse, ignore.case = TRUE)
      }
    uniqueTable$guideFoundForward <- grepl(guideRNA,
                                           uniqueTable$Forward,
                                           ignore.case = TRUE)
    uniqueTable$guideFoundReverse <- grepl(configTable[i, "RguideRNA"],
                                           uniqueTable$Reverse,
                                           ignore.case = TRUE)
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

    IDuniqueTable[c("Forward", "Reverse")] <-
      lapply(IDuniqueTable[c("Forward", "Reverse")], toupper)
    IDuniqueTable$Reverse <-
      as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(IDuniqueTable$Reverse)))

    if (dim(IDuniqueTable)[1] > 0) {
      for (r in seq_len(dim(IDuniqueTable)[1])) {
        # The return of gotoh is a string of five parts with the
        # separator '++++'
        # [1] is the verbose alignment result
        # [2] is a string with the ins and del and mm
        # [3] is a string with the ins and del and mm, but from
        #  the subject (amplicon) coordinates
        # [4] is the alignment of the pattern.
        # [5] is the alignment of the subject.
        alignForward <- c("", "", "", "", "")
        if (fastqfiles != 2) {
          alignForward <- gRCPP(IDuniqueTable[r, "Forward"],
                                amplicon,
                                scoring_matrix,
                                gap_opening,
                                gap_extension,
                                gap_ending,
                                far_indels)
          alignForward <- unlist(strsplit(alignForward, "++++", fixed = TRUE))
        }
        alignReverse <- c("", "", "", "", "")
        if (fastqfiles != 1) {
          alignReverse <- gRCPP(IDuniqueTable[r, "Reverse"],
                                amplicon,
                                scoring_matrix,
                                gap_opening,
                                gap_extension,
                                gap_ending,
                                far_indels)
          alignReverse <- unlist(strsplit(alignReverse, "++++", fixed = TRUE))
        }
        # Write the alignments
        if (write_alignments >= 2) {
          algn_det_file_con <- c(algn_det_file_con,
                                 paste("ID:", currentID,
                                       "Count:",
                                       toString(IDuniqueTable$Total[r]), "\n"),
                                 "FORWARD AND AMPLICON:", alignForward[1],
                                 "REVERSE AND AMPLICON:", alignReverse[1])
        }
        alignForward <- gsub("\n", "", alignForward)
        alignReverse <- gsub("\n", "", alignReverse)
        if (write_alignments >= 1) {
          algn_file_con <- c(algn_file_con,
                             paste("ID:", currentID,
                                   "Count:",
                                   toString(IDuniqueTable$Total[r])),
                             alignForward[4],
                             alignReverse[4])
        }

        forwardData <- getEventInfo(alignForward[3],
                                    currentID, "+",
                                    alignForward[4])
        reverseData <- getEventInfo(alignReverse[3],
                                    currentID, "-",
                                    alignReverse[4])

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

        # Filter our deletions on end and beginings
        forwardDataFiltered <-
          forwardData[!(forwardData$type == "deletion" &
                          end(forwardData) > nchar(amplicon))]
        reverseDataFiltered <-
          reverseData[!(reverseData$type == "deletion" &
                          start(reverseData) == 1)]
        # Filter insertions on end and beginings
        forwardDataFiltered <-
          forwardData[!(forwardData$type == "insertion" &
                          start(forwardData) > nchar(amplicon))]
        reverseDataFiltered <-
          reverseData[!(reverseData$type == "insertion" &
                          start(reverseData) == 1)]
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
        overlapFd <- subsetByOverlaps(ranges(forwardDataFiltered[forwardDataFiltered$type == "deletion"]), cutSites)
        overlapRe <- subsetByOverlaps(ranges(reverseDataFiltered[reverseDataFiltered$type == "deletion"]), cutSites)
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


    # Alignment verbosity levels
    if (write_alignments >= 1) {
      algn_file <- file.path(currentIDFolderName, "alignments.txt")
      if (file.exists(algn_file)) {
        file.remove(algn_file)
      }
      masterFileConn <- file(algn_file, open = "at")
      writeLines(algn_file_con, masterFileConn)
      close(masterFileConn)
    }
    if (write_alignments >= 2) {
      algn_detail_file <- file.path(currentIDFolderName, "detailed_alignments.txt")
      if (file.exists(algn_detail_file)) {
        file.remove(algn_detail_file)
      }
      uberAlignmentFD <- file(algn_detail_file, open = "at")
      writeLines(algn_det_file_con, uberAlignmentFD)
      close(uberAlignmentFD)
    }

    alignmentRangesBar <- c(alignmentRangesBar, alignmentRanges)
  }

  if (length(alignmentRangesBar) > 0) {
    utils::write.csv(as.data.frame(alignmentRangesBar),
                     file.path(resultsFolder,
                               paste0(barcode, "_alignment_ranges.csv")),
                     row.names = FALSE)
  }

  # Warnings
  sublog_file <- file.path(resultsFolder, paste0(currentID, "_SUBLOG.txt"))
  if (file.exists(sublog_file)) {
    file.remove(sublog_file)
  }
  if (!is.null(sublog_file_con)) {
    logFileConn <- file(sublog_file, open = "at")
    writeLines(sublog_file_con, logFileConn)
    close(logFileConn)
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
