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
  finalGR <- GRanges(
    ranges = x,
    strand = Rle(rep(strand_info, length(x))),
    seqnames = Rle(rep(ID, length(x))))

  finalGR$originally = originally
  finalGR$replacement = replacement
  finalGR$type = type

  return(finalGR)
}


#' This function takes alignments and gives back the events coordinates.
#'
#' @param align (PairwiseAlignmentsSingleSubject)
#' @param ID (string)
#' @param ampl_shift (numeric) Shift events additionaly by this value.
#' PairwiseAlignmentsSingleSubject returns truncated alignments.
#' @param strand_info (string) Either '+', '-' or default '*'
#' @return (GRanges) Object with meta-data for insertion, deletion, mismatch
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#'
getEventInfo <- function(align, ID, ampl_shift, strand_info = "+") {

  if (ampl_shift < 1) stop("Amplicon shift can't be less than 1.")

  del <- deletion(align)[[1]]
  ins <- insertion(align)[[1]]
  mm <- mismatchSummary(summary(align))$subject
  finalGR <- GRanges()

  if (length(ins) > 0) {
    ins_shift <- shift(ins, c(0, cumsum(width(ins))[-length(ins)]))
    ins_seq <-
      substr(rep(as.character(pattern(align)), length(ins)),
             start = start(ins_shift), stop = end(ins_shift))
    finalGR <- defGR(ins, ID, strand_info, "insertion", "", ins_seq)
  }

  if (length(del) > 0) {
    del <- shift(del, c(0, cumsum(width(del))[-length(del)]))
    # make deletions to be relative to the subject
    if (length(ins) > 0) { # shift when insertions are before deletions
      for (i in seq_along(del)) {
        ins_before <- start(del)[i] > start(ins)
        if (any(ins_before)) {
          del[i] <- shift(del[i], -1 * sum(width(ins[ins_before])))
        }
      }
    }
    finalGR <- c(finalGR, defGR(del, ID, strand_info))
  }

  if (dim(mm)[1] > 0) {
    finalGR <- c(finalGR,
                 defGR(IRanges(mm$SubjectPosition, width = 1),
                       ID,
                       strand_info,
                       "mismatch",
                       as.character(mm$Subject),
                       as.character(mm$Pattern)))
  }

  # for some reason mismatches are relative to the full subject sequence
  # insertions and deletions are not
  ins_del <- finalGR$type %in% c("insertion", "deletion")
  if (strand_info == "+") {
    finalGR[ins_del] <- shift(finalGR[ins_del], ampl_shift - 1)
  } else {
    subj_bas <- stringr::str_count(as.character(subject(align)), "[ATCG]")
    finalGR[ins_del] <- shift(finalGR[ins_del], ampl_shift - subj_bas)
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
    start = which(stringr::str_detect(strsplit(candidate, "")[[1]],
                                      "[[:upper:]]")),
    width = 1
  )))
}


#' Reverse complement events that have amplicons with direction 1.
#'
#' @param idR (data.frame) Loaded events.
#' @param cfgT (data.frame) Loaded configuration file.
#' @return (data.frame) Returns input idR, but events for amplicons with
#' direction 1 reverse complemented, "+" and "-" swapped.
#'
flipRanges <- function(idR, cfgT) {

  is_dir <- as.logical(cfgT$Direction)
  to_flip <- cfgT[is_dir, "ID"]
  to_flip <- idR$seqnames %in% to_flip

  if (any(to_flip)) {
    ampl_lengths <- nchar(as.character(cfgT[is_dir, "Amplicon"]))
    ampl_ids <- as.character(cfgT[is_dir, "ID"])

    ids_mapping <- match(idR[to_flip, "seqnames"], ampl_ids)
    ampl_lengths <- ampl_lengths[ids_mapping]

    idR[to_flip, "originally"] <- revComp(idR[to_flip, "originally"])
    idR[to_flip, "replacement"] <- revComp(idR[to_flip, "replacement"])

    old_starts <- idR[to_flip, "start"]
    idR[to_flip, "start"] <- ampl_lengths - idR[to_flip, "end"] + 1
    idR[to_flip, "end"] <- ampl_lengths - old_starts + 1

    strand <- idR[to_flip, "strand"]
    strand_minus <-strand == "-"
    strand[strand == "+"] <- "-"
    strand[strand_minus] <- "+"
    idR[to_flip, "strand"] <- strand
  }
  return(idR)
}


#' Make alignments helper.
#'
#' Main functionality of the package, aligning reads to the amplicon.
#' @param cfgT config file as data table
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
makeAlignment <- function(cfgT,
                          resultsFolder,
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
  barcode <- cfgT$Barcode[1]
  almRBar <- GRanges()
  GenomeInfoDb::seqlevels(almRBar) <- unique(cfgT$ID)

  message(paste0("Aligning reads for ", barcode))

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
                             experiment_count = cfgT$ExperimentsCount[1],
                             read_count = length(goodReads),
                             bad_base_quality = sum(!goodq),
                             bad_average_quality = sum(!avrq),
                             bad_alphabet = sum(!nucq),
                             filtered_read_count = sum(goodReads))

  fwdT <- fwdT[goodReads]
  rveT <- rveT[goodReads]

  # Unique reads
  unqT <- data.frame(
    if (fastqfiles == 2) "" else as.character(sread(fwdT)),
    if (fastqfiles == 1) "" else as.character(sread(rveT)))
  colnames(unqT) <- c("Forward", "Reverse")
  unqT$Total <- paste0(unqT$Forward, unqT$Reverse)
  unqT <- stats::aggregate(Total ~ Forward + Reverse, unqT, length)
  unqT$BarcodeFrequency <- unqT$Total / sum(unqT$Total)
  unqT <- unqT[order(unqT$Forward, unqT$Reverse), ]
  unqT$Asigned <- FALSE
  unqT$PRIMER_DIMER <- FALSE
  unqT[c("Forward", "Reverse")] <- lapply(unqT[c("Forward", "Reverse")],
                                          function(x) toupper(as.character(x)))
  barcodeTable$unique_reads <- nrow(unqT)

  # for each experiment
  for (i in seq_len(dim(cfgT)[1])) {

    almR <- GRanges()
    thisID <- cfgT[i, "ID"]
    # Primers and amplicon info
    fwdPrimer <- toupper(cfgT[i, "Forward_Primer"])
    rvePrimer <- toupper(cfgT[i, "Reverse_Primer"])
    gRNA <- toupper(cfgT[i, "guideRNA"])
    rgRNA <- toupper(cfgT[i, "RguideRNA"])
    amplicon <- cfgT[i, "Amplicon"]

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
    unqT$guideFoundForward <- stringr::str_detect(unqT$Forward, gRNA)
    unqT$guideFoundReverse <- stringr::str_detect(unqT$Reverse, rgRNA)
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

    if (dim(IDunqT)[1] > 0) {

      if (fastqfiles != 2) {
        alignFwd <- Biostrings::pairwiseAlignment(
          Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Forward"]),
                             start = IDunqT$fwdPrInReadPos),
          Biostrings::subseq(toupper(amplicon),
                             start = cfgT$fwdPrPos[i]),
          type = "global",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension)
      }
      if (fastqfiles != 1) {
        alignRve <- Biostrings::pairwiseAlignment(
          Biostrings::reverseComplement(
            Biostrings::subseq(Biostrings::DNAStringSet(IDunqT[, "Reverse"]),
                               start = IDunqT$rvePrInReadPos)),
          Biostrings::subseq(toupper(amplicon),
                             start = 1,
                             end = cfgT$rvePrPosEnd[i]),
          type = "global",
          substitutionMatrix =  scoring_matrix,
          gapOpening = gap_opening,
          gapExtension = gap_extension)
      }
      # Write the alignments
      if (write_alignments >= 1) {
        algn_file <-
          file.path(resultsFolder, paste0(thisID, "_", barcode, ".txt"))
        if (file.exists(algn_file)) {
          file.remove(algn_file)
        }
        algn_file_con <- file(algn_file, open = "at")
        writeLines(as.vector(rbind(
          paste("ID:", thisID,
                "Count:", format(IDunqT$Total)),
          if (fastqfiles != 2) rbind(as.character(pattern(alignFwd)),
                                     as.character(subject(alignFwd))),
          if (fastqfiles < 1) "",
          if (fastqfiles != 1) rbind(
            as.character(pattern(alignRve)),
            as.character(subject(alignRve))),
          ""
        )),  algn_file_con)
        close(algn_file_con)
      }

      for (r in seq_len(dim(IDunqT)[1])) {

        fwdD <- if (fastqfiles != 2) {
          getEventInfo(alignFwd[r], thisID, cfgT$fwdPrPos[i], "+")
        } else {
          GRanges()
        }

        rveD <- if (fastqfiles != 1) {
          getEventInfo(alignRve[r], thisID, cfgT$rvePrPosEnd[i], "-")
        } else {
          GRanges()
        }

        # Filter PRIMER DIMERS and sum how many
        PD_cutoff <- nchar(amplicon) -
          (nchar(fwdPrimer) + nchar(rvePrimer) + PRIMER_DIMER)
        isPD <- any(c(width(fwdD),
                      width(rveD)) > PD_cutoff)
        IDunqT[r, "PRIMER_DIMER"] <- isPD
        cfgT$PRIMER_DIMER[i] <- cfgT$PRIMER_DIMER[i] +
          isPD * IDunqT$Total[r]
        if (isPD) {
          next
        }

        # Filters
        fwdDFtd <-
          fwdD[!(fwdD$type %in% div &
                   end(fwdD) >= cfgT$rvePrimerPosition[i])]
        rveDFtd <-
          rveD[!(rveD$type %in% div &
                   end(rveD) >= cfgT$rvePrimerPosition[i])]
        fwdDFtd <-
          fwdD[!(fwdD$type %in% div &
                   start(fwdD) <= cfgT$fwdPrPosEnd[i])]
        rveDFtd <-
          rveD[!(rveD$type %in% div &
                   start(rveD) <= cfgT$fwdPrPosEnd[i])]

        # Frameshift table
        frameshift <- FALSE
        if (fastqfiles == 0 | fastqfiles == 0.5) {
          frameshift <-
            sum(width(fwdDFtd[fwdDFtd$type == "insertion" |
                                fwdDFtd$type == "deletion"])) %%
            3 !=  0 &
            sum(width(rveDFtd[rveDFtd$type == "insertion" |
                                rveDFtd$type == "deletion"])) %%
            3 != 0
        } else if (fastqfiles == 1) {
          frameshift <-
            sum(width(fwdDFtd[fwdDFtd$type == "insertion" |
                                fwdDFtd$type == "deletion"])) %%
            3 != 0
        } else {
          frameshift <-
            sum(width(rveD[rveDFtd$type == "insertion" |
                             rveDFtd$type == "deletion"])) %%
            3 != 0
        }
        cfgT$Frameshift[i] <- cfgT$Frameshift[i] +
          frameshift * IDunqT$Total[r]

        # definitions
        if (length(fwdD) > 0) {
          fwdD$cut <- FALSE
          fwdD$count <- IDunqT$Total[r]
          fwdD$frequency <- 0  #prepare field
          fwdD$read_id <- r
        }

        if (length(rveD) > 0) {
          rveD$cut <- FALSE
          rveD$count <- IDunqT$Total[r]
          rveD$frequency <- 0  #prepare field
          rveD$read_id <- r
        }

        # cut assessment
        overlapFd <- subsetByOverlaps(
          ranges(fwdDFtd[fwdDFtd$type == "deletion"]),
          cfgT$cutSites[[i]])
        overlapRe <- subsetByOverlaps(
          ranges(rveDFtd[rveDFtd$type == "deletion"]),
          cfgT$cutSites[[i]])
        if (fastqfiles == 0 | fastqfiles == 0.5) {
          # forward and reverse have to agree on deletion
          overlapFd <- overlapFd[!is.na(match(overlapFd, overlapRe))]
          overlapRe <- overlapRe[!is.na(match(overlapRe, overlapFd))]
          fwdD$cut[!is.na(match(ranges(fwdD), overlapFd))] <- TRUE
          rveD$cut[!is.na(match(ranges(rveD), overlapRe))] <- TRUE
        } else if (fastqfiles == 1) {
          fwdD$cut[!is.na(match(ranges(fwdD), overlapFd))] <- TRUE
        } else {
          rveD$cut[!is.na(match(ranges(rveD), overlapRe))] <- TRUE
        }

        cfgT$Cut[i] <- cfgT$Cut[i] +
          if (fastqfiles == 2) {
            any(rveD$cut) * IDunqT$Total[r]
          } else {
            any(fwdD$cut) * IDunqT$Total[r]
          }

        if (length(fwdD) > 0 | length(rveD) > 0) {
          almR <- c(almR,
                    fwdD,
                    rveD)
        }
      }

      cfgT$Reads[i] <- sum(IDunqT$Total)
      # fill frequency
      if (length(almR) > 0) {
        reads_count_filtered <- cfgT$Reads[i] - cfgT$PRIMER_DIMER[i]
        almR$frequency <- almR$count / reads_count_filtered
      }
    }

    almRBar <- c(almRBar, almR)
  }

  if (length(almRBar) > 0) {
    almRBar <- flipRanges(as.data.frame(almRBar),
                          cfgT)
    if (fastqfiles == 2 | fastqfiles == 1) {
      # when only one strand treat everything as on forward strand
      almRBar$strand <- "+"
    }
    utils::write.csv(almRBar,
                     file.path(resultsFolder,
                               paste0(barcode, "_alignment_ranges.csv")),
                     row.names = FALSE)
  }

  barcodeTable$unassigned_reads <- sum(!unqT$Asigned)
  barcodeTable$assigned_reads <- sum(unqT$Asigned)
  utils::write.csv(unqT[!unqT$Asigned, ],
                   file = file.path(resultsFolder,
                                    "unassigned_sequences",
                                    paste0(barcode, "_unassigned_reads.csv")),
                   quote = FALSE, row.names = FALSE)

  revDir <- cfgT[, "Direction"] == 1
  # revert guides to 5'-3'
  cfgT[revDir, "guideRNA"] <- revComp(cfgT[revDir, "guideRNA"])
  utils::write.csv(cfgT[c("ID",
                          "Barcode",
                          "Forward_Reads_File",
                          "Reverse_Reads_File",
                          "Group",
                          "guideRNA",
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
