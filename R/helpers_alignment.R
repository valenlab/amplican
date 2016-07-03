#' This function takes a lite string and gives back the events coordinates
#'
#' For the given example:
#' 0         1         2         3         4         5
#' 012345678901234567890123456789012345678901234567890123456789
#' NNNNNNNNNNN-----NNNNNNN-------NNNNNNNN--------NNNNNNNN------    Subject
#' ------NNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNN---NNNNNNNNNN    Pattern
#'
#' The RCPP library returns the information regarding the start and ends of the events in the
#' following way:
#' (Notice that there is a +1 shift in the index)
#'
#' #returned: @4@12,16*24,30*39,46*55,60*!@3@1,6*31,33*48,50*!@0@
#' #example: @0@!@1@152,207*!@1@25,t,A*
#' @param liteString (string)
#' @param ID (string)
#' @param strand (string) Either "+", "-" or default "*"
#' @return (GRanges) Object with metadata for insertion, deletion, missmatch
#' @import GenomicRanges S4Vectors
#'
getEventInfo <- function(liteString, ID, strand = "*"){

    if(is.na(liteString) | nchar(liteString) < 12){
      return(GRanges())
    }

    alignmentData <- unlist(strsplit(liteString, "!", fixed = TRUE))
    insertionsData <- alignmentData[1]
    deletionsData <- alignmentData[2]
    missmatchesData <- alignmentData[3]

    insertionsR <- GRanges()
    insCount <- as.numeric(unlist(strsplit(insertionsData, "@", fixed = TRUE))[2])
    if(insCount > 0){
      insertions <- unlist(strsplit(insertionsData, "@", fixed = TRUE))[3]
      insertions <- unlist(strsplit(insertions, "*", fixed = TRUE))
      insertionsPairs <- as.numeric(unlist(strsplit(insertions, ",", fixed = TRUE)))
      insertionsR <- GRanges(ranges = IRanges(start = insertionsPairs[c(T, F)], end = insertionsPairs[c(F, T)]),
                            strand = S4Vectors::Rle(rep(strand, insCount)),
                            seqnames = Rle(rep(ID, insCount)))
      insertionsR$mm_originally = ""
      insertionsR$mm_replacement = ""
      insertionsR$type = "insertion"
    }

    deletionsR <- GRanges()
    delCount <- as.numeric(unlist(strsplit(deletionsData , "@", fixed = TRUE))[2])
    if(delCount > 0){
      deletions <- unlist(strsplit(deletionsData, "@", fixed = TRUE))[3]
      deletions <- unlist(strsplit(deletions, "*", fixed = TRUE))
      deletionsPairs <- as.numeric(unlist(strsplit(deletions, ",", fixed = TRUE)))
      deletionsR <- GRanges(ranges = IRanges(start = deletionsPairs[c(T, F)], end = deletionsPairs[c(F, T)]),
                             strand = Rle(rep(strand, delCount)),
                             seqnames = Rle(rep(ID, delCount)))
      deletionsR$mm_originally = ""
      deletionsR$mm_replacement = ""
      deletionsR$type = "deletion"
    }

    missmatchesR <- GRanges()
    mmCount <- as.numeric(unlist(strsplit(missmatchesData, "@", fixed = TRUE))[2])
    if(mmCount > 0){
      missmatches <- unlist(strsplit(missmatchesData, "@", fixed = TRUE))[3]
      missmatches <- unlist(strsplit(missmatches, "*", fixed = TRUE))
      missmatches <- unlist(strsplit(missmatches, ",", fixed = TRUE))
      missmatchesR <- GRanges(ranges = IRanges(start = as.numeric(missmatches[c(T, F, F)]), width = 1),
                              strand = Rle(rep(strand, mmCount)),
                              seqnames = Rle(rep(ID, mmCount)))
      missmatchesR$mm_originally = missmatches[c(F, T, F)]
      missmatchesR$mm_replacement = missmatches[c(F, F, T)]
      missmatchesR$type = "missmatch"
    }

    return(c(insertionsR, deletionsR, missmatchesR))
}


#' For a given string, detect how many groups of uppercases is inside, where are
#' they, and how long they are.
#'
#' For example:
#'    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
#'
#' Has 4 groups of uppercases of length 7, 4, 1 and 3.
#' @param candidate (string) A string with the nucleotide sequence.
#' @import IRanges seqinr
#' @return (Ranges) A Ranges object with uppercases groups for given candidate string
#'
upperGroups <- function(candidate){
  return(reduce(IRanges(start = c(which(grepl("[[:upper:]]", s2c(candidate)))),
                        width = 1)))
}


#' Make alignments
#'
#' @param skip_bad_nucleotides (bool)
#' @param average_quality (int)
#' @param min_quality (int)
#' @param write_alignments (bool)
#' @param scoring_matrix (char)
#' @param gap_opening (bool)
#' @param gap_extension (int)
#' @param gap_ending (bool)
#' @param far_indels (bool)
#' @param configTable (bool)
#' @param resultsFolder (char)
#' @param fastqfiles (char)
#' @param PRIMER_DIMER (numeric)
#' @param cut_buffer (numeric)
#' @import ShortRead seqinr GenomicRanges
#' @importFrom utils write.table read.table
#' @return GRanges object with insertions, deletions and missmatches for each of the ID
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
                          PRIMER_DIMER = 10,
                          cut_buffer = 5){

  alignmentRanges <- GRanges()
  barcode <- toString(configTable$Barcode[1])
  message(paste0("Aligning reads for ", barcode))

  #Read Reads for this Barcode
  forwardsTable <- if (fastqfiles == 2) NULL else readFastq(configTable$Forward_Reads_File[1])
  reversesTable <- if (fastqfiles == 1) NULL else readFastq(configTable$Reverse_Reads_File[1])
  if (fastqfiles == 1) { rewersesTable <- rep(T, length(forwardsTable)) }
  if (fastqfiles == 2) { forwardsTable <- rep(T, length(reversesTable)) }

  #Filter Reads
  goodq <- goodBaseQuality(forwardsTable, min = min_quality) & goodBaseQuality(reversesTable, min = min_quality)
  avrq <- goodAvgQuality(forwardsTable, avg = average_quality) & goodAvgQuality(reversesTable, avg = average_quality)
  nucq <- alphabetQuality(forwardsTable) & alphabetQuality(reversesTable)
  goodReads <- goodq & avrq & nucq

  barcodeTable <- data.frame(barcode = barcode,
                             read_count = length(goodReads),
                             bad_base_quality = sum(!goodq),
                             bad_average_quality = sum(!avrq),
                             bad_alphabet = sum(!nucq),
                             filtered_read_count = sum(goodReads))

  forwardsTable <- forwardsTable[goodReads]
  reversesTable <- reversesTable[goodReads]

  #Unique reads
  uniqueTable <- data.frame( if (fastqfiles == 2) "" else as.character(forwardsTable@sread),
                             if (fastqfiles == 1) "" else as.character(reversesTable@sread))
  colnames(uniqueTable) <- c("Forward", "Reverse")
  uniqueTable$Total <- paste0(uniqueTable$Forward, uniqueTable$Reverse)
  uniqueTable <- aggregate(Total ~  Reverse + Forward, uniqueTable, length)
  uniqueTable$BarcodeFrequency <- uniqueTable$Total / sum(uniqueTable$Total)
  uniqueTable <- uniqueTable[order(uniqueTable$Forward, uniqueTable$Reverse),]
  uniqueTable$Asigned <- F
  uniqueTable$PRIMER_DIMER <- F

  barcodeTable$unique_reads <- nrow(uniqueTable)

  for(i in 1:dim(configTable)[1]){ #for each id
    currentID <- configTable[i, "ID"]
    # Primers and amplicon info
    forwardPrimer <- toString(configTable[i, "Forward_Primer"])
    reversePrimer <- toString(configTable[i, "Reverse_Primer"])
    guideRNA <- toString(configTable[i, "Target_Primer"])
    if (configTable[i, "Strand"] == 1) {
      guideRNA  <-  c2s(rev(comp(s2c(guideRNA), forceToLower = F, ambiguous = T)))
    }
    amplicon <- toString(configTable[i, "Amplicon"])
    # Names of files and folders
    currentIDFolderName <- paste0(resultsFolder, "/", currentID, "_", barcode)
    if (!dir.exists(currentIDFolderName)) {
      dir.create(file.path(currentIDFolderName), showWarnings = F)
    }

    # Alignment verbosity levels
    if (write_alignments >= 1) {
      algn_file <- paste0(currentIDFolderName, "/alignments.txt")
      if (file.exists(algn_file)) {file.remove(algn_file)}
      masterFileConn <- file(algn_file, open = "at")
    }
    if (write_alignments >= 2) {
      algn_detail_file <- paste0(currentIDFolderName, "/detailed_alignments.txt")
      if (file.exists(algn_detail_file)) {file.remove(algn_detail_file)}
      uberAlignmentFD <- file(paste0(currentIDFolderName, "/detailed_alignments.txt"), open = "at")
    }

    #Warnings
    sublog_file <- paste0(resultsFolder, "/", barcode, "_SUBLOG.txt")
    if (file.exists(sublog_file)) {file.remove(sublog_file)}
    logFileConn <- file(sublog_file, open = "at")
    configTable$Found_Guide[i]  <- checkTarget(guideRNA, amplicon, currentID, barcode, logFileConn)
    configTable$Found_Primers[i] <- checkPrimers(forwardPrimer, c2s(rev(comp(s2c(reversePrimer)))),
                                                 amplicon, currentID, barcode, logFileConn)
    cutSites <- upperGroups(amplicon) + cut_buffer
    if (length(cutSites) == 0) {
      message("Warning: Aligment position was not found in the amplicon. Find more information in the log file.")
      writeLines(paste0("Couldn't find upper case groups (PAM) position in amplicon.",
                        "/nFor ID: ", currentID, " and barcode: ", barcode, "/n"), logFileConn)
      configTable$Found_PAM[i] <- 0
      cutSites <- IRanges(start = 1, width = nchar(amplicon))
    }
    close(logFileConn)

    # Search for the forward, reverse and targets
    uniqueTable$forwardFound <- if (fastqfiles == 2) {F} else grepl(forwardPrimer, uniqueTable$Forward,
                                                                    ignore.case=TRUE)
    uniqueTable$reverseFound <- if (fastqfiles == 1) {F} else grepl(reversePrimer, uniqueTable$Reverse,
                                                                    ignore.case=TRUE)
    uniqueTable$guideFoundForward <- grepl(guideRNA, uniqueTable$Forward, ignore.case=TRUE)
    uniqueTable$guideFoundReverse <- grepl(c2s(rev(comp(s2c(guideRNA)))), uniqueTable$Reverse, ignore.case=TRUE)
    primersFound <- uniqueTable$forwardFound & uniqueTable$reverseFound
    if (fastqfiles == 1) { primersFound <- uniqueTable$forwardFound }
    if (fastqfiles == 2) { primersFound <- uniqueTable$reverseFound }
    uniqueTable$Asigned <- uniqueTable$Asigned | primersFound
    IDuniqueTable <- uniqueTable[primersFound, ]

    if (dim(IDuniqueTable)[1] > 0) {
      for (r in 1:dim(IDuniqueTable)[1]) {
        forwardString   <- toupper(toString(IDuniqueTable[r, "Forward"]))
        reverseString   <- toupper(c2s(rev(comp(s2c(toString(IDuniqueTable[r, "Reverse"]))))))
        # The return of gotoh is a string divided in five parts with the separator "++++"
        # -- [1] is the verbose alignment result
        # -- [2] is a string with the insertions and deletions and missmatches
        # -- [3] is a string with the insertions and deletions and missmatches,
        #    but from the subject (amplicon) coordinates
        # -- [4] is the alignment of the pattern.
        # -- [5] is the alignment of the subject.
        alignForward <- c("", "", "", "", "")
        if (fastqfiles != 2) {
          alignForward <- gRCPP(forwardString,
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
          alignReverse <- gRCPP(reverseString,
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
          writeLines(c(paste("ID:", currentID, "Count:", toString(IDuniqueTable$Total[r]), "\n"),
                       "FORWARD AND AMPLICON:", alignForward[1],
                       "REVERSE AND AMPLICON:", alignReverse[1]), uberAlignmentFD)
        }
        alignForward <- gsub("\n", "", alignForward)
        alignReverse <- gsub("\n", "", alignReverse)
        if (write_alignments >= 1) {
          writeLines(c(paste("ID:", currentID, "Count:", toString(IDuniqueTable$Total[r])),
                       alignForward[4],
                       alignReverse[4]), masterFileConn)
        }

        forwardData <- getEventInfo(alignForward[3], currentID)
        reverseData <- getEventInfo(alignReverse[3], currentID)

        #Filter PRIMER DIMERS and sum how many
        PD_cutoff <- nchar(amplicon) - (nchar(forwardPrimer) + nchar(reversePrimer) + PRIMER_DIMER)
        isPD <- any(c(width(forwardData), width(reverseData)) > PD_cutoff)
        IDuniqueTable[r, "PRIMER_DIMER"] <- isPD
        configTable$PRIMER_DIMER[i] <- configTable$PRIMER_DIMER[i] + isPD
        if (isPD) { next }

        #Frameshift table
        configTable$Frameshift_Forward[i] <- configTable$Frameshift_Forward[i] +
          sum(width(forwardData[forwardData$type != "missmatch"])) %% 3 != 0
        configTable$Frameshift_Reverse[i] <- configTable$Frameshift_Reverse[i] +
          sum(width(reverseData[reverseData$type != "missmatch"])) %% 3 != 0

        #Cut definition
        forwardData$cut <- F
        reverseData$cut <- F
        overlapFd <- subsetByOverlaps(ranges(forwardData[forwardData$type == "deletion"]), cutSites)
        overlapRe <- subsetByOverlaps(ranges(reverseData[reverseData$type == "deletion"]), cutSites)
        if (fastqfiles == 0) {
          #forward and reverse have to agree on deletion
          overlapFd <- overlapFd[!is.na(match(overlapFd, overlapRe))]
          overlapRe <- overlapRe[!is.na(match(overlapRe, overlapFd))]
          forwardData$cut[!is.na(match(ranges(forwardData), overlapFd))] <- T
          reverseData$cut[!is.na(match(ranges(reverseData), overlapRe))] <- T
        } else if (fastqfiles == 1) {
          forwardData$cut[!is.na(match(ranges(forwardData), overlapFd))] <- T
        } else {
          reverseData$cut[!is.na(match(ranges(reverseData), overlapRe))] <- T
        }

        #Strand Count Frequency
        strand(forwardData) <- "+"
        forwardData$count <- IDuniqueTable$Total[r]
        forwardData$frequency <- IDuniqueTable$Total[r]/sum(IDuniqueTable$Total)

        strand(reverseData) <- "-"
        reverseData$count <- IDuniqueTable$Total[r]
        reverseData$frequency <- IDuniqueTable$Total[r]/sum(IDuniqueTable$Total)

        alignmentRanges <- c(alignmentRanges, forwardData, reverseData)

        # events counts
        configTable$Cut_Forward[i] <- configTable$Cut_Forward[i] +
          any(forwardData$cut) * IDuniqueTable$Total[r]
        configTable$Cut_Reverse[i] <- configTable$Cut_Reverse[i] +
          any(reverseData$cut) * IDuniqueTable$Total[r]
      }
      configTable$Reads[i] <- sum(IDuniqueTable$Total)
    }

    if (isOpen(masterFileConn)) {
      close(masterFileConn)
    }
    if (isOpen(uberAlignmentFD)) {
      close(uberAlignmentFD)
    }
  }

  barcodeTable$unassigned_reads <- sum(!uniqueTable$Asigned)
  barcodeTable$assigned_reads <- sum(uniqueTable$Asigned)
  write.table(uniqueTable[!uniqueTable$Asigned, ],
              file = paste(resultsFolder, "/unassigned_sequences/", barcode, "_unassigned_reads.txt", sep = ''),
              quote = FALSE, sep = "\t", row.names = F)

  write.table(configTable, paste0(resultsFolder, "/", barcode, "_configFile_results") , sep="\t", row.names = F)
  write.table(barcodeTable, paste0(resultsFolder, "/", barcode, "_reads_filters.csv"), sep="\t", row.names = F)
  return(alignmentRanges)
}
