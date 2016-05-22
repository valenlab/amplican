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
#' @param liteString (list)
#' @import GenomicRanges
#' @return (GRanges) Object with metadata for insertion, deletion, missmatch
getEventInfo <- function(liteString, ID, strand = "*"){

    if(is.na(liteString) | length(liteString) < 12){
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
                            strand = Rle(rep(strand, insCount)),
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

#' For a given string, detect how many groups of uppercases are there, where are
#' they, and how long they are.
#'
#' For example:
#'    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
#'
#' Has 4 groups of uppercases of length 7, 4, 1 and 3.
#' @param candidate (String) A string with the nucleotide sequence.
#' @import IRanges
#' @return (Ranges) A Ranges object with uppercases groups for given candidate string
#'
upperGroups <- function(candidate){
  return(reduce(IRanges(start=c(which(grepl("[[:upper:]]", s2c(candidate)))), width = 1)))
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
#' @param alignmentFolder (char)
#' @param procesID (int)
#' @param temp_folder (char)
#' @param fastqfiles (char)
#' @import ShortRead, seqinr, GenomicRanges
#' @return no clue
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
                          far_indels = TRUE){

  barcode <- configTable$Barcode[1]

  #Read Reads for this Barcode
  forwardsTable <- readFastq(configTable$Forward_Reads_File[1])
  reversesTable <- readFastq(configTable$Reverse_Reads_File[1])

  #Filter Reads
  goodq <- goodBaseQuality(forwardsTable, min = min_quality) & goodBaseQuality(reversesTable, min = min_quality)
  avrq <- goodAvgQuality(forwardsTable, avg = average_quality) & goodAvgQuality(reversesTable, avg = average_quality)
  nucq <- alphabetQuality(forwardsTable) & alphabetQuality(reversesTable)
  goodReads <- goodq & avrq & nucq

  barcodeTable <- data.frame(Total_Pre_Filter = length(goodReads), Total_Post_Filter = sum(goodReads))

  forwardsTable <- forwardsTable[goodReads]
  reversesTable <- reversesTable[goodReads]

  #Unique reads
  uniqueTable <- data.frame(as.character(forwardsTable@sread), as.character(reversesTable@sread))
  colnames(uniqueTable) <- c("Forward","Reverse")
  uniqueTable$Total <- paste0(uniqueTable$Forward, uniqueTable$Reverse)
  uniqueTable <- aggregate(Total ~  Reverse + Forward, uniqueTable, length)
  uniqueTable$Frequency <- uniqueTable$Total / sum(uniqueTable$Total)
  uniqueTable <- uniqueTable[order(uniqueTable$Forward, uniqueTable$Reverse),]
  uniqueTable$asigned <- F

  barcodeTable$Experiment_Unique_Sequences <- nrow(uniqueTable)
  barcodeTable$Total_Experiment_Sequences  <- sum(uniqueTable$Total)

  for(i in dim(configTable)[1]){ #for each id
    currentID <- configTable[i, "ID"]
    # Primers and amplicon info
    forwardPrimer <- toString(configTable[i, "Forward_Primer"])
    reversePrimer <- toString(configTable[i, "Reverse_Primer"])
    targetPrimer <- toString(configTable[i, "Target_Primer"])
    # Update the target Primer if necessary
    reverseAmplicon <- configTable[i, "Strand"]
    if (reverseAmplicon == 1) {
      targetPrimer  <-  c2s(rev(comp(s2c(targetPrimer), forceToLower = F, ambiguous = T)))
    }
    amplicon <- toString(configTable[i, "Amplicon"])
    # Names of files and folders
    currentIDFolderName <- paste0(resultsFolder,"/", currentID, "_", barcode)
    if(!dir.exists(currentIDFolderName)) {
      dir.create(file.path(currentIDFolderName), showWarnings = TRUE)
    }
    masterAlignmentFilePath <- paste0(currentIDFolderName, "/alignments.txt")
    uberAlignmentFilePath <- paste0(currentIDFolderName, "/verbose.txt")
    deletionFilePath <- paste0(currentIDFolderName, "/", currentID, "_", barcode, "_deletions.txt")
    deletionRelativeFilePath <- paste0(currentIDFolderName, "/", currentID, "_", barcode, "_deletionsRelative.txt")
    uniqueTablePath <- paste0(currentIDFolderName, "/", currentID, "_", barcode, "_uniques.txt")

    cutSites <- upperGroups(amplicon)

    # If we need to write the alignment, prepare the proper file descriptors
    if (write_alignments >= 1) {
      # Prepare the summary alignment file
      masterFileConn <- file(masterAlignmentFilePath, open = "at")
      writeLines(c(toString(amplicon), "\n"), masterFileConn)
    }
    if (write_alignments >= 2) {
      # Start the file where the alignments are accumulated
      uberAlignmentFD <- file(uberAlignmentFilePath, open = "at")
    }

    # Prepare the unique table where all the sequences are going.
    uniqueTableFD <- file(uniqueTablePath, open = "at")
    writeLines(paste("Total", "Frequency", "Forward", "Reverse", "Forward_Alignment_String",
                     "Reverse_Alignment_String", "Forward_Alignment_Genome_Coordinates_String",
                     "Reverse_Alignment_Genome_Coordinates_String", "Forward_Alignment_Length",
                     "Reverse_Alignment_Length", "Alignment", "Forward_Alignment_Positions",
                     "Reverse_Alignment_Positions", "Forward_Wide_Position", "Reverse_Wide_Position",
                     "Target_Forward_Found", "Target_Reverse_Found", "Insertions_Forward", "Insertions_Reverse",
                     "Deletions_Forward", "Deletions_Reverse", "Missmatches_Forward", "Missmatches_Reverse",
                     sep = '\t'), uniqueTableFD)

    #Warnings
    logFileConn <- file(paste0(resultsFolder, "/", barcode, "_SUBLOG.txt"), open = "at")
    configTable$Found_Target[i]  <- checkTarget(targetPrimer, amplicon, currentID, barcode, logFileConn)
    configTable$Found_Primers[i] <- checkPrimers(forwardPrimer, c2s(rev(comp(s2c(reversePrimer)))), amplicon, currentID, barcode, logFileConn)
    configTable$Found_AP[i] <- checkPositions(start(cutSites), amplicon, currentID, barcode, logFileConn)
    close(logFileConn)

    # Search for the forward, reverse and targets
    uniqueTable$forwardFound <- grepl(forwardPrimer, uniqueTable$Forward, ignore.case=TRUE)
    uniqueTable$reverseFound <- grepl(reversePrimer, uniqueTable$Reverse, ignore.case=TRUE)
    uniqueTable$targetFoundForward <- grepl(targetPrimer, uniqueTable$Forward, ignore.case=TRUE)
    uniqueTable$targetFoundReverse <- grepl(c2s(rev(comp(s2c(targetPrimer)))), uniqueTable$Reverse, ignore.case=TRUE)
    primersFound <- uniqueTable$forwardFound & uniqueTable$reverseFound
    if (fastqfiles == 1) {
      primersFound <- uniqueTable$forwardFound
    }
    if (fastqfiles == 2) {
      primersFound <- uniqueTable$reverseFound
    }
    uniqueTable$asigned <- uniqueTable$asigned | primersFound
    IDuniqueTable <- uniqueTable[primersFound, ]

    for (r in dim(IDuniqueTable)[1]) {
      forwardString   <- toupper(toString(IDuniqueTable[r, "Forward"]))
      reverseString   <- toupper(c2s(rev(comp(s2c(toString(IDuniqueTable[r, "Reverse"]))))))
      # The return of gotoh is a string divided in five parts with the separator "++++"
      # -- [1] is the verbose alignment result
      # -- [2] is a string with the insertions and deletions and missmatches
      # -- [3] is a string with the insertions and deletions and missmatches,
      #    but from the subject (amplicon) coordinates
      # -- [4] is the alignment of the pattern.
      # -- [5] is the alignment of the subject.
      alignForward <- ""
      if (fastqfiles == 0 || fastqfiles == 1) {
        alignForward <- gRCPP(forwardString,
                              amplicon,
                              scoring_matrix,
                              gap_opening,
                              gap_extension,
                              gap_ending,
                              far_indels)
        alignForward <- gsub("\n", "", unlist(strsplit(alignForward, "++++", fixed = TRUE)))
      }
      alignReverse <- ""
      if (fastqfiles == 0 || fastqfiles == 2) {
        alignReverse <- gRCPP(reverseString,
                              amplicon,
                              scoring_matrix,
                              gap_opening,
                              gap_extension,
                              gap_ending,
                              far_indels)
        alignReverse <- gsub("\n", "", unlist(strsplit(alignReverse, "++++", fixed = TRUE)))
      }
      # Write the alignments
      if (write_alignments >= 1) {
        writeLines(c(paste(currentID, toString(uniqueTable$Total[r])),
                     alignForward[4],
                     alignReverse[4]), masterFileConn)
      }
      if (write_alignments >= 2) {
        writeLines(c(paste0(currentID, "_", r), "\n FORWARD AND AMPLICON: \n",
                     alignForward[1], "\n REVERSE AND AMPLICON: \n",
                     alignReverse[1], "\n"), uberAlignmentFD)
      }

      forwardData <- getEventInfo(alignForward[2], currentID)
      reverseData <- getEventInfo(alignReverse[2], currentID)

      # Write the sequence info
      forwardAllPositions <- upperGroups(toString(alignForward[5]))
      reverseAllPositions <- upperGroups(toString(alignReverse[5]))

      writeLines(paste(as.numeric(r),
                       as.numeric(IDuniqueTable$Total[r]),
                       as.numeric(IDuniqueTable$Frequency[r]),
                       IDuniqueTable$Forward[r],
                       IDuniqueTable$Reverse[r],
                       gsub("\n", "", alignForward[2]),
                       gsub("\n", "", alignReverse[2]),
                       gsub("\n", "", alignForward[3]),
                       gsub("\n", "", alignReverse[3]),
                       as.numeric(ifelse(fastqfiles == 0 || fastqfiles == 1, nchar(alignForward[4]), 0)),
                       as.numeric(ifelse(fastqfiles == 0 || fastqfiles == 2, nchar(alignReverse[4]), 0)),
                       paste0(currentID, ".txt"),
                       start(forwardAllPositions),
                       start(reverseAllPositions),
                       width(forwardAllPositions),
                       width(reverseAllPositions),
                       IDuniqueTable$ForwardFound[r],
                       IDuniqueTable$ReverseFound[r],
                       sum(forwardData$type == "insertion"),
                       sum(reverseData$type == "insertion"),
                       sum(forwardData$type == "deletion"),
                       sum(reverseData$type == "deletion"),
                       sum(forwardData$type == "missmatch"),
                       sum(forwardData$type == "missmatch"),
                       sep = '\t'), uniqueTableFD)

      configTable$Sum_Insertions_Forward[i] <- configTable$Sum_Insertions_Forward[i] +
                                               sum(forwardData$type == "insertion") * IDuniqueTable$Total[i]
      configTable$Sum_Insertions_Reverse[i] <- configTable$Sum_Insertions_Reverse[i] +
                                               sum(reverseData$type == "insertion") * IDuniqueTable$Total[i]
      configTable$Sum_Deletions_Forward[i] <- configTable$Sum_Deletions_Forward[i] +
                                              sum(forwardData$type == "deletion") * IDuniqueTable$Total[i]
      configTable$Sum_Deletions_Reverse[i] <- configTable$Sum_Deletions_Reverse[i] +
                                              sum(reverseData$type == "deletion") * IDuniqueTable$Total[i]
      configTable$Sum_Missmatches_Forward[i] <- configTable$Sum_Missmatches_Forward[i] +
                                                sum(forwardData$type == "missmatch") * IDuniqueTable$Total[i]
      configTable$Sum_Missmatches_Reverse[i] <- configTable$Sum_Missmatches_Reverse[i] +
                                                sum(reverseData$type == "missmatch") * IDuniqueTable$Total[i]
    }
    configTable$Sum_Target_Forward_Found[i] <- sum(uniqueTable$forwardFound)
    configTable$Sum_Target_Reverse_Found[i] <- sum(uniqueTable$reverseFound)
  }

  write.table(uniqueTable[!uniqueTable$asigned, ],
              file = paste(resultsFolder, "/", barcode, "_unassigned.txt", sep = '', collapse = ''),
              quote = FALSE, sep = "\t")

  if (masterFileConnOpen) {
    close(masterFileConn)
  }
  if (uberAlignmentFDOpen) {
    close(uberAlignmentFD)
  }
  close(uniqueTableFD)

  write.table(configTable, paste0(alignmentFolder, "/", barcode, "_configFile_results") , sep="\t")
  write.table(barcodeTable, paste0(alignmentFolder, "/", barcode, "_barcodeFile_results"), sep="\t")
  return()
}
