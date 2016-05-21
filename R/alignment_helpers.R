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
#' (Please notice that there is a +1 shift in the index)
#'
#' #returned: @4@ 12,16 * 24,30 * 39,46 * 55,60 * ! @3@ 1,6 * 31,33  * 48,50 *
#' @param liteString (list)
#' @return A list with 6 elements:
#'     - The element 1 and element 4 is the total of insertions and deletions (in that order)
#'     - The element 2 and element 5 is where the insertions and deletions starts (in that order)
#'     - The element 3 and element 6 is where the insertions and deletion ends in (that order)
#'     - Elements 2,3,4,5 are integer arrays.
#'     - If the total of insertions/deletions is 0, the integer arrays will be equal to NULL
getEventInfo <- function(liteString){
  if(is.na(liteString) == FALSE){
    # Divide the string in two pieces; one for the insertions and the other one for the deletions
    alignmentData <- unlist(strsplit(liteString, "!", fixed = TRUE))
    alignmentInsertionsData <- alignmentData[1]
    alignmentDeletionsData  <- alignmentData[2]
    alignmentMissmatchData  <- alignmentData[3]

    # The totals of events are after the first @ and before the second @
    alignmentInsertionsTotal  <- as.numeric(unlist(strsplit(alignmentInsertionsData,"@",
                                                            fixed = TRUE))[2])
    alignmentDeletionsTotal   <- as.numeric(unlist(strsplit(alignmentDeletionsData ,"@",
                                                            fixed = TRUE))[2])
    alignmentMissmatchesTotal <- as.numeric(unlist(strsplit(alignmentMissmatchData ,"@",
                                                            fixed = TRUE))[2])

    # We need 7 arrays, 2 for the insertions, 2 for the deletions, 3 for the missmatches
    #
    # For deletions and insertions:
    # Each of those two arrays, will represent the starts and the ends of each event.
    # So for event 1, we need to look where does it start in the start_array[1] and where does
    # ends in end_array[1]
    #
    # For missmatches
    # One array represent the original nucleotide
    # One array represent the new      nucleotide
    # One array represent the position in coordinates

    # Create an array for the boundaries of the event insertions
    alignmentInsertionsStarts <- rep(0,alignmentInsertionsTotal)
    alignmentInsertionsEnds <-   rep(0,alignmentInsertionsTotal)

    # Create an array for the boundaries of the event deletions
    alignmentDeletionsStarts <- rep(0,alignmentDeletionsTotal)
    alignmentDeletionsEnds <-   rep(0,alignmentDeletionsTotal)

    # Create the arrays for the mutations
    alignmentMutationsOriginal <- rep('x', alignmentMissmatchesTotal)
    alignmentMutationsMutated  <- rep('x', alignmentMissmatchesTotal)
    alignmentMutationsPosition <- rep(0  , alignmentMissmatchesTotal)

    # Fill the arrays for the event insertions
    alignmentInsertionsBoundariesData  <- unlist(strsplit(alignmentInsertionsData,"@",
                                                          fixed = TRUE))[3]
    alignmentInsertionsBoundariesArray <- unlist(strsplit(alignmentInsertionsBoundariesData,"*",
                                                          fixed = TRUE))
    if (alignmentInsertionsTotal>0) {
      for(i in 1:alignmentInsertionsTotal){
        indelPair <- unlist(strsplit(alignmentInsertionsBoundariesArray[i],",",fixed = TRUE))
        alignmentInsertionsStarts[i] <- as.numeric(indelPair[1])
        alignmentInsertionsEnds[i]   <- as.numeric(indelPair[2])
      }
    } else {
      alignmentInsertionsStarts <- NULL
      alignmentInsertionsEnds   <- NULL
    }

    # Fill the arrays for the event deletions
    alignmentDeletionsBoundariesData <- unlist(strsplit(alignmentDeletionsData,"@",
                                                        fixed = TRUE))[3]
    alignmentDeletionsBoundariesArray <- unlist(strsplit(alignmentDeletionsBoundariesData,"*",
                                                         fixed = TRUE))
    if (alignmentDeletionsTotal>0) {
      for(l in 1:alignmentDeletionsTotal){
        indelPair <- unlist(strsplit(alignmentDeletionsBoundariesArray[l],",",fixed = TRUE))
        alignmentDeletionsStarts[l] <- as.numeric(indelPair[1])
        alignmentDeletionsEnds[l]   <- as.numeric(indelPair[2])
      }
    } else {
      alignmentDeletionsStarts <- NULL
      alignmentDeletionsEnds   <- NULL
    }

    # Fill the arrays for the missmatches
    alignmentMissmatchBoundariesData <- unlist(strsplit(alignmentMissmatchData,"@",
                                                        fixed = TRUE))[3]
    alignmentMissmatchBoundariesArray <- unlist(strsplit(alignmentMissmatchBoundariesData,"*",
                                                         fixed = TRUE))
    if (alignmentMissmatchesTotal>0) {
      for(l in 1:alignmentMissmatchesTotal){
        missMatchTrio <- unlist(strsplit(alignmentMissmatchBoundariesArray[l],",",fixed = TRUE))
        alignmentMutationsOriginal[l]  <- missMatchTrio[2]
        alignmentMutationsMutated[l]   <- missMatchTrio[3]
        alignmentMutationsPosition [l] <- as.numeric(missMatchTrio[1])

      }
    } else {
      alignmentMutationsOriginal <- NULL
      alignmentMutationsMutated  <- NULL
      alignmentMutationsPosition <- NULL
    }

    return(list(alignmentInsertionsTotal, alignmentInsertionsStarts, alignmentInsertionsEnds,
                 alignmentDeletionsTotal,  alignmentDeletionsStarts, alignmentDeletionsEnds,
                 alignmentMutationsOriginal, alignmentMutationsMutated, alignmentMutationsPosition))

  } else {
    return(list(0, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL))
  }
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
#' @import ShortRead, seqinr
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
    writeLines(paste("Total", "Frequency", "Forward", "Reverse", "Forward_Deletions_String",
                     "Reverse_Deletions_String", "Forward_Deletions_Genome_Coordinates_String",
                     "Reverse_Deletions_Genome_Coordinates_String", "Forward_Alignment_Length",
                     "Reverse_Alignment_Length", "Alignment", "Forward_Alignment_Positions",
                     "Reverse_Alignment_Positions", "Forward_Wide_Position","Reverse_Wide_Position",
                     "Target_Forward_Found", "Target_Reverse_Found", "Deletions_Forward", "Deletions_Reverse",
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
      # Prepare the variables as a string. This is necessary for the input of the aligning function.
      forwardString   <- toString(IDuniqueTable[r, "Forward"])
      reverseString   <- c2s(rev(comp(s2c(toString(IDuniqueTable[r, "Reverse"])))))
      #------------------------------------------------------
      # Align the candidates with the genome
      #------------------------------------------------------
      # The return of gotoh is a string divided in five parts with the separator "++++"
      # -- [1] is the verbose alignment result
      # -- [2] is a string with the insertions and deletions and missmatches
      # -- [3] is a string with the insertions and deletions and missmatches,
      #    but from the subject (amplicon) coordinates
      # -- [4] is the alignment of the pattern.
      # -- [5] is the alignment of the subject.

      # If we are using both FASTQ files, or only the forward, align the
      # forward sequence with the amplicon.
      if(fastqfiles == 0 || fastqfiles == 1){
        alignForward <- gRCPP(forwardString,
                              amplicon,
                              scoring_matrix,
                              gap_opening,
                              gap_extension,
                              gap_ending,
                              far_indels)
      }
      # If we are using both FASTQ files, or only the reverse, align the
      # reverse sequence with the amplicon
      if(fastqfiles == 0 || fastqfiles == 2){
        alignReverse <- gRCPP(reverseString,
                              amplicon,
                              scoring_matrix,
                              gap_opening,
                              gap_extension,
                              gap_ending,
                              far_indels)
      }
      # Extract the information from the alignments (each part to a different variable)
      #------------------------------------------------------
      partForward <- gsub("\n","", unlist(strsplit(alignForward, "++++", fixed = TRUE)))
      partReverse <- gsub("\n","", unlist(strsplit(alignForward, "++++", fixed = TRUE)))
      # Extract the verbose summaries and the lite indel representation
      forwardSummaryString <- partForward[1]
      forwardLiteString <- partForward[2]
      reverseSummaryString <- partReverse[1]
      reverseLiteString <- partReverse[2]
      # Extract the lite indel representation that is set respect the genome coordinates
      forwardGenomeRelativeLiteString <- partForward[3]
      reverseGenomeRelativeLiteString <- partReverse[3]
      # Extract each part of the alignment forward/reverse in different strings
      forwardPatternAlignment <- partForward[4]
      forwardSubjectAlignment <- partForward[5]
      reversePatternAlignment <- partReverse[4]
      reverseSubjectAlignment <- partReverse[5]

      # We need to count the lenght of the alignment for the forward and the
      # reverse; but only if they are actually relevant. Otherwise they are
      # of lenght 0.
      forwardAlignmentLength <- ifelse(fastqfiles == 0 || fastqfiles == 1, nchar(forwardPatternAlignment), 0)
      reverseAlignmentLength <- ifelse(fastqfiles == 0 || fastqfiles == 2, nchar(reversePatternAlignment), 0)
      #------------------------------------------------------
      #Get the forward and reverse events coordinates
      #------------------------------------------------------
      forwardData <- getEventInfo(forwardLiteString)
      reverseData <- getEventInfo(reverseLiteString)
      forwardInsertionsTotal  <- forwardData[[1]]
      forwardInsertionsStarts <- forwardData[[2]]
      forwardInsertionsEnds   <- forwardData[[3]]
      forwardDeletionsTotal   <- forwardData[[4]]
      forwardDeletionsStarts  <- forwardData[[5]]
      forwardDeletionsEnds    <- forwardData[[6]]
      reverseInsertionsTotal  <- reverseData[[1]]
      reverseInsertionsStarts <- reverseData[[2]]
      reverseInsertionsEnds   <- reverseData[[3]]
      reverseDeletionsTotal   <- reverseData[[4]]
      reverseDeletionsStarts  <- reverseData[[5]]
      reverseDeletionsEnds    <- reverseData[[6]]

      #------------------------------------------------------
      # Extract the information from the lite string in genome coordinates
      #------------------------------------------------------
      forwardGenomeData <- getEventInfo(forwardGenomeRelativeLiteString)
      reverseGenomeData <- getEventInfo(reverseGenomeRelativeLiteString)
      forwardGenomeInsertionsTotal  <- forwardGenomeData[[1]]
      forwardGenomeInsertionsStarts <- forwardGenomeData[[2]]
      forwardGenomeInsertionsEnds   <- forwardGenomeData[[3]]
      forwardGenomeDeletionsTotal   <- forwardGenomeData[[4]]
      forwardGenomeDeletionsStarts  <- forwardGenomeData[[5]]
      forwardGenomeDeletionsEnds    <- forwardGenomeData[[6]]
      reverseGenomeInsertionsTotal  <- reverseGenomeData[[1]]
      reverseGenomeInsertionsStarts <- reverseGenomeData[[2]]
      reverseGenomeInsertionsEnds   <- reverseGenomeData[[3]]
      reverseGenomeDeletionsTotal   <- reverseGenomeData[[4]]
      reverseGenomeDeletionsStarts  <- reverseGenomeData[[5]]
      reverseGenomeDeletionsEnds    <- reverseGenomeData[[6]]

      forwardDeletions <- forwardDeletionsTotal
      reverseDeletions <- reverseDeletionsTotal

      Forward_Deletions_String <- forwardLiteString
      Reverse_Deletions_String <- reverseLiteString
      Forward_Deletions_Genome_Coordinates_String <- forwardGenomeRelativeLiteString
      Reverse_Deletions_Genome_Coordinates_String <- reverseGenomeRelativeLiteString

      # Lets find out the offsets of the alignment position with respect the alignment.
      # Get the position of the alignment and how wide it is
      forwardAllPositions       <- upperGroups(toString(forwardSubjectAlignment))
      forwardAlignmentPositions <- start(forwardAllPositions)
      forwardWidePosition       <- width(forwardWidePosition)

      reverseAllPositions       <- upperGroups(toString(reverseSubjectAlignment))
      reverseAlignmentPositions <- start(reverseAllPositions)
      reverseWidePosition       <- width(reverseWidePosition)
      #--------------------------------------------------
      # Write the alignments into disk
      #--------------------------------------------------
      # If we have an alignment, write it into a TXT file
      # The name of the alignment file is keep regarless of the option of writing on disk or not.
      # Is a nice visual aid to locate uniques with both primers without looking at the TRUE/FALSE of the primers.
      Alignment_File <- paste(uniqueIndex,".txt",sep='')

      if (write_alignments >= 1) {
        # Write the summary master alignment
        writeLines(c(paste(uniqueIndex, toString(uniqueTable$Total[uniqueIndex])),
                     forwardPatternAlignment,
                     reversePatternAlignment), masterFileConn)
      }
      if (write_alignments >= 2) {
        # Write the verbose alignment into an individual file
        alignmentName <- paste(uniqueIndex,".txt",sep='')
        # Write the uber alignment file
        writeLines(c(alignmentName, "\n FORWARD AND AMPLICON: \n",
                     forwardSummaryString, "\n REVERSE AND AMPLICON: \n",
                     reverseSummaryString, "\n"), uberAlignmentFD)
      }
      if (write_alignments == 3) {
        Alignment_File_Path <- paste(currentAlignmentsResultsFolderName, "/", alignmentName, sep='')
        fileConn<-file(Alignment_File_Path)
        writeLines(c(forwardSummaryString,"\n + \n",reverseSummaryString), fileConn)
        close(fileConn)
      }
      #--------------------------------------------------
      # Fill the deletion table DEPRECRATED
      #--------------------------------------------------
      if(forwardDeletionsTotal > 0){
        for(m in 1:forwardDeletionsTotal){
          # TODO:
          # VERY VERY IMPORTANT
          # The Gotoh has a bug where it returns an interval of [start,end-1] instead of [start,ends] if ends is the
          # end of the amplicon. This piece of code correct that BUT THAT MUST BE CORRECTED IN GOTOH.CPP!!!!
          addjustedEndForward <- forwardDeletionsEnds[m]
          if(forwardDeletionsEnds[m] == nchar(amplicon)-1){
            addjustedEndForward <- addjustedEndForward + 1
          }
          addjustedRelativeEndForward <- forwardGenomeDeletionsEnds[m]
          if(forwardGenomeDeletionsEnds[m] == nchar(amplicon)-1){
            addjustedRelativeEndForward <- addjustedRelativeEndForward + 1
          }
        }
      }

      if(reverseDeletionsTotal > 0){
        for(m in 1:reverseDeletionsTotal){
          # TODO:
          # VERY VERY IMPORTANT
          # The Gotoh has a bug where it returns an interval of [start,end-1] instead of [start,ends] if ends is the
          # end of the amplicon. This piece of code correct that BUT THAT MUST BE CORRECTED IN GOTOH.CPP!!!!
          addjustedEndReverse <- reverseDeletionsEnds[m]
          if(reverseDeletionsEnds[m] == nchar(amplicon)-1){
            addjustedEndReverse <- addjustedEndReverse + 1
          }
          addjustedRelativeEndReverse <- reverseGenomeDeletionsEnds[m]
          if(reverseGenomeDeletionsEnds[m] == nchar(amplicon)-1){
            addjustedRelativeEndReverse <- addjustedRelativeEndReverse + 1
          }
        }
      }
      # Write the sequence into disk
      # -- Each folder has a comprehensive list of all the sequences with some stats
      # -- This part add the sequence to its file.
      #--------------------------------------------------
      writeLines( paste(as.numeric(uniqueIndex),
                        as.numeric(uniqueTable$Total[uniqueIndex]),
                        as.numeric(uniqueTable$Frequency[uniqueIndex]),
                        uniqueTable$Forward[uniqueIndex],
                        uniqueTable$Reverse[uniqueIndex],
                        gsub("\n", "", Forward_Deletions_String),
                        gsub("\n", "", Reverse_Deletions_String),
                        gsub("\n", "", Forward_Deletions_Genome_Coordinates_String),
                        gsub("\n", "", Reverse_Deletions_Genome_Coordinates_String),
                        as.numeric(forwardAlignmentLength),
                        as.numeric(reverseAlignmentLength),
                        Alignment_File,
                        forwardAlignmentPositions,
                        reverseAlignmentPositions,
                        forwardWidePosition,
                        reverseWidePosition,
                        targetFoundForward,
                        targetFoundReverse,
                        as.numeric(forwardDeletionsTotal),
                        as.numeric(reverseDeletionsTotal),
                        sep = '\t') , uniqueTableFD)
    }

    # Add the currect statistics to the config table
    configTable$Sum_Target_Forward_Found[i] <- sum(uniqueTable$forwardFound)
    configTable$Sum_Target_Reverse_Found[i] <- sum(uniqueTable$reverseFound)
    ###???
    configTable$Sum_Deletions_Forward[i] <- configTable$Sum_Target_Reverse_Found[i] + (forwardDeletions   * uniqueTable$Total[i])
    configTable$Sum_Deletions_Reverse[i] <- sumDeletionsReverse
    configTable$Sum_Is_Sequence[i] <- configTable$Sum_Is_Sequence[i]  + uniqueTable$Total[i]
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

  write.table(configTable, paste(alignmentFolder, "/", barcode, "_configFile_results", sep = '') , sep="\t")
  write.table(barcodeTable, paste(alignmentFolder, "/", barcode, "_barcodeFile_results", sep = ''), sep="\t")
  return()
}
