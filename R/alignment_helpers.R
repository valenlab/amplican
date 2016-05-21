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

  logFileConn    <- file(paste0(resultsFolder, "/", barcode, "_SUBLOG.txt"), open = "at")
  tempFileConn <- file(paste0(resultsFolder, "/", barcode, "_TEMPORALFILES.txt"), open = "at")
  unnasignedFileName <- paste(resultsFolder, "/", barcode, "_unassigned.txt", sep = '', collapse = '')

  #Read Reads for this Barcode
  forwardsTable <- readFastq(configTable$Forward_Reads_File[1])
  reversesTable <- readFastq(configTable$Reverse_Reads_File[1])
  #Filter Reads
  goodq <- goodBaseQuality(forwardsTable, min = min_quality) & goodBaseQuality(reversesTable, min = min_quality)
  avrq <- goodAvgQuality(forwardsTable, avg = average_quality) & goodAvgQuality(reversesTable, avg = average_quality)
  nucq <- alphabetQuality(forwardsTable) & alphabetQuality(reversesTable)
  goodReads <- goodq & avrq & nucq

  barcodeTable <- data.frame(Total_Pre_N_Filter = length(goodReads), Total_Post_N_Filter = sum(goodReads))

  forwardsTable <- forwardsTable[goodReads]
  reversesTable <- reversesTable[goodReads]
  #Unique reads
  uniqueTable <- data.frame(as.character(forwardsTable@sread), as.character(reversesTable@sread))
  colnames(uniqueTable) <- c("Forward","Reverse")
  uniqueTable$Total <- paste0(uniqueTable$Forward, uniqueTable$Reverse)
  uniqueTable <- aggregate(Total ~  Reverse + Forward, uniqueTable, length)
  uniqueTable$Frequency <- uniqueTable$Total / sum(uniqueTable$Total)
  uniqueTable <- uniqueTable[order(uniqueTable$Forward, uniqueTable$Reverse),]

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
    reversePrimerRC  <- c2s(rev(comp(s2c(reversePrimer))))
    amplicon <- toString(configTable[i, "Amplicon"])
    forwardPrimerLength <- nchar(forwardPrimer)
    reversePrimerLength <- nchar(reversePrimer)
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

  }

      cutSites <- upperGroups(amplicon)
      alignmentPositions <- start(cutSites)
      widePosition <- width(cutSites)

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
targetFound  <- checkTarget(targetPrimer, amplicon, currentID, barcode, logFileConn)
        # If the primers are not in the amplicon tell it to the user.
primersFound <- checkPrimers(forwardPrimer, reversePrimerRC, amplicon, currentID, barcode, logFileConn)
        # If the amplicon has no CUT sites specified.
apFound <- checkPositions(start(cutSites), amplicon, currentID, barcode, logFileConn)


    #------------------------------------------------------
    # PREPARE THE UNASSIGNED VARIABLES
    #------------------------------------------------------
      # Also, we are going to keep a list for the unique sequences that are not in use
      # We are going to keep TRUE those which are unnasigned, so we don't need to complement this array later
      unnasignedSequences <- rep(TRUE, totalUniqueSequences)
      forwardFound        <- rep(FALSE, totalUniqueSequences)
      reverseFound        <- rep(FALSE, totalUniqueSequences)

      for (uniqueIndex in dim(uniqueTable)[1]) {

      }
        # Search for the forward, reverse and targets
          # Get the sequences
        candidateForwardSequence <- toString(uniqueTable$Forward[uniqueIndex])
        candidateReverseSequence <- toString(uniqueTable$Reverse[uniqueIndex])

        forwardFound <- grepl(forwardPrimer, candidateForwardSequence, ignore.case=TRUE)
        reverseFound <- grepl(reversePrimer, candidateReverseSequence, ignore.case=TRUE)
        targetFoundForward <- grepl(targetPrimer, candidateForwardSequence, ignore.case=TRUE)
        targetFoundReverse <- grepl(c2s(rev(comp(s2c(targetPrimer)))), candidateReverseSequence, ignore.case=TRUE)

        # Addjust the finding of primers. If we are not going to use either the
        # forward sequences, or the reverse sequences, finding or not that primer
        # is irrelevant.
        if (fastqfiles == 1) {
          reverseFound <- TRUE
        }
        if (fastqfiles == 2) {
            forwardFound <- TRUE
        }

        # If the forward and reverse are found, do the alignments and everything else
        if (forwardFound && reverseFound) {
          # Mark this one as assigned
          unnasignedSequences[uniqueIndex] <- FALSE

          # Find the reverse complementary of the reverse sequence for the multiple alignment
          candidateComplementarySequence <- c2s(rev(comp(s2c(candidateReverseSequence))))

          # Prepare the variables as a string. This is necessary for the input of the aligning function.
          forwardString   <- candidateForwardSequence
          reverseString   <- toString(candidateComplementarySequence)[1]
          ampliconString  <- amplicon

          # Here is where the alignments are going to be stored
          alignForward <- ""
          alignReverse <- ""

          #------------------------------------------------------
          # Align the candidates with the genome
          #------------------------------------------------------
            # The return of gotoh is a string divided in five parts with the separator "++++"
            # -- The first part is the verbose alignment result
            # -- The second part is a string with the insertions and deletions and missmatches
            # -- The third part is a string with the insertions and deletions and missmatches but from the subject (amplicon) coordinates
            # -- The fourth part is the alignment of the pattern.
            # -- The fifth part is the alignment of the subject.

            # If we are using both FASTQ files, or only the forward, align the
            # forward sequence with the amplicon.
            if(fastqfiles == 0 || fastqfiles == 1){
              alignForward <- gRCPP(forwardString,
                                    ampliconString,
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
                                    ampliconString,
                                    scoring_matrix,
                                    gap_opening,
                                    gap_extension,
                                    gap_ending,
                                    far_indels)
            }
          # Extract the information from the alignments (each part to a different variable)
          #------------------------------------------------------
            # Extract the verbose summaries and the lite indel representation
            forwardSummaryString <- unlist(strsplit(alignForward, "++++", fixed = TRUE))[1]
            forwardLiteString <- unlist(strsplit(alignForward, "++++", fixed = TRUE))[2]
            reverseSummaryString <- unlist(strsplit(alignReverse, "++++", fixed = TRUE))[1]
            reverseLiteString <- unlist(strsplit(alignReverse, "++++", fixed = TRUE))[2]

            # Extract the lite indel representation that is set respect the genome coordinates
            forwardGenomeRelativeLiteString <- unlist(strsplit(alignForward, "++++", fixed = TRUE))[3]
            reverseGenomeRelativeLiteString <- unlist(strsplit(alignReverse, "++++", fixed = TRUE))[3]

            # Extract each part of the alignment forward/reverse in different strings
            forwardPatternAlignment <- gsub("\n","",unlist(strsplit(alignForward, "++++", fixed = TRUE))[4])
            forwardSubjectAlignment <- gsub("\n","",unlist(strsplit(alignForward, "++++", fixed = TRUE))[5])
            reversePatternAlignment <- gsub("\n","",unlist(strsplit(alignReverse, "++++", fixed = TRUE))[4])
            reverseSubjectAlignment <- gsub("\n","",unlist(strsplit(alignReverse, "++++", fixed = TRUE))[5])

            # We need to count the lenght of the alignment for the forward and the
            # reverse; but only if they are actually relevant. Otherwise they are
            # of lenght 0.
            if (fastqfiles == 0 || fastqfiles == 1) {
              forwardAlignmentLength <- nchar(forwardPatternAlignment)
            } else {
              forwardAlignmentLength <- 0
            }

            if (fastqfiles == 0 || fastqfiles == 2) {
              reverseAlignmentLength <- nchar(reversePatternAlignment)
            } else {
              reverseAlignmentLength <- 0
            }
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
          # TODO: Adjust for more than one
          forwardAllPositions       <- upperGroups(toString(forwardSubjectAlignment))
          forwardAlignmentPositions <- forwardAllPositions[[2]][1]
          forwardWidePosition       <- forwardAllPositions[[3]][1]

          reverseAllPositions       <- upperGroups(toString(reverseSubjectAlignment))
          reverseAlignmentPositions <- reverseAllPositions[[2]][1]
          reverseWidePosition       <- reverseAllPositions[[3]][1]
          #--------------------------------------------------
          # Write the alignments into disk
          #--------------------------------------------------
            # If we have an alignment, write it into a TXT file
            # The name of the alignment file is keep regarless of the option of writing on disk or not.
            # Is a nice visual aid to locate uniques with both primers without looking at the TRUE/FALSE of the primers.
            Alignment_File <- paste(uniqueIndex,".txt",sep='')

            if (write_alignments > 0) {
              if (write_alignments >= 1) {
                # Write the summary master alignment
                writeLines(c(paste(uniqueIndex,toString(uniqueTable$Total[uniqueIndex])),forwardPatternAlignment,reversePatternAlignment), masterFileConn)
              }
              if (write_alignments >= 2) {
                # Write the verbose alignment into an individual file
                alignmentName <- paste(uniqueIndex,".txt",sep='')
                # Write the uber alignment file
                writeLines(c(alignmentName,       "\n FORWARD AND AMPLICON: \n",
                             forwardSummaryString,"\n REVERSE AND AMPLICON: \n",
                             reverseSummaryString,"\n"), uberAlignmentFD)
              }
              if (write_alignments == 3) {
                Alignment_File_Path <- paste(currentAlignmentsResultsFolderName, "/", alignmentName, sep='')
                fileConn<-file(Alignment_File_Path)
                writeLines(c(forwardSummaryString,"\n + \n",reverseSummaryString), fileConn)
                close(fileConn)
              }
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
                if(forwardDeletionsEnds[m] == nchar(ampliconString)-1){
                  addjustedEndForward <- addjustedEndForward + 1
                }
                addjustedRelativeEndForward <- forwardGenomeDeletionsEnds[m]
                if(forwardGenomeDeletionsEnds[m] == nchar(ampliconString)-1){
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
                if(reverseDeletionsEnds[m] == nchar(ampliconString)-1){
                  addjustedEndReverse <- addjustedEndReverse + 1
                }
                addjustedRelativeEndReverse <- reverseGenomeDeletionsEnds[m]
                if(reverseGenomeDeletionsEnds[m] == nchar(ampliconString)-1){
                  addjustedRelativeEndReverse <- addjustedRelativeEndReverse + 1
                }
              }
            }
          # Write the sequence into disk
          # -- Each folder has a comprehensive list of all the sequences with some stats
          # -- This part add the sequence to its file.
          #--------------------------------------------------
            Forward_Deletions_String                    <- gsub("\n", "", Forward_Deletions_String)
            Reverse_Deletions_String                    <- gsub("\n", "", Reverse_Deletions_String)
            Forward_Deletions_Genome_Coordinates_String <- gsub("\n", "", Forward_Deletions_Genome_Coordinates_String)
            Reverse_Deletions_Genome_Coordinates_String <- gsub("\n", "", Reverse_Deletions_Genome_Coordinates_String)

            writeLines( paste(as.numeric(uniqueIndex),
                              as.numeric(uniqueTable$Total[uniqueIndex]),
                              as.numeric(uniqueTable$Frequency[uniqueIndex]),
                              uniqueTable$Forward[uniqueIndex],
                              uniqueTable$Reverse[uniqueIndex],
                              Forward_Deletions_String,
                              Reverse_Deletions_String,
                              Forward_Deletions_Genome_Coordinates_String,
                              Reverse_Deletions_Genome_Coordinates_String,
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

          #--------------------------------------------------
          # Update the config statistics
          #--------------------------------------------------
            sumTargetForwardFound     <- sumTargetForwardFound + (targetFoundForward * uniqueTable$Total[uniqueIndex])
            sumTargetReverseFound     <- sumTargetReverseFound + (targetFoundReverse * uniqueTable$Total[uniqueIndex])
            sumDeletionsForward       <- sumDeletionsForward   + (forwardDeletions   * uniqueTable$Total[uniqueIndex])
            sumDeletionsReverse       <- sumDeletionsReverse   + (reverseDeletions   * uniqueTable$Total[uniqueIndex])
            sumIsSequence             <- sumIsSequence         + uniqueTable$Total[uniqueIndex]
          # Go forward one sequence
          uniqueIndex <- uniqueIndex + 1
        } else { # If either forward or reverse is not found
          # Make sure that we are still under the primer scope. Otherwise take the next primer
          # This part is quite complicated to understand; so just trust me.
          forwardPrimerLength <- nchar(forwardPrimer)
          reversePrimerLength <- nchar(reversePrimer)
          forwardSequenceSample <- toupper(substr(candidateForwardSequence, 1, forwardPrimerLength))
          reverseSequenceSample <- toupper(substr(candidateReverseSequence, 1, reversePrimerLength))
          forwardSequenceLength <- nchar(candidateForwardSequence)
          reverseSequenceLength <- nchar(candidateReverseSequence)
          forwardPrimerUpper    <- paste(forwardPrimer,
                                         paste(rep("T", forwardSequenceLength - forwardPrimerLength),
                                               sep='', collapse=''), sep='', collapse='')
          forwardPrimerLower    <- paste(forwardPrimer,
                                         paste(rep("A", forwardSequenceLength - forwardPrimerLength),
                                               sep='', collapse=''), sep='', collapse='')
          reversePrimerUpper    <- paste(reversePrimer,
                                         paste(rep("T", reverseSequenceLength - reversePrimerLength),
                                               sep='', collapse=''), sep='', collapse='')
          reversePrimerLower    <- paste(reversePrimer,
                                         paste(rep("A", reverseSequenceLength - reversePrimerLength),
                                               sep='', collapse=''), sep='', collapse='')
          # If the forwardSequence is bigger than the forward primer upper
          # If the forwardUpper is equal than the forwardSequence AND if the reverseSequence is bigger than reverse primer upper
          # If this was the last sequence in the unique table
          insideScope <- TRUE

          if ((toupper(forwardPrimerUpper) == candidateForwardSequence &&
               toupper(reversePrimerUpper) < candidateReverseSequence) ||
            (toupper(forwardPrimerUpper) <  candidateForwardSequence) || uniqueIndex == nrow(uniqueTable)) {
            insideScope <- FALSE
          }

          # If we are inside the scope of the primers
          if (insideScope == TRUE) {
            # Keep track on which primer was not found
            forwardFound[uniqueIndex] <- forwardFound
            reverseFound[uniqueIndex] <- reverseFound

            # Pick the next one
            uniqueIndex <- uniqueIndex + 1
          } else { # If we passed the scope of the current primers, take the next two primers

            # Add the currect statistics to the config table
            configTable$Sum_Target_Forward_Found[configIndex] <- sumTargetForwardFound
            configTable$Sum_Target_Reverse_Found[configIndex] <- sumTargetReverseFound
            configTable$Sum_Deletions_Forward[configIndex]    <- sumDeletionsForward
            configTable$Sum_Deletions_Reverse[configIndex]    <- sumDeletionsReverse
            configTable$Sum_Is_Sequence[configIndex]          <- sumIsSequence

            # Barcode stuff
            configTable$Total_Pre_N_Filter[configIndex]          <- totalSequencesPrefilter
            configTable$Total_Post_N_Filter[configIndex]         <- totalSequencesPostfilter
            configTable$Experiment_Unique_Sequences[configIndex] <- totalUniqueSequences
            configTable$Total_Experiment_Sequences[configIndex]  <- totalSequences

            # Warnings
            configTable$Found_Target[configIndex]  <- targetFound
            configTable$Found_Primers[configIndex] <- primersFound
            configTable$Found_AP[configIndex]      <- apFound

            # Take the next one
            configSubsetIndex <- configSubsetIndex + 1
            configIndex       <- configIndex       + 1

            # But be carefull!! ; if the next primers are the same as the previous one, we need to revert
            # the uniqueIndex to the previous state. If we still have a next primer, check that
            if(configSubsetIndex <= subConfigLength){
              # Save the old primers
              lastForwardPrimer <- forwardPrimer
              lastRevesrePrimer <- reversePrimer

              # Get the new primers
              forwardPrimer <- toString(configSubset$Forward_Primer[configSubsetIndex])
              reversePrimer <- toString(configSubset$Reverse_Primer[configSubsetIndex])

              # If they are the same, revert the unique Index
              if(lastForwardPrimer  == forwardPrimer && lastRevesrePrimer  == reversePrimer){
                print("Repeat primers, they are the same")
                uniqueIndex <- equalPrimersIndex
              } else { # Otherwise, update the primer index to the current unique
                equalPrimersIndex <- uniqueIndex
              }
            }

            # You also need to restart the variables for the next row, if we have next row.
            # If you look inside here, you will notice that we close the FDs from the previous one.
            # Don't worry for the last one, because they are closed outside the loop, just like the
            # first init is done outside the loop too.
            if(configSubsetIndex <= subConfigLength){
              #------------------------------------------------------
              # INITIALIZE THE CONFIG ROW VARIABLES
              # The same as before, we in here we are also closing the FD if needed
              #------------------------------------------------------
                # Progress feedback for the user
                print(paste("The processor",barcode,"is doing line",configIndex,"of",nrow(configTable),
                            "Progress of: ",round(configIndex/nrow(configTable)*100,2),"%"))
                #------------------------
                # CHECK FOR BASIC INFO
                #------------------------
                  # Variables
                  totalPreNFilter           <- 0
                  totalPostNFilter          <- 0
                  experimentUniqueSequences <- 0
                  totalExperimentSequences  <- 0
                  sumForwardFound           <- 0
                  sumReverseFound           <- 0
                  sumTargetForwardFound     <- 0
                  sumTargetReverseFound     <- 0
                  sumDeletionsForward       <- 0
                  sumDeletionsReverse       <- 0
                  sumIsSequence             <- 0

                  # Specials
                  currentID                          <- configTable[configIndex,"ID"]
                  barcode                     <- configTable[configIndex,"Barcode"] # Should be equal to the variable: barcode

                  forwardPrimer                      <- toString(configTable[configIndex,"Forward_Primer"])
                  reversePrimer                      <- toString(configTable[configIndex,"Reverse_Primer"])
                  targetPrimer                       <- toString(configTable[configIndex,"Target_Primer"])
                  reversePrimerRC  <- reverseComplement(reversePrimer)
                  amplicon                           <- toString(configTable[configIndex,"Genome"])
                  forwardPrimerLength                <- nchar(forwardPrimer)
                  reversePrimerLength                <- nchar(reversePrimer)

                  currentIDFolderName     <- paste(alignmentFolder,"/",currentID,"_",barcode,sep = '')
                  currentAlignmentsResultsFolderName <- paste(currentIDFolderName,"/alignments",sep='')
                  masterAlignmentFilePath            <- paste(currentIDFolderName,"/alignments.txt",sep='')
                  uberAlignmentFilePath              <- paste(currentIDFolderName,"/verbose.txt",sep='')
                  deletionFilePath                   <- paste(currentIDFolderName,"/",currentID,"_",
                                                              barcode,"_deletions.txt",sep = '')
                  deletionRelativeFilePath           <- paste(currentIDFolderName,"/",currentID,"_",
                                                              barcode,"_deletionsRelative.txt",sep = '')
                  uniqueTablePath                    <- paste(currentIDFolderName,"/",currentID,"_",
                                                              barcode,"_uniques.txt",sep = '')

                  targetFound  <- FALSE
                  primersFound <- FALSE
                  apFound      <- FALSE
                  #------------------------
                  # CLOSE PREVIOUS FILE DESCRIPTORS
                  #------------------------
                    # Close all the file descriptors that we had open so far at barcode level
                    # These three are open for sure
                    #                   close(deletionFD)
                    #                   close(deletionRelativeFD)
                    close(uniqueTableFD)
                    # These two might not be needed, check it out and close if necessary
                    if (isOpen(masterFileConn)) {
                      close(masterFileConn)
                    }
                    if (isOpen(uberAlignmentFD)) {
                      close(uberAlignmentFD)
                    }

                # Create folders
                dir.create(file.path(currentIDFolderName), showWarnings = TRUE)
                dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)
                # Update the target Primer if necessary
                reverseAmplicon <- configTable[configIndex, "Strand"]
                #               print(reverseAmplicon)
                if(reverseAmplicon == 1){
                  targetPrimer  <- reverseComplement(targetPrimer)
                }
                # Get the position of the alignment and how wide it is
                # TODO: Adjust for more than one
                allPositions       <- upperGroups(amplicon)
                alignmentPositions <- allPositions[[2]][1]
                widePosition       <- allPositions[[3]][1]
                #------------------------
                # PREPARE THE FILE DESCRIPTORS
                #------------------------
                  # If we need to write the alignment, prepare the proper file descriptors
                  if (write_alignments > 0) {
                    if (write_alignments >= 1) {
                      # Prepare the summary alignment file
                      masterFileConn     <- file(masterAlignmentFilePath, open="at")
                      writeLines(c(toString(amplicon),"\n"), masterFileConn)
                    }
                    if (write_alignments >= 2) {
                      # Start the file where the alignments are accumulated
                      uberAlignmentFD     <- file(uberAlignmentFilePath, open="at")
                    }
                  }
                  #                 deletionFD         <-  file(deletionFilePath, open="at")
                  #                 deletionRelativeFD <-  file(deletionRelativeFilePath, open="at")
                  uniqueTableFD      <- file(uniqueTablePath, open="at")
                  #                 writeLines( paste("Gene_ID", "Hash_ID", "Start","End","Freq","Type","Location", "Wide",sep = '\t')
                  #                             , deletionFD)
                  #                 writeLines( paste("Gene_ID", "Hash_ID", "Start","End","Freq","Type","Location", "Wide",sep = '\t')
                  #                             , deletionRelativeFD)
                  writeLines( paste("Total", "Frequency", "Forward", "Reverse", "Forward_Deletions_String",
                                    "Reverse_Deletions_String", "Forward_Deletions_Genome_Coordinates_String",
                                    "Reverse_Deletions_Genome_Coordinates_String",
                                    "Forward_Alignment_Length", "Reverse_Alignment_Length",  "Alignment",
                                    "Forward_Alignment_Positions", "Reverse_Alignment_Positions", "Forward_Wide_Position",
                                    "Reverse_Wide_Position", "Target_Forward_Found",
                                    "Target_Reverse_Found", "Deletions_Forward", "Deletions_Reverse",
                                    sep = '\t'), uniqueTableFD)
                #------------------------
                # CHECK FOR WARNINGS
                #------------------------
                  # If the target is not in the amplicon, tell it to the user.
                  targetFound  <- checkTarget(targetPrimer, amplicon, currentID, barcode, configFilePath,
                                              logFileConn)
                  # If the primers are not in the amplicon tell it to the user.
                  primersFound <- checkPrimers(forwardPrimer, reversePrimerRC, amplicon,
                                               currentID, barcode, configFilePath, logFileConn)
                  # If the amplicon has no alignment positions, tell the user.
                  apFound      <- checkPositions(alignmentPositions, amplicon, currentID, barcode,
                                                 configFilePath, logFileConn)
              }
            # Reboot the statistics
            totalPreNFilter           <- 0
            totalPostNFilter          <- 0
            experimentUniqueSequences <- 0
            totalExperimentSequences  <- 0
            sumForwardFound           <- 0
            sumReverseFound           <- 0
            sumTargetForwardFound     <- 0
            sumTargetReverseFound     <- 0
            sumDeletionsForward       <- 0
            sumDeletionsReverse       <- 0
            sumIsSequence             <- 0
          }
        }
      } # While loop
    } # If unique is not empty


    uniqueTable["Forward_Found"] <- forwardFound
    uniqueTable["Reverse_Found"] <- reverseFound
    uniqueTable <- uniqueTable[unnasignedSequences, ]
    write.table(uniqueTable, file = unnasignedFileName, quote = FALSE, sep = "\t")

    # Close all the file descriptors that we had open so far at barcode level
    close(uniqueTableFD)

    # These two might not be needed, check it out and close if necessary
    if (masterFileConnOpen) {
      close(masterFileConn)
    }
    if (uberAlignmentFDOpen) {
      close(uberAlignmentFD)
    }

  close(logFileConn)
  close(tempFileConn)

  write.table(configTable, paste(alignmentFolder, "/", barcode, "_configFile_results", sep = '') , sep="\t")
  write.table(barcodeTable, paste(alignmentFolder, "/", barcode, "_barcodeFile_results", sep = ''), sep="\t")
  return()
}
