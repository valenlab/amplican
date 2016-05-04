#' This function return the reverse complementary dna string of a given dna
#' string.
#
#' The function respect uppercases. Any other character that is not
#' a,c,g,t,A,C,G,T is ignored.
#'
#' @param dna (String) A string with the nucleotide sequence.
#' @return (String) A string with the nucleotides reversed and complemented
reverseComplement <- function(dna){
  myResult <- dna
  # SWAP a <-> t
  myResult <- gsub("a", "1", myResult)
  myResult <- gsub("t", "a", myResult)
  myResult <- gsub("1", "t", myResult)
  # SWAP c <-> g
  myResult <- gsub("c", "1", myResult)
  myResult <- gsub("g", "c", myResult)
  myResult <- gsub("1", "g", myResult)
  # SWAP A <-> T
  myResult <- gsub("A", "1", myResult)
  myResult <- gsub("T", "A", myResult)
  myResult <- gsub("1", "T", myResult)
  # SWAP C <-> G
  myResult <- gsub("C", "1", myResult)
  myResult <- gsub("G", "C", myResult)
  myResult <- gsub("1", "G", myResult)
  myResult <- paste( rev(strsplit(myResult, NULL)[[1]]) , collapse='')
  return (myResult)
}

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
    if(alignmentInsertionsTotal>0){
      for(i in 1:alignmentInsertionsTotal){
        indelPair <- unlist(strsplit(alignmentInsertionsBoundariesArray[i],",",fixed = TRUE))
        alignmentInsertionsStarts[i] <- as.numeric(indelPair[1])
        alignmentInsertionsEnds[i]   <- as.numeric(indelPair[2])
      }}
    else{
      alignmentInsertionsStarts <- NULL
      alignmentInsertionsEnds   <- NULL
    }

    # Fill the arrays for the event deletions
    alignmentDeletionsBoundariesData <- unlist(strsplit(alignmentDeletionsData,"@",
                                                        fixed = TRUE))[3]
    alignmentDeletionsBoundariesArray <- unlist(strsplit(alignmentDeletionsBoundariesData,"*",
                                                         fixed = TRUE))
    if(alignmentDeletionsTotal>0){
      for(l in 1:alignmentDeletionsTotal){
        indelPair <- unlist(strsplit(alignmentDeletionsBoundariesArray[l],",",fixed = TRUE))
        alignmentDeletionsStarts[l] <- as.numeric(indelPair[1])
        alignmentDeletionsEnds[l]   <- as.numeric(indelPair[2])
      }}
    else{
      alignmentDeletionsStarts <- NULL
      alignmentDeletionsEnds   <- NULL
    }

    # Fill the arrays for the missmatches
    alignmentMissmatchBoundariesData <- unlist(strsplit(alignmentMissmatchData,"@",
                                                        fixed = TRUE))[3]
    alignmentMissmatchBoundariesArray <- unlist(strsplit(alignmentMissmatchBoundariesData,"*",
                                                         fixed = TRUE))
    if(alignmentMissmatchesTotal>0){
      for(l in 1:alignmentMissmatchesTotal){

        missMatchTrio <- unlist(strsplit(alignmentMissmatchBoundariesArray[l],",",fixed = TRUE))
        alignmentMutationsOriginal[l]  <- missMatchTrio[2]
        alignmentMutationsMutated[l]   <- missMatchTrio[3]
        alignmentMutationsPosition [l] <- as.numeric(missMatchTrio[1])

      }}
    else{
      alignmentMutationsOriginal <- NULL
      alignmentMutationsMutated  <- NULL
      alignmentMutationsPosition <- NULL
    }

    return (list(alignmentInsertionsTotal, alignmentInsertionsStarts, alignmentInsertionsEnds,
                 alignmentDeletionsTotal,  alignmentDeletionsStarts, alignmentDeletionsEnds,
                 alignmentMutationsOriginal, alignmentMutationsMutated, alignmentMutationsPosition))

  } else {
    return (list(0, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL))
  }
}

#' For a given string, detect how many groups of uppercases are there, where are
#' they, and how long are they.
#'
#' For example:
#'
#'    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
#'
#' Has 4 groups of uppercases, and they have length 7,4,1, and 3 respectevily.
#' The function returns a list with the length of each one. So the length of the
#' list represent.
#' The function returns a list with 3 elements:
#' -- How many groups we found. If will be 1 if we found nothing.
#' -- The starting position for each AP. It would be -1 if we found nothing.
#' -- The wide for each AP. It will be 1 if we found nothing.
#' @param candidate: (String) A string with the nucleotide sequence.
#' @return (List) A list with three element. The first element is an integer with the
#'          total of upper case group; lets say is X. The second element is a
#'          list of integer of size X with the start position of each group. The
#'          third element is a list of integer with the wide of each group.
#'
countUppercaseGroups <- function(candidate){
  # Add a tiny lowercase letter at the end of candidate. That will trick gregexpr, and avoid making
  #weird cases later.
  candidate <- paste(candidate,"a",sep='')

  alignmentPositionsStarts <- gregexpr("(?<![A-Z])[A-Z]",  candidate,
                                       ignore.case = FALSE, perl = TRUE, fixed = FALSE)[[1]]
  alignmentPositionsEnds   <- gregexpr("([A-Z])([a-z,-])", candidate,
                                       ignore.case = FALSE, perl = TRUE, fixed = FALSE)[[1]]

  total <- length(alignmentPositionsStarts)

  length <- (alignmentPositionsEnds - alignmentPositionsStarts) + 1

  return(list(total,alignmentPositionsStarts,length))
}


#' make alignments
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
#' @param resultFolder (char)
#' @param alignmentFolder (char)
#' @param procesID (int)
#' @param temp_folder (char)
#' @param fastqfiles (char)
#' @return no clue
makeAlignment <- function(configTable,
                          resultFolder,
                          temp_folder,
                          fastqfiles,
                          skip_bad_nucleotides = TRUE,
                          average_quality = 0,
                          min_quality = 0,
                          write_alignments = TRUE,
                          scoring_matrix = "NUC44",
                          gap_opening = 50,
                          gap_extension = 0,
                          gap_ending = FALSE,
                          far_indels = TRUE,
                          processID = 0){

  # Prepare the unnasigned folder
  currentUnassignedFolderName <- paste0(resultsFolder, "/unassigned_sequences")
  dir.create(file.path(currentUnassignedFolderName), showWarnings = F)

  # In this variable we are going to keep track of the last config row. This is important
  # because it is going to be sorted. You need to be very careful with this function
  configIndex <- 1

  #------------------------------------------------------
  # SET THE FILES DESCRIPTOR
  #------------------------------------------------------
    # Make your own file descriptor for the errors log
    subLogFileName <- paste(resultFolder, "/", processID, "_SUBLOG.txt", sep = '')
    logFileConn    <- file(subLogFileName, open = "at")

    # Make your own file descriptor for the temporal files
    tempFileName <- paste(resultFolder, "/", processID, "_TEMPORALFILES.txt", sep = '')
    tempFileConn <- file(tempFileName, open = "at")

    # These are for the configs later one; we need to keep track of them here
    # (Actually we don't, because these variables stay alive and R doesn't care about variable scope
    #  but I'm going to do it anyway; it should be easy to adapt this into another language later on)
    masterFileConn     <- NULL
    uberAlignmentFD    <- NULL
    uniqueTableFD      <- NULL

    # Depending on the flags, these two FDs might not open; we keep track if they are open or not
    masterFileConnOpen     <- FALSE
    uberAlignmentFDOpen    <- FALSE

  # Sort the config table by barcode first, forward primer second, and reverse primer third
  configTable <- configTable[order(configTable$Barcode, configTable$Forward_Primer, configTable$Reverse_Primer),]

  # Add a bunch of columns to the dataframe of the config file that will serve as a summary
  # Please refer to the /doc folder for a verbose explanation for each one of these
    # Barcode stuff
    configTable$Total_Pre_N_Filter <- 0
    configTable$Total_Post_N_Filter <- 0

    configTable$Experiment_Unique_Sequences <- 0
    configTable$Total_Experiment_Sequences  <- 0

    # Several statistics about deletions, cuts and reads
    configTable$Sum_Target_Forward_Found <- 0
    configTable$Sum_Target_Reverse_Found <- 0

    # Deletions, valid or not
    configTable$Sum_Deletions_Forward <- 0
    configTable$Sum_Deletions_Reverse <- 0

    configTable$Sum_Is_Sequence    <- 0

    # Warnings
    configTable$Found_Target  <- 0
    configTable$Found_Primers <- 0
    configTable$Found_AP      <- 0
  # Get the barcode subdataframe, this is the barcode and FWD and RVS files
  # Simplify the rows so we don't have duplicated.
  barcodeTable   <- unique(configTable[,c("Barcode","Forward_Reads_File", "Reverse_Reads_File")])

  # Add a bunch of columns to the barcode table. Although this info is already on the config, it actually belongs
  # to the barcodes
    # Barcode stuff
    barcodeTable$Total_Pre_N_Filter <- 0
    barcodeTable$Total_Post_N_Filter <- 0

    barcodeTable$Experiment_Unique_Sequences <- 0
    barcodeTable$Total_Experiment_Sequences  <- 0

  #------------------------------------------------------
  # PREPARE SOME TIMING VARIABLES
  #------------------------------------------------------
    # We are going to create another dataframe for timing stuff.
    # Each row of this dataFrame represent each row of the dataframe
    configTiming <- data.frame(matrix(0, nrow = nrow(configTable), ncol = 15))
    colnames(configTiming) <- c("Gene_ID", "Barcode", "Total_Time","Alignment_Time", "Search_Primers",
                                "Write_Sequence", "Write_Alignments", "Filling_Time", "Filtering", "Unique_Time",
                                "Extract_Info", "Extract_Events", "Extract_Lite", "Fill_Deletions",
                                "Update_Variables")

  # For each barcode
  for (i in 1:nrow(barcodeTable)){
    #------------------------------------------------------
    # INITIALIZE VARIABLES
    #------------------------------------------------------
      # Regular variables
      absolutePathForward      <- ""
      absolutePathReverse      <- ""
      forwardReadsFileName     <- ""
      reverseReadsFileName     <- ""
      dataFolderPath           <- ""
      emptyUniques             <- FALSE
      totalSequenceLost        <- 0
      totalSequencesPrefilter  <- 0
      totalSequencesPostfilter <- 0
      totalUniqueSequences     <- 0
      totalSequences           <- 0
      barcode                  <- barcodeTable$Barcode[i]
      equalPrimersIndex        <- 1
      subConfigLength          <- 0
      uniqueIndex              <- 1
      configSubsetIndex        <- 1


      # Dataframes and matrixes
      keepThisRows             <- NULL
      uniqueTable              <- NULL
      forwardsTable            <- NULL
      reversesTable            <- NULL
      configSubset             <- NULL
      forwardFound             <- NULL
      reverseFound             <- NULL
      asignedSequences         <- NULL

    #------------------------------------------------------
    # INITIALIZE THE CONFIG ROW VARIABLES
    # We need to do this everytime we change the config row. This will be repeated in the next loop
    #------------------------------------------------------
      # Progress feedback for the user
      print(paste("The processor",processID,"is doing line",configIndex,"of",nrow(configTable),
                  "Progress of: ",round(configIndex/nrow(configTable)*100,2),"%"))
      print(paste(i,"/",nrow(barcodeTable)))

      #------------------------
      # CHECK FOR BASIC INFO
      #------------------------
        # Initialize Variables
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

        # Specials IDs
        currentID                          <- configTable[configIndex,"ID"]
        currentBarcode                     <- configTable[configIndex,"Barcode"] # Should be equal to the variable: barcode

        # Primers and amplicon info
        forwardPrimer                      <- toString(configTable[configIndex,"Forward_Primer"])
        reversePrimer                      <- toString(configTable[configIndex,"Reverse_Primer"])
        targetPrimer                       <- toString(configTable[configIndex,"Target_Primer"])
        reversePrimerReverseComplementary  <- reverseComplement(reversePrimer)
        amplicon                           <- toString(configTable[configIndex,"Genome"])
        ampliconLength                     <- nchar(amplicon)
        forwardPrimerLength                <- nchar(forwardPrimer)
        reversePrimerLength                <- nchar(reversePrimer)

        # Names of files and folders
        currentIDCodeResultsFolderName     <- paste(alignmentFolder,"/", currentID,"_",currentBarcode,sep = '')
        currentAlignmentsResultsFolderName <- paste(currentIDCodeResultsFolderName,"/alignments",sep='')
        masterAlignmentFilePath            <- paste(currentIDCodeResultsFolderName,"/alignments.txt",sep='')
        uberAlignmentFilePath              <- paste(currentIDCodeResultsFolderName,"/verbose.txt",sep='')
        deletionFilePath                   <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                    currentBarcode,"_deletions.txt",sep = '')
        deletionRelativeFilePath           <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                    currentBarcode,"_deletionsRelative.txt",sep = '')
        uniqueTablePath                    <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                    currentBarcode,"_uniques.txt",sep = '')
        # Flags
        targetFound  <- FALSE
        primersFound <- FALSE
        apFound      <- FALSE

        # File descriptors
        masterFileConn     <- NULL
        uberAlignmentFD    <- NULL
        uniqueTableFD      <- NULL

      # Create folders
      dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
      dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)

      # Update the target Primer if necessary
      reverseAmplicon <- configTable[configIndex,"Strand"]
      if(reverseAmplicon == 1){
        targetPrimer  <- reverseComplement(targetPrimer)
      }

      # Get the position of the alignment and how wide it is
      # TODO: Adjust for more than one
      allPositions       <- countUppercaseGroups(amplicon)
      alignmentPositions <- allPositions[[2]][1]
      widePosition       <- allPositions[[3]][1]
      #------------------------
      # PREPARE THE FILE DESCRIPTORS
      #------------------------
        # If we need to write the alignment, prepare the proper file descriptors
        if(write_alignments > 0){
          if(write_alignments >= 1){
            # Prepare the summary alignment file
            masterFileConn     <- file(masterAlignmentFilePath, open="at")
            masterFileConnOpen <- TRUE
            writeLines(c(toString(amplicon),"\n"), masterFileConn)
          }
          if(write_alignments >= 2){
            # Start the file where the alignments are accumulated
            uberAlignmentFD     <- file(uberAlignmentFilePath, open="at")
            uberAlignmentFDOpen <- TRUE
          }
        }

        # Prepare the unique table where all the sequences are going.
        uniqueTableFD      <- file(uniqueTablePath, open="at")

        writeLines( paste("Total", "Frequency", "Forward", "Reverse", "Forward_Deletions_String",
                          "Reverse_Deletions_String", "Forward_Deletions_Genome_Coordinates_String",
                          "Reverse_Deletions_Genome_Coordinates_String",
                          "Forward_Alignment_Length", "Reverse_Alignment_Length",	"Alignment",
                          "Forward_Alignment_Positions", "Reverse_Alignment_Positions", "Forward_Wide_Position",
                          "Reverse_Wide_Position", "Target_Forward_Found",
                          "Target_Reverse_Found", "Deletions_Forward", "Deletions_Reverse",
                          sep = '\t')
                    , uniqueTableFD)
      #------------------------
      # CHECK FOR WARNINGS
      #------------------------
        # If the target is not in the amplicon, tell it to the user.
        targetFound  <- checkTarget(targetPrimer, amplicon, currentID, currentBarcode, configFilePath, logFileConn)

        # If the primers are not in the amplicon tell it to the user.
        primersFound <- checkPrimers(forwardPrimer, reversePrimerReverseComplementary, amplicon, currentID,
                                     currentBarcode, configFilePath, logFileConn)

        # If the amplicon has no alignment positions, tell the user.
        apFound      <- checkPositions(alignmentPositions, amplicon, currentID, currentBarcode, configFilePath,
                                       logFileConn)
    #------------------------------------------------------
    # GET THE FILES NAMES
    # - We need to adjust for relative or absolute paths
    # - We are just checkig for files names for both forwards and reverse
    #------------------------------------------------------
      # Check that we are in an absolute path or relative path
      # (Not ./folder/folder/file.txt path, but /folfer/file.txt or file.txt , something not unixy)
      absolutePathForward <- length(grep("/", barcodeTable$"Forward_Reads_File"[i], ignore.case = TRUE)) >= 1
      absolutePathReverse <- length(grep("/", barcodeTable$"Reverse_Reads_File"[i], ignore.case = TRUE)) >= 1

      # Get the names of the forward read and reverse read files
      if (absolutePathForward == TRUE) {
        forwardReadsFileName <- barcodeTable$"Forward_Reads_File"[i]
      } else {
        dataFolderPath <- getConfigFolder(configFilePath)
        forwardReadsFileName <- paste(dataFolderPath, "/",barcodeTable$"Forward_Reads_File"[i],sep = '')
      }

      if (absolutePathReverse == TRUE) {
        reverseReadsFileName <- barcodeTable$"Reverse_Reads_File"[i]
      } else {
        dataFolderPath <- getConfigFolder(configFilePath)
        reverseReadsFileName <- paste(dataFolderPath, "/",barcodeTable$"Reverse_Reads_File"[i],sep = '')
      }
    #------------------------------------------------------
    # GET THE CONTENT FOR THE FORWARD AND REVERSE FILES
    # - The content from the files will be in forwardsTable and reversesTable
    #------------------------------------------------------
      #---------------------------------------------------------------------------------
      # Get the actual FORWARD read file, unzip it, read it, an put it into a dataframe.
      #---------------------------------------------------------------------------------
      forwardsTable <- getReadsFile(forwardReadsFileName, temp_folder, tempFileConn, 0, 0)
      #---------------------------------------------------------------------------------
      # Get the actual REVERSE read file, unzip it, read it, an put it into a dataframe.
      #---------------------------------------------------------------------------------
      reversesTable <- getReadsFile(reverseReadsFileName, temp_folder, tempFileConn, 0, 0)
    #------------------------------------------------------
    # APPLY THE FILTERS
    #
    # Now, here we are going to apply the filters.
    # We are going to delete the rows in Forward and Reverse that didn't follow our filtering rules
    #------------------------------------------------------
      keepThisRows <- rep(TRUE,nrow(forwardsTable)) #In here we mark the bad rows to FALSE

      # Filter the sequences with bad nucleotides
      if (skip_bad_nucleotides == TRUE) {
        keepThisRows <- filterByNucleotides(forwardsTable, reversesTable, keepThisRows)
      }

      # Filter the sequences with bad quality
      if (min_quality > 0) {
        keepThisRows <- filterByQuality(forwardsTable, reversesTable, keepThisRows, min_quality, average_quality)
      }

      # Count how many sequences we have lost to the filters fascist uberzealous control system.
      totalSequenceLost        <- nrow(forwardsTable) - sum(keepThisRows)
      totalSequencesPrefilter  <- nrow(forwardsTable)
      totalSequencesPostfilter <- sum(keepThisRows)

      # Filter away the bad rows from the tables
      forwardsTable <- forwardsTable[keepThisRows,]
      reversesTable <- reversesTable[keepThisRows,]
    #------------------------------------------------------
    # MAKE THE UNIQUE TABLE
    #
    # Make a new table where we have the unique combination of forward and reverse reads.
    # In this table we also have the original total of those combination and the relative frequency.
    # Lastly; we are going to keep track of which config line uses this combo.
    #
    # "Total","Frequency","Forward", "Reverse"
    #------------------------------------------------------
      # Get the sequence columns from the forward and reverse table
      forwardsSubset <- subset(forwardsTable, select = c("Sequence"))
      reversesSubset <- subset(reversesTable, select = c("Sequence"))

      # Put them toguether
      uniqueTable <- cbind(forwardsSubset, reversesSubset)
      colnames(uniqueTable) <- c("Forward","Reverse")

      # Find out the unique and the total for each one; if we have no candidates, we can skip the whole function almost completely
      if (nrow(uniqueTable) > 0) {
        uniqueTable <- ddply(uniqueTable,.(uniqueTable$"Forward", uniqueTable$"Reverse"),nrow)
      } else {
        emptyUniques <- TRUE
        uniqueTable <- rbind(uniqueTable, c(NA,NA,1)) # We are doing a rbind of size 1 because magic
      }
      colnames(uniqueTable) <- c("Forward","Reverse","Total")

      # Prepare the relative sequence column, the user ID, and swap columns around
      uniqueTable$Frequency      <- 0
      uniqueTable <- uniqueTable[c("Total","Frequency","Forward", "Reverse")]

      # Sort the unique table by forward sequence and reverse sequence
      uniqueTable <- uniqueTable[order(uniqueTable$Forward, uniqueTable$Reverse),]

      # How many uniques we have is not the sum of the totals, but how many each unique without repetition we have
      totalUniqueSequences <- nrow(uniqueTable)

      # How many sequences in total we have
      totalSequences <- sum(uniqueTable$Total)

      #CANT MAKE THIS WORK WITH MUTATE/TRANSFORM!! WHY!????
      for (k in 1:nrow(uniqueTable)){
        uniqueTable$Frequency[k] <- uniqueTable$Total[k] / totalSequences
      }
    #------------------------------------------------------
    # Update the barcode table
    #------------------------------------------------------
      # Barcode stuff
      barcodeTable$Total_Pre_N_Filter[i]          <- totalSequencesPrefilter
      barcodeTable$Total_Post_N_Filter[i]         <- totalSequencesPostfilter
      barcodeTable$Experiment_Unique_Sequences[i] <- totalUniqueSequences
      barcodeTable$Total_Experiment_Sequences[i]  <- totalSequences
    #------------------------------------------------------
    # PREPARE THE UNASSIGNED VARIABLES
    #------------------------------------------------------
      # Also, we are going to keep a list for the unique sequences that are not in use
      # We are going to keep TRUE those which are unnasigned, so we don't need to complement this array later
      unnasignedSequences <- rep(TRUE, totalUniqueSequences)
      forwardFound        <- rep(FALSE, totalUniqueSequences)
      reverseFound        <- rep(FALSE, totalUniqueSequences)

      # We want to save the info in the barcode folder under the name <barcode>_unassigned.txt
      # Create the file and write the columns names
      unnasignedFileName <- paste(currentUnassignedFolderName,"/",barcode,"_unassigned.txt",sep = '', collapse = '')
    #------------------------------------------------------
    # PREPARE THE CONFIG LINES FOR THAT BARCODE
    #------------------------------------------------------
      # Get the rows in the config table with this barcode; remember that is sorted by barcode and primers
      configSubset    <- subset(configTable, Barcode == barcode)
      subConfigLength <- nrow(configSubset)

      # Now, keep tracking of where are we. We have two list to track, the unique table, and the config Table
      uniqueIndex       <- 1
      configSubsetIndex <- 1

      # Remember that we have this variable. It may happen that a config table have two duplicate entries
      # with, the same barcode and the same primers. If that happen, we need to pass throught the unique
      # table twice. So we use this variable to keep track of that, and roll back to a previous index
      # if so is needed.
      equalPrimersIndex <- uniqueIndex

      # Finally, we have another index that takes care of the overall config table. We start at that point
      # and we finish at that point plus the size of configSubset
    #------------------------------------------------------
    # LOOK FOR PRIMERS
    #
    # Repeat until we reach and process the end of the configSubset
    # Esentially, we are transversing the uniqueTable once.
    # However, because reasons (later), instead of making a for loop for each sequence in uniqueTable
    # we are using this while.
    # In most cases, there is a giagantic chuck of unnasigned that we don't need to process.
    #
    # You also need to take care of the emptyUniques, it might happens (rare) that we don't have
    # any sequence in the files, or none of them passed the filters. TODO
    #------------------------------------------------------
    if(emptyUniques == FALSE){
      while(configSubsetIndex <= subConfigLength) {
        # print(paste("Unique Index: ",uniqueIndex))
        # Search for the forward, reverse and targets
          # Initialize the variables for that sequence
          forwardFound       <- FALSE
          reverseFound       <- FALSE
          targetFoundForward <- FALSE
          targetFoundReverse <- FALSE
          # Get the sequences
          candidateForwardSequence <- toString(uniqueTable$Forward[uniqueIndex])
          candidateReverseSequence <- toString(uniqueTable$Reverse[uniqueIndex])
          # Search for the forward primer at the begginning of the forward read sequence
          forwardPrimerPositions <- regexpr(forwardPrimer, candidateForwardSequence, ignore.case=TRUE, fixed=FALSE)[1]
          if(forwardPrimerPositions > 0){
            forwardFound <- TRUE
          }
          # Search for the reverse primer at the beggining of the reverse read sequence
          reversePrimerPositions <- regexpr(reversePrimer, candidateReverseSequence, ignore.case=TRUE, fixed=FALSE)[1]
          if(reversePrimerPositions > 0){
            reverseFound <- TRUE
          }
          # Search for the target in the forward read sequence.
          targetPrimerPositions <- regexpr(targetPrimer, candidateForwardSequence, ignore.case=TRUE, fixed=FALSE)[1]
          if(targetPrimerPositions > 0 ){
            targetFoundForward <- TRUE
          }
          # Search for the target in the reverse read sequence.
          targetPrimerPositions <- regexpr(reverseComplement(targetPrimer), candidateReverseSequence, ignore.case=TRUE, fixed=FALSE)[1]
          if(targetPrimerPositions > 0){
            targetFoundReverse <- TRUE
          }
        # Addjust the finding of primers. If we are not going to use either the
        # forward sequences, or the reverse sequences, finding or not that primer
        # is irrelevant.
        if (fastqfiles == 1) {
          reverseFound <- TRUE
        } else {
          if(fastqfiles == 2){
            forwardFound <- TRUE
          }
        }

        # If the forward and reverse are found, do the alignments and everything else
        if (forwardFound && reverseFound) {
          # print("primers found")
          # Mark this one as assigned
          unnasignedSequences[uniqueIndex] <- FALSE

          # Find the reverse complementary of the reverse sequence for the multiple alignment
          candidateComplementarySequence <- reverseComplement(candidateReverseSequence)

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
          forwardAllPositions       <- countUppercaseGroups(toString(forwardSubjectAlignment))
          forwardAlignmentPositions <- forwardAllPositions[[2]][1]
          forwardWidePosition       <- forwardAllPositions[[3]][1]

          reverseAllPositions       <- countUppercaseGroups(toString(reverseSubjectAlignment))
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
                print(paste("The processor",processID,"is doing line",configIndex,"of",nrow(configTable),
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
                  currentBarcode                     <- configTable[configIndex,"Barcode"] # Should be equal to the variable: barcode

                  forwardPrimer                      <- toString(configTable[configIndex,"Forward_Primer"])
                  reversePrimer                      <- toString(configTable[configIndex,"Reverse_Primer"])
                  targetPrimer                       <- toString(configTable[configIndex,"Target_Primer"])
                  reversePrimerReverseComplementary  <- reverseComplement(reversePrimer)
                  amplicon                           <- toString(configTable[configIndex,"Genome"])
                  ampliconLength                     <- nchar(amplicon)
                  forwardPrimerLength                <- nchar(forwardPrimer)
                  reversePrimerLength                <- nchar(reversePrimer)

                  currentIDCodeResultsFolderName     <- paste(alignmentFolder,"/",currentID,"_",currentBarcode,sep = '')
                  currentAlignmentsResultsFolderName <- paste(currentIDCodeResultsFolderName,"/alignments",sep='')
                  masterAlignmentFilePath            <- paste(currentIDCodeResultsFolderName,"/alignments.txt",sep='')
                  uberAlignmentFilePath              <- paste(currentIDCodeResultsFolderName,"/verbose.txt",sep='')
                  deletionFilePath                   <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                              currentBarcode,"_deletions.txt",sep = '')
                  deletionRelativeFilePath           <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                              currentBarcode,"_deletionsRelative.txt",sep = '')
                  uniqueTablePath                    <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",
                                                              currentBarcode,"_uniques.txt",sep = '')

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
                    if(masterFileConnOpen){
                      close(masterFileConn)
                    }
                    if(uberAlignmentFDOpen){
                      close(uberAlignmentFD)
                    }
                  # File descriptors
                  masterFileConn     <- NULL
                  uberAlignmentFD    <- NULL
                  #                 deletionFD         <- NULL
                  #                 deletionRelativeFD <- NULL
                  uniqueTableFD      <- NULL
                # Create folders
                dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
                dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)
                # Update the target Primer if necessary
                reverseAmplicon <- configTable[configIndex,"Strand"]
                #               print(reverseAmplicon)
                if(reverseAmplicon == 1){
                  targetPrimer  <- reverseComplement(targetPrimer)
                }
                # Get the position of the alignment and how wide it is
                # TODO: Adjust for more than one
                allPositions       <- countUppercaseGroups(amplicon)
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
                      masterFileConnOpen <- TRUE
                      writeLines(c(toString(amplicon),"\n"), masterFileConn)
                    }
                    if (write_alignments >= 2) {
                      # Start the file where the alignments are accumulated
                      uberAlignmentFD     <- file(uberAlignmentFilePath, open="at")
                      uberAlignmentFDOpen <- TRUE
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
                  targetFound  <- checkTarget(targetPrimer, amplicon, currentID, currentBarcode, configFilePath,
                                              logFileConn)
                  # If the primers are not in the amplicon tell it to the user.
                  primersFound <- checkPrimers(forwardPrimer, reversePrimerReverseComplementary, amplicon,
                                               currentID, currentBarcode, configFilePath, logFileConn)
                  # If the amplicon has no alignment positions, tell the user.
                  apFound      <- checkPositions(alignmentPositions, amplicon, currentID, currentBarcode,
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

    #------------------------------------------------------
    # WRITE DOWN THE UNNASIGNED
    #------------------------------------------------------
      # Now that we are finish we know which sequences are unnasigned and which don't
      # We just need to write a subset of the uniqueTable
      # But first we are going to include the trackking of the forward and reverse primers that were found
      uniqueTable["Forward_Found"] <- forwardFound
      uniqueTable["Reverse_Found"] <- reverseFound
      uniqueTable <- uniqueTable[unnasignedSequences,]
      write.table(uniqueTable, file = unnasignedFileName, quote = FALSE, sep = "\t")

    # Close all the file descriptors that we had open so far at barcode level
    # These one is open for sure
    close(uniqueTableFD)

    # These two might not be needed, check it out and close if necessary
    if (masterFileConnOpen) {
      close(masterFileConn)
    }
    if (uberAlignmentFDOpen) {
      close(uberAlignmentFD)
    }
  }

  # Close all the file descriptors that we had open so far at function level
  close(logFileConn)
  close(tempFileConn)

  #------------------------------------------------------
  # WRITE THE BARCODE AND CONFIG TABLE
  #------------------------------------------------------
    write.table(configTable, paste(alignmentFolder, "/", processID, "_configFile_results", sep = '') , sep="\t")
    write.table(barcodeTable, paste(alignmentFolder, "/", processID, "_barcodeFile_results", sep = ''), sep="\t")
  return(0)
}
