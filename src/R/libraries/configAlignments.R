
source("libraries/tools2.R") # Minor stuff like the reverse complement of a DNA sequence.


makeAlignment5 <- function(SKIP_BAD_NUCLEOTIDES = TRUE, AVERAGE_QUALITY = 0,
                           MIN_QUALITY = 0, WRITE_ALIGNMENTS = TRUE,
                           SCORING_MATRIX = "NUC44", GAP_OPENING = 50,
                           GAP_EXTENSION = 0, GAP_ENDING = FALSE,
                           FAR_INDELS = TRUE, configDataframe = NULL,
                           resultFolder = NULL, alignmentFolder = NULL, 
                           processID = 0, configFilePath = "", TIMING){

  # Feedback with the processor that is working now
  print(paste("Nice to meet you; I'm mighty processor number: ",processID))
  
  # Prepare the unnasigned folder
  sourceCpp("libraries/gotoh3.cpp")
  currentUnassignedFolderName <- paste(alignmentFolder,"/unassigned_sequences",sep = '') # This dir is created previosly
  
  # In this variable we are going to keep track of the last config row. This is important
  # because it is going to be sorted. You need to be very careful with this function
  configIndex <- 1
  
  #------------------------------------------------------
  # SET THE FILES DESCRIPTOR
  #------------------------------------------------------
  {
    # Make your own file descriptor for the errors log
    subLogFileName <- paste(resultFolder,"/", processID,"_SUBLOG.txt",sep = '')
    logFileConn    <-  file(subLogFileName, open="at")
    
    # Make your own file descriptor for the final config file, and add the columns name
#     subConfigFileName <- paste(alignmentFolder,"/", processID,"_SUBCONFIG.txt",sep = '')
#     configFileConn    <-  file(subConfigFileName, open="at")
#     writeLines("'ID'\t'Barcode'\t'Forward_Reads_File'\t'Reverse_Reads_File'\t'Experiment_Type'\t'Target_Primer'\t'Forward_Primer'\t'Reverse_Primer'\t'Strand'\t'Amplicon'\t'Total_Pre_N_Filter'\t'Total_Post_N_Filter'\t'Experiment_Unique_Sequences'\t'Total_Experiment_Sequences'\t'Sum_Forward_Found'\t'Sum_Reverse_Found'\t'Sum_Target_Forward_Found'\t'Sum_Target_Reverse_Found'\t'Sum_Deletions_Forward'\t'Sum_Deletions_Reverse'\t'Sum_Is_Sequence'", configFileConn)    
    
    # These are for the configs later one; we need to keep track of them here
    # (Actually we don't, because these variables stay alive and R doesn't care about variable scope
    #  but I'm going to do it anyway; it should be easy to adapt this into another language later on)
    masterFileConn     <- NULL
    uberAlignmentFD    <- NULL
#     deletionFD         <- NULL
#     deletionRelativeFD <- NULL
    uniqueTableFD      <- NULL
    
    # Depending on the flags, these two FDs might not open; we keep track if they are open or not
    masterFileConnOpen     <- FALSE
    uberAlignmentFDOpen    <- FALSE
    
  }

  # Lets rename the configDataframe to configTable (so we don't make a lot of copypastes for the time being)
  configTable <- configDataframe

  # Sort the config table by barcode first, forward primer second, and reverse primer third
  configTable <- configTable[order(configTable$Barcode, configTable$Forward_Primer, configTable$Reverse_Primer),]

  # Add a bunch of columns to the dataframe of the config file that will serve as a summary
  # Please refer to the /doc folder for a verbose explanation for each one of these
  {
    # Barcode stuff
    configTable$Total_Pre_N_Filter <- 0
    configTable$Total_Post_N_Filter <- 0
    
    configTable$Experiment_Unique_Sequences <- 0
    configTable$Total_Experiment_Sequences  <- 0
    
    # Several statistics about deletions, cuts and reads 
#     configTable$Sum_Forward_Found        <- 0
#     configTable$Sum_Reverse_Found        <- 0
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
  }
  
  # Get the barcode subdataframe, this is the barcode and FWD and RVS files
  # Simplify the rows so we don't have duplicated.
  barcodeTable   <- unique(configTable[,c("Barcode","Forward_Reads_File", "Reverse_Reads_File")])

  # Add a bunch of columns to the barcode table. Although this info is already on the config, it actually belongs
  # to the barcodes
  {
    
    # Barcode stuff
    barcodeTable$Total_Pre_N_Filter <- 0
    barcodeTable$Total_Post_N_Filter <- 0
    
    barcodeTable$Experiment_Unique_Sequences <- 0
    barcodeTable$Total_Experiment_Sequences  <- 0
    
  }
  
  #------------------------------------------------------
  # PREPARE SOME TIMING VARIABLES
  #------------------------------------------------------
  {
#     # Get the barcode timing as well. In here we are going to keep track of barcode wise timing.
#     barcodeTiming <- data.frame(matrix(NA, nrow = nrow(barcodeTable), ncol = 4))
#     colnames(barcodeTiming) <- c("Barcode", "Total_Time","Sorting_Time", "Something_Time")
#     
#     # If we want to time the run, initialize the timing dataframe
#     if(TIMING == TRUE){
#       barcodeTiming$Barcode        <- barcodeTable$Barcode
#       barcodeTiming$Total_Time     <- 0
#       barcodeTiming$Sorting_Time   <- 0
#     }
    
    # We are going to create another dataframe for timing stuff.
    # Each row of this dataFrame represent each row of the dataframe
    configTiming <- data.frame(matrix(0, nrow = nrow(configTable), ncol = 15))
    colnames(configTiming) <- c("Gene_ID", "Barcode", "Total_Time","Alignment_Time", "Search_Primers",
                                "Write_Sequence", "Write_Alignments", "Filling_Time", "Filtering", "Unique_Time",
                                "Extract_Info", "Extract_Events", "Extract_Lite", "Fill_Deletions",
                                "Update_Variables")
    
    # If we want to time the run, initialize the timing dataframe
    if(TIMING == TRUE){
      configTiming$Gene_ID          <- configTable$ID
      configTiming$Barcode          <- configTable$Barcode
      configTiming$Total_Time       <- 0
      configTiming$Alignment_Time   <- 0
      configTiming$Search_Primers   <- 0
      configTiming$Write_Sequence   <- 0
      configTiming$Write_Alignments <- 0
      configTiming$Filling_Time     <- 0
      configTiming$Filtering        <- 0
      configTiming$Unique_Time      <- 0
      configTiming$Extract_Info     <- 0
      configTiming$Extract_Events   <- 0
      configTiming$Extract_Lite     <- 0
      configTiming$Fill_Deletions   <- 0
      configTiming$Update_Variables <- 0

    }
    
    t1 <- Sys.time()
    t2 <- Sys.time()
    t3 <- Sys.time()
    t4 <- Sys.time()
    t5 <- Sys.time()
    t6 <- Sys.time()
    t7 <- Sys.time()
    t8 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    td2 <- as.numeric(t4-t3, units = "secs")
    
  }
  
  # For each barcode
  for (i in 1:nrow(barcodeTable)){
    
    #------------------------------------------------------
    # INITIALIZE VARIABLES
    #------------------------------------------------------
    {
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
      
      # Timing variables
#       searchPrimersTime        <- 0
      
    } 
    
    #------------------------------------------------------
    # INITIALIZE THE CONFIG ROW VARIABLES
    # We need to do this everytime we change the config row. This will be repeated in the next loop
    #------------------------------------------------------
    {
     
      # Progress feedback for the user
      print(paste("The processor",processID,"is doing line",configIndex,"of",nrow(configTable),
                  "Progress of: ",round(configIndex/nrow(configTable)*100,2),"%"))
      print(paste(i,"/",nrow(barcodeTable)))
      
      #------------------------
      # CHECK FOR BASIC INFO
      #------------------------
      {
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
        
        # File descriptors
        masterFileConn     <- NULL
        uberAlignmentFD    <- NULL
#         deletionFD         <- NULL
#         deletionRelativeFD <- NULL
        uniqueTableFD      <- NULL
        
        # Timing
        t3                   <- Sys.time()
        totalTime            <- 0
        alignTime            <- 0
        searchPrimers        <- 0
        writeSequence        <- 0
        writeAlignment       <- 0
        fillingTime          <- 0
        filterTime           <- 0
        uniqueTime           <- 0
        extractAlignmentInfo <- 0
        extractEvents        <- 0
        extractLite          <- 0
        fillDeletionFile     <- 0
        updateVariables      <- 0
        
      }
      
#       #------------------------
#       # DEBUG
#       #------------------------
#       {
#         print(currentID)
#       }
#       print("A")
      # Create folders
      dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
      dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)
      
      # Update the target Primer if necessary
      reverseAmplicon <- configTable[configIndex,"Strand"]      
      if(reverseAmplicon == 1){
        targetPrimer  <- reverseComplement(targetPrimer)
      }
#       print("B")
      # Get the position of the alignment and how wide it is
      # TODO: Adjust for more than one
      allPositions       <- countUppercaseGroups(amplicon)
      alignmentPositions <- allPositions[[2]][1]
      widePosition       <- allPositions[[3]][1]
      print("C")
      #------------------------
      # PREPARE THE FILE DESCRIPTORS
      #------------------------
      {
      
        # If we need to write the alignment, prepare the proper file descriptors
        if(WRITE_ALIGNMENTS > 0){
          
          if(WRITE_ALIGNMENTS >= 1){
            
            # Prepare the summary alignment file
            masterFileConn     <- file(masterAlignmentFilePath, open="at")
            masterFileConnOpen <- TRUE
            writeLines(c(toString(amplicon),"\n"), masterFileConn)
            
          }
          
          if(WRITE_ALIGNMENTS >= 2){
            
            # Start the file where the alignments are accumulated            
            uberAlignmentFD     <- file(uberAlignmentFilePath, open="at")
            uberAlignmentFDOpen <- TRUE
          }
          
        }
        
#         deletionFD         <- file(deletionFilePath, open="at")
#         deletionRelativeFD <- file(deletionRelativeFilePath, open="at")
        uniqueTableFD      <- file(uniqueTablePath, open="at")
        
#         writeLines( paste("Gene_ID", "Hash_ID", "Start","End","Freq","Type","Location", "Wide",sep = '\t') 
#                     , deletionFD)
#         writeLines( paste("Gene_ID", "Hash_ID", "Start","End","Freq","Type","Location", "Wide",sep = '\t') 
#                     , deletionRelativeFD)
        writeLines( paste("Total", "Frequency", "Forward", "Reverse", "Forward_Deletions_String",
                          "Reverse_Deletions_String", "Forward_Deletions_Genome_Coordinates_String",
                          "Reverse_Deletions_Genome_Coordinates_String",
                          "Forward_Alignment_Length", "Reverse_Alignment_Length",	"Alignment",
                          "Forward_Alignment_Positions", "Reverse_Alignment_Positions", "Forward_Wide_Position",
                          "Reverse_Wide_Position", "Target_Forward_Found",
                          "Target_Reverse_Found", "Deletions_Forward", "Deletions_Reverse",
                          sep = '\t') 
                    , uniqueTableFD)
        
      }      
#       print("D")
      #------------------------
      # CHECK FOR WARNINGS
      #------------------------
      {
        # If the target is not in the amplicon, tell it to the user.
        targetFound  <- checkTarget(targetPrimer, amplicon, currentID, currentBarcode, configFilePath, logFileConn)
        
        # If the primers are not in the amplicon tell it to the user.
        primersFound <- checkPrimers(forwardPrimer, reversePrimerReverseComplementary, amplicon, currentID,
                                     currentBarcode, configFilePath, logFileConn)
        
        # If the amplicon has no alignment positions, tell the user.
        apFound      <- checkPositions(alignmentPositions, amplicon, currentID, currentBarcode, configFilePath,
                                       logFileConn)      
      }
      
    }
    
    #------------------------------------------------------
    # GET THE FILES NAMES
    #------------------------------------------------------
    {
      # Check that we are in an absolute path or relative path
      # (Not ./folder/folder/file.txt path, but /folfer/file.txt or file.txt , something not unixy)
      absolutePathForward <- ( length(grep("/", barcodeTable$"Forward_Reads_File"[i], ignore.case = TRUE)) >= 1 )
      absolutePathReverse <- ( length(grep("/", barcodeTable$"Reverse_Reads_File"[i], ignore.case = TRUE)) >= 1 )
      
      # Get the names of the forward read and reverse read files
      if(absolutePathForward == TRUE){
        forwardReadsFileName <- barcodeTable$"Forward_Reads_File"[i]
      }
      else{
        dataFolderPath <- getConfigFolder(configFilePath)
        forwardReadsFileName <- paste(dataFolderPath, "/",barcodeTable$"Forward_Reads_File"[i],sep = '')
      }
      
      if(absolutePathReverse == TRUE){
        reverseReadsFileName <- barcodeTable$"Reverse_Reads_File"[i]
      }
      else{
        dataFolderPath <- getConfigFolder(configFilePath)
        reverseReadsFileName <- paste(dataFolderPath, "/",barcodeTable$"Reverse_Reads_File"[i],sep = '')
      }
    }
#     print("E")
    #------------------------------------------------------
    # GET THE CONTENT FOR THE FORWARD AND REVERSE FILES
    #------------------------------------------------------
    {
      # Time how much does it take to read the forward and reverse files
      t1 <- Sys.time()
      
      #---------------------------------------------------------------------------------
      # Get the actual FORWARD read file, unzip it, read it, an put it into a dataframe.
      #---------------------------------------------------------------------------------
      forwardsTable <- getReadsFile(forwardReadsFileName)
      
      #---------------------------------------------------------------------------------
      # Get the actual REVERSE read file, unzip it, read it, an put it into a dataframe.
      #---------------------------------------------------------------------------------
      reversesTable <- getReadsFile(reverseReadsFileName)
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      fillingTime <- td
    }
#     print("F")
    #------------------------------------------------------
    # APPLY THE FILTERS
    #
    # Now, here we are going to apply the filters.
    # We are going to delete the rows in Forward and Reverse that didn't follow our filtering rules
    #------------------------------------------------------
    {
      # Time how much does it takes to filter the sequences
      t1 <- Sys.time()        
      
      keepThisRows <- rep(TRUE,nrow(forwardsTable)) #In here we mark the bad rows to FALSE
      
      # Filter the sequences with bad nucleotides
      if(SKIP_BAD_NUCLEOTIDES == TRUE){ 
        keepThisRows <- filterByNucleotides(forwardsTable, reversesTable, keepThisRows)
      }
      
      # Filter the sequences with bad quality
      if(MIN_QUALITY>0){
        keepThisRows <- filterByQuality(forwardsTable, reversesTable, keepThisRows, MIN_QUALITY, AVERAGE_QUALITY)
      }
      
      # Count how many sequences we have lost to the filters fascist uberzealous control system.
      totalSequenceLost        <- nrow(forwardsTable) - sum(keepThisRows)
      totalSequencesPrefilter  <- nrow(forwardsTable)
      totalSequencesPostfilter <- sum(keepThisRows)
      
      # Filter away the bad rows from the tables
      forwardsTable <- forwardsTable[keepThisRows,]
      reversesTable <- reversesTable[keepThisRows,]
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      filterTime <- td
    }
#     print("G")
    #------------------------------------------------------
    # MAKE THE UNIQUE TABLE
    #
    # Make a new table where we have the unique combination of forward and reverse reads.
    # In this table we also have the original total of those combination and the relative frequency.
    # Lastly; we are going to keep track of which config line uses this combo.
    #
    # "Total","Frequency","Forward", "Reverse"
    #------------------------------------------------------
    {      
      t1 <- Sys.time()        
      
      # Get the sequence columns from the forward and reverse table
      forwardsSubset <- subset(forwardsTable, select = c("Sequence"))
      reversesSubset <- subset(reversesTable, select = c("Sequence"))
      
      # Put them toguether
      uniqueTable <- cbind(forwardsSubset, reversesSubset)
      colnames(uniqueTable) <- c("Forward","Reverse")
      
      # Find out the unique and the total for each one; if we have no candidates, we can skip the whole function almost completely
      if(nrow(uniqueTable) > 0){
        uniqueTable <- ddply(uniqueTable,.(uniqueTable$"Forward", uniqueTable$"Reverse"),nrow)
      }
      else{
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
      
      write.table(uniqueTable,  paste(currentResultsFolderName,"/",processID,"_DEBUGUNIQUE.txt",sep = '') , sep="\t")
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      uniqueTime <- td
      
    }    
#     print("H")
    #------------------------------------------------------
    # Update the barcode table
    #------------------------------------------------------
    {
    
      # Barcode stuff
      barcodeTable$Total_Pre_N_Filter[i]          <- totalSequencesPrefilter
      barcodeTable$Total_Post_N_Filter[i]         <- totalSequencesPostfilter
      barcodeTable$Experiment_Unique_Sequences[i] <- totalUniqueSequences
      barcodeTable$Total_Experiment_Sequences[i]  <- totalSequences
    
    }
#     print("I")
    #------------------------------------------------------
    # PREPARE THE UNASSIGNED VARIABLES
    #------------------------------------------------------
    {
      # Also, we are going to keep a list for the unique sequences that are not in use
      # We are going to keep TRUE those which are unnasigned, so we don't need to complement this array later
      unnasignedSequences <- rep(TRUE,totalUniqueSequences)
      forwardFound        <- rep(FALSE,totalUniqueSequences)
      reverseFound        <- rep(FALSE,totalUniqueSequences)
      
      # We want to save the info in the barcode folder under the name <barcode>_unassigned.txt
      # Create the file and write the columns names
      unnasignedFileName <- paste(currentUnassignedFolderName,"/",barcode,"_unassigned.txt",sep = '', collapse = '')      
    }
#     print("J")
    #------------------------------------------------------
    # PREPARE THE CONFIG LINES FOR THAT BARCODE
    #------------------------------------------------------
    {
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
    }
#     print("K")
    #------------------------------------------------------
    # DEBUG STUFF
    #------------------------------------------------------
    {
      print("DEBUG")
      print(paste("barcode", barcode))
#       print(paste("Forward file", forwardReadsFileName))
#       print(paste("Reverse file", reverseReadsFileName))
      
#       print(paste("Sequences Lost", totalSequenceLost))
#       print(paste("Prefilter", totalSequencesPrefilter))
#       print(paste("Postfilter", totalSequencesPostfilter))
      print(paste("Total Uniques", totalUniqueSequences))
#       print(paste("Total Sequences", totalSequences))
      
      print(paste("Config Subset Size", subConfigLength))
      print(paste("First config table row", configIndex))
      print(paste("Last config table row", configIndex + subConfigLength - 1))
      print(paste("", 0))
    }

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
    # any sequence in the files, or none of them passed the filters.
    #------------------------------------------------------
    if(emptyUniques == FALSE){
    while(configSubsetIndex <= subConfigLength) {
#       print(paste("Unique Index: ",uniqueIndex))
      # Search for the forward, reverse and targets
      {
        
        t1 <- Sys.time()
        
        # Initialize the variables for that sequence
        forwardFound       <- FALSE
        reverseFound       <- FALSE
        targetFoundForward <- FALSE
        targetFoundReverse <- FALSE
#         print("Look for candidates")
        # Get the sequences
        candidateForwardSequence <- toString(uniqueTable$Forward[uniqueIndex])
        candidateReverseSequence <- toString(uniqueTable$Reverse[uniqueIndex])
#         print("Look for primers")
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
#         print("Look for targets")
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

        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        searchPrimers <- searchPrimers + td
      }
      
      # If the forward and reverse are found, do the alignments and everything else
      if(forwardFound && reverseFound){
#         print("primers found")
        # Mark this one as assigned
        unnasignedSequences[uniqueIndex] <- FALSE

        # Find the reverse complementary of the reverse sequence for the multiple alignment
        candidateComplementarySequence <- reverseComplement(candidateReverseSequence)
        
        # Prepare the variables as a string. This is necessary for the input of the aligning function.
        forwardString   <- candidateForwardSequence
        reverseString   <- toString(candidateComplementarySequence)[1]
        ampliconString  <- amplicon
        
        #------------------------------------------------------
        # Align the candidates with the genome
        #------------------------------------------------------
        {
          # The return of gotoh is a string divided in five parts with the separator "++++"
          # -- The first part is the verbose alignment result
          # -- The second part is a string with the insertions and deletions and missmatches
          # -- The third part is a string with the insertions and deletions and missmatches but from the subject (amplicon) coordinates
          # -- The fourth part is the alignment of the pattern.
          # -- The fifth part is the alignment of the subject.
          t1 <- Sys.time()
          alignForward <- gotoh(forwardString, ampliconString , SCORING_MATRIX, GAP_OPENING, GAP_EXTENSION, GAP_ENDING, FAR_INDELS)
          alignReverse <- gotoh(reverseString, ampliconString , SCORING_MATRIX, GAP_OPENING, GAP_EXTENSION, GAP_ENDING, FAR_INDELS)
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          alignTime <- alignTime + td
        }

        #------------------------------------------------------
        # Extract the information from the alignments (each part to a different variable)
        #------------------------------------------------------
        {
          t1 <- Sys.time()
          
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
          
          forwardAlignmentLength <- nchar(forwardPatternAlignment)
          reverseAlignmentLength <- nchar(reversePatternAlignment)
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          extractAlignmentInfo <- extractAlignmentInfo + td
        }
        
        #------------------------------------------------------
        #Get the forward and reverse events coordinates
        #------------------------------------------------------
        {
          t1 <- Sys.time()
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
          
#           # Missmatches test
#           fMO <- forwardData[[7]]
#           fMM <- forwardData[[8]]
#           fMP <- forwardData[[9]]
#           print("Miss example")
#           print(forwardLiteString)
#           print(forwardData[[7]])
#           print(forwardData[[8]])
#           print(forwardData[[9]])
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          extractEvents <- extractEvents + td
        }
        
        #------------------------------------------------------
        # Extract the information from the lite string in genome coordinates
        #------------------------------------------------------
        {
          t1 <- Sys.time()
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
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          extractLite <- extractLite + td
        }
        
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
        {
          t1 <- Sys.time()
          
          # If we have an alignment, write it into a TXT file
          # The name of the alignment file is keep regarless of the option of writing on disk or not.
          # Is a nice visual aid to locate uniques with both primers without looking at the TRUE/FALSE of the primers.
          Alignment_File <- paste(uniqueIndex,".txt",sep='')  
          
          if(WRITE_ALIGNMENTS > 0){
            
            if(WRITE_ALIGNMENTS >= 1){
              
              # Write the summary master alignment
              writeLines(c(paste(uniqueIndex,toString(uniqueTable$Total[uniqueIndex])),forwardPatternAlignment,reversePatternAlignment), masterFileConn)                
              
            }
            
            if(WRITE_ALIGNMENTS >= 2){
              
              # Write the verbose alignment into an individual file
              alignmentName <- paste(uniqueIndex,".txt",sep='')
              
              # Write the uber alignment file
              writeLines(c(alignmentName,       "\n FORWARD AND AMPLICON: \n",
                           forwardSummaryString,"\n REVERSE AND AMPLICON: \n",
                           reverseSummaryString,"\n"), uberAlignmentFD)
              
            }
            
            if(WRITE_ALIGNMENTS == 3){
              
              Alignment_File_Path <- paste(currentAlignmentsResultsFolderName, "/", alignmentName, sep='')            
              fileConn<-file(Alignment_File_Path)
              writeLines(c(forwardSummaryString,"\n + \n",reverseSummaryString), fileConn)
              close(fileConn)
              
            }
            
          }
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          writeAlignment <- writeAlignment + td
          
        }
        
        #--------------------------------------------------
        # Fill the deletion table DEPRECRATED
        #--------------------------------------------------
        {
          t1 <- Sys.time()
          
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
              
#               writeLines( paste(currentID,
#                                 uniqueIndex,
#                                 as.numeric(forwardDeletionsStarts[m]),
#                                 as.numeric(addjustedEndForward),
#                                 as.numeric(uniqueTable[uniqueIndex,2]),
#                                 "deletion",
#                                 "Forward",
#                                 addjustedEndForward - forwardDeletionsStarts[m],
#                                 sep = '\t') , deletionFD)
#               
#               writeLines( paste(currentID,
#                                 uniqueIndex,
#                                 as.numeric(forwardGenomeDeletionsStarts[m]),
#                                 as.numeric(addjustedRelativeEndForward),
#                                 as.numeric(uniqueTable[uniqueIndex,2]),
#                                 "deletion",
#                                 "Forward",
#                                 addjustedRelativeEndForward - forwardGenomeDeletionsStarts[m],
#                                 sep = '\t') , deletionRelativeFD)
              
            }}
          
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
              
#               writeLines( paste(currentID,
#                                 k,
#                                 as.numeric(reverseDeletionsStarts[m]),
#                                 as.numeric(addjustedEndReverse),
#                                 -as.numeric(uniqueTable[uniqueIndex,2]),
#                                 "deletion",
#                                 "Reverse",
#                                 addjustedEndReverse - reverseDeletionsStarts[m],
#                                 sep = '\t') , deletionFD)
#               
#               writeLines( paste(currentID,
#                                 k,
#                                 as.numeric(reverseGenomeDeletionsStarts[m]),
#                                 as.numeric(addjustedRelativeEndReverse),
#                                 -as.numeric(uniqueTable[uniqueIndex,2]),
#                                 "deletion",
#                                 "Reverse",
#                                 addjustedRelativeEndReverse - reverseGenomeDeletionsStarts[m],
#                                 sep = '\t') , deletionRelativeFD)
              
            }}
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          fillDeletionFile <- fillDeletionFile + td
          
        }       
        
        #--------------------------------------------------
        # Write the sequence into disk
        # -- Each folder has a comprehensive list of all the sequences with some stats
        # -- This part add the sequence to its file.
        #--------------------------------------------------
        {
        
          t1 <- Sys.time()
          
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
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          
          writeSequence  <- writeSequence + td
        
        }
        
        #--------------------------------------------------
        # Update the config statistics
        #--------------------------------------------------
        {
          t1 <- Sys.time()
          
          sumTargetForwardFound     <- sumTargetForwardFound + (targetFoundForward * uniqueTable$Total[uniqueIndex])
          sumTargetReverseFound     <- sumTargetReverseFound + (targetFoundReverse * uniqueTable$Total[uniqueIndex])
          sumDeletionsForward       <- sumDeletionsForward   + (forwardDeletions   * uniqueTable$Total[uniqueIndex])
          sumDeletionsReverse       <- sumDeletionsReverse   + (reverseDeletions   * uniqueTable$Total[uniqueIndex])
          sumIsSequence             <- sumIsSequence         + uniqueTable$Total[uniqueIndex]        
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          updateVariables <- updateVariables + td
        }

        # Go forward one sequence
        uniqueIndex <- uniqueIndex + 1
        
      }
      # If either forward or reverse is not found
      else{
#         print("Not found")
#         print(configSubset$ID)
#         print("Sequences")
#         print(candidateForwardSequence)
#         print(candidateReverseSequence)
        # Make sure that we are still under the primer scope. Otherwise take the next primer
        # This part is quite complicated to understand; so just trust me.
        forwardPrimerLength <- nchar(forwardPrimer)
        reversePrimerLength <- nchar(reversePrimer)
#         print("Getting samples")
        forwardSequenceSample <- toupper(substr(candidateForwardSequence, 1, forwardPrimerLength))
        reverseSequenceSample <- toupper(substr(candidateReverseSequence, 1, reversePrimerLength))
#         print("Getting lengths")
        forwardSequenceLength <- nchar(candidateForwardSequence)
        reverseSequenceLength <- nchar(candidateReverseSequence)
#         print("Transform to upper")
#         print(forwardPrimer)
#         print(reversePrimer)
#         print(forwardSequenceLength)
#         print(reverseSequenceLength)
#         print(forwardPrimerLength)
#         print(reversePrimerLength)
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
#         print(forwardPrimerUpper)
#         print(forwardPrimerLower)
#         print(reversePrimerUpper)
#         print(reversePrimerLower)
#         print(toupper(forwardPrimerUpper))
#         print(toupper(forwardPrimerLower))
#         print(toupper(reversePrimerUpper))
#         print(toupper(reversePrimerLower))
        insideScope <- TRUE

        if (
          (toupper(forwardPrimerUpper) == candidateForwardSequence && toupper(reversePrimerUpper) < candidateReverseSequence)
          ||
          (toupper(forwardPrimerUpper) <  candidateForwardSequence)
          ||
          uniqueIndex == nrow(uniqueTable)
          ){
#           print("Out of scope")
          insideScope <- FALSE
          
        }
        
        # If we are inside the scope of the primers
        if(insideScope == TRUE){
#           print("Inside Scope")
          # Keep track on which primer was not found
          forwardFound[uniqueIndex] <- forwardFound
          reverseFound[uniqueIndex] <- reverseFound
          
          # Pick the next one
          uniqueIndex <- uniqueIndex + 1
        }
        
        # If we passed the scope of the current primers, take the next two primers
        else{
#           print("Recording configs")
          # Finish the timing for the last row, and add the timing variables to the timing dataframes          
          t4 <- Sys.time()
          totalTime <- as.numeric(t4-t3, units = "secs")

          configTiming$Total_Time[configIndex]       <- totalTime          
          configTiming$Alignment_Time[configIndex]   <- alignTime
          configTiming$Search_Primers[configIndex]   <- searchPrimers
          configTiming$Write_Sequence[configIndex]   <- writeSequence
          configTiming$Write_Alignments[configIndex] <- writeAlignment
          configTiming$Filling_Time[configIndex]     <- fillingTime
          configTiming$Filtering[configIndex]        <- filterTime
          configTiming$Unique_Time[configIndex]      <- uniqueTime
          configTiming$Extract_Info[configIndex]     <- extractAlignmentInfo
          configTiming$Extract_Events[configIndex]   <- extractEvents
          configTiming$Extract_Lite[configIndex]     <- extractLite
          configTiming$Fill_Deletions[configIndex]   <- fillDeletionFile
          configTiming$Update_Variables[configIndex] <- updateVariables
                  
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
              
            }
            # Otherwise, update the primer index to the current unique
            else{
              
              equalPrimersIndex <- uniqueIndex
              
            }
            
          }
          
          # You also need to restart the variables for the next row, if we have next row.
          # If you look inside here, you will notice that we close the FDs from the previous one.
          # Don't worry for the last one, because they are close outside the loop, just like the
          # first init is done outside the loop too.
          if(configSubsetIndex <= subConfigLength){
#             print("Another row to come")
            #------------------------------------------------------
            # INITIALIZE THE CONFIG ROW VARIABLES
            # The same as before, we in here we are also closing the FD if needed
            #------------------------------------------------------
            {
              
              # Progress feedback for the user
              print(paste("The processor",processID,"is doing line",configIndex,"of",nrow(configTable),
                          "Progress of: ",round(configIndex/nrow(configTable)*100,2),"%"))
#               print("M")
              #------------------------
              # CHECK FOR BASIC INFO
              #------------------------
              {
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
                
                # Timing
                t3                   <- Sys.time()
                totalTime            <- 0
                alignTime            <- 0
                searchPrimers        <- 0
                writeSequence        <- 0
                writeAlignment       <- 0
                fillingTime          <- 0
                filterTime           <- 0
                uniqueTime           <- 0
                extractAlignmentInfo <- 0
                extractEvents        <- 0
                extractLite          <- 0
                fillDeletionFile     <- 0
                updateVariables      <- 0
                
                #------------------------
                # CLOSE PREVIOUS FILE DESCRIPTORS
                #------------------------
                {
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
                }
            
                # File descriptors
                masterFileConn     <- NULL
                uberAlignmentFD    <- NULL
#                 deletionFD         <- NULL
#                 deletionRelativeFD <- NULL
                uniqueTableFD      <- NULL
              }
#               print("N")
#               print(currentIDCodeResultsFolderName)
#               print(currentAlignmentsResultsFolderName)
              # Create folders
              dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
              dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)

              # Update the target Primer if necessary
              reverseAmplicon <- configTable[configIndex,"Strand"]      
#               print(reverseAmplicon)
              if(reverseAmplicon == 1){
                targetPrimer  <- reverseComplement(targetPrimer)
              }
#               print(reverseAmplicon)
#               print(amplicon)
              # Get the position of the alignment and how wide it is
              # TODO: Adjust for more than one
              allPositions       <- countUppercaseGroups(amplicon)
              alignmentPositions <- allPositions[[2]][1]
              widePosition       <- allPositions[[3]][1]
#               print("O")
              #------------------------
              # PREPARE THE FILE DESCRIPTORS
              #------------------------
              {
                
                # If we need to write the alignment, prepare the proper file descriptors
                if(WRITE_ALIGNMENTS > 0){
                  
                  if(WRITE_ALIGNMENTS >= 1){
                    
                    # Prepare the summary alignment file
                    masterFileConn     <- file(masterAlignmentFilePath, open="at")
                    masterFileConnOpen <- TRUE
                    writeLines(c(toString(amplicon),"\n"), masterFileConn)
                    
                  }
                  
                  if(WRITE_ALIGNMENTS >= 2){
                    
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
                                  sep = '\t') 
                            , uniqueTableFD)                
                
              }      
#               print("P")
              #------------------------
              # CHECK FOR WARNINGS
              #------------------------
              {
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
  
            }
          
          }
#           print("Q")
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
        
#         print("Unique Index")
#         print(uniqueIndex)
#         
#         print("Forward and reverse")
#         print(toString(uniqueTable$Forward[uniqueIndex]))
#         print(toString(uniqueTable$Reverse[uniqueIndex]))
#         
#         print("Whole line")
#         print(unnasignedSequence)
        
      }
      
    } # While loop
    } # If unique is not empty

    #------------------------------------------------------
    # WRITE DOWN THE UNNASIGNED
    #------------------------------------------------------
    {
      # Now that we are finish we know which sequences are unnasigned and which don't
      # We just need to write a subset of the uniqueTable
      # But first we are going to include the trackking of the forward and reverse primers that were found
      uniqueTable["Forward_Found"] <- forwardFound
      uniqueTable["Reverse_Found"] <- reverseFound
      uniqueTable <- uniqueTable[unnasignedSequences,]
      write.table(uniqueTable, file = unnasignedFileName, quote = FALSE, sep = "\t")
    }

    # Close all the file descriptors that we had open so far at barcode level
    # These three are open for sure
#     close(deletionFD)
#     close(deletionRelativeFD)
    close(uniqueTableFD)
    # These two might not be needed, check it out and close if necessary
    if(masterFileConnOpen){
      close(masterFileConn)
    }
    if(uberAlignmentFDOpen){
      close(uberAlignmentFD)
    }

    # Write down the timing for the last row
    #configTiming$Total_Time[configIndex] <- totalTime

  }
  
  # Close all the file descriptors that we had open so far at function level
  close(logFileConn)
#   close(configFileConn)

  #------------------------------------------------------
  # WRITE THE TIMING VARIABLES
  #------------------------------------------------------
  {
    if(TIMING == TRUE){
      write.table(configTiming, paste(currentResultsFolderName,"/",processID,
                                      "_configFile_Timing",sep = ''), sep="\t")
#       write.table(barcodeTiming, paste(currentResultsFolderName,"/",processID,
#                                       "_barcodeFile_Timing",sep = ''), sep="\t")
    }
    
    
    
  }

  #------------------------------------------------------
  # WRITE THE BARCODE AND CONFIG TABLE
  #------------------------------------------------------
  {
    write.table(configTable,  paste(currentResultsFolderName,"/",processID,"_configFile_results",sep = '') , sep="\t")
    write.table(barcodeTable, paste(currentResultsFolderName,"/",processID,"_barcodeFile_results",sep = ''), sep="\t")
  }
  
  
  
  return (0)
}


