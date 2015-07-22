analyzeConfig <- function(ALL_IN_MODE = 3, ALIGNMENT_DISTANCE = 5, SAME_CUTS = TRUE, PRIMERDIMER = TRUE,
                          DIMERBASE = 10, BOTHDIMERS = FALSE, configDataframe = NULL,
                          analysisPath = NULL, processID = 0, resultsFolder = NULL, TIMING = FALSE ){
  

  
  # Get some timing variables ready
  {
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
    allTime <- 0                  # This variable is how long does it takes to process one row in the config table
    fillingTime <- 0             
    filterTime <- 0
    uniqueTime <- 0
    warningChecking <- 0
    readBasicData <- 0
    expandUniques <- 0
    processRow <- 0              # This is what to take to process one row, after the initialization is complete
    
    # Everything bellow this (including this) is the breakdown of what it is doing inside the row processing
    initializeVariables <- 0
    searchPrimers <- 0
    processRow_B <- 0            # If the primers are found, do more stuff, everything bellow breakdown this time
    
    alignTime <- 0
    extractAlignmentInfo <- 0
    extractEvents <- 0
    extractLite <- 0
    allInTiming <- 0
    sameTotalTiming <- 0
    sameCutTiming <- 0
    primerDimerTiming <- 0
    frameShiftTest <- 0
    applyCriterias <- 0
    writeAlignments <- 0
    
    fillUniqueTable <- 0
    archplotTiming <- 0
    histogramPlottingTiming <- 0
    recordStats <- 0
  }
  
  # Feedback with the processor that is working now
  print(paste("Nice to meet you; I'm mighty processor number: ",processID))

  # Lets rename the configDataframe to configTable (so we don't make a lot of copypastes for the time being)
  configTable <- configDataframe

  # Make your own file descriptor for the logs
  subLogFileName <- paste(analysisPath,"/", processID,"_SUBLOG.txt",sep = '')
  logFileConn <-  file(subLogFileName, open="at")

  # Prepare the timing dataframe
  configTiming <- data.frame(matrix(0, nrow = nrow(configTable), ncol = 15))
  colnames(configTiming) <- c("Gene_ID", "Barcode", "Total_Time","Time_ReadBasic", "Time_ProcessRow",
                              "Time_RecordStat", "Time_ExpandUniques", "Time_InitVariables",
                              "Time_ExtractAlignmentInfo", "Time_ALLIN", "Time_SAMECUT",
                              "Time_PRIMERDIMER", "Time_Frameshift", "Time_ApplyCriterias",
                              "Time_FillUnique")
  
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
  
  # Add a bunch of columns to the dataframe of the config file that will serve as a summary
  # Please refer to the /doc folder for a verbose explanation for each one of these
  {
  
    # Read Rate, what it is commonly call Cut Rate
    configTable$Cut_Rate <- 0
    
    # Deletions, valid or not
    configTable$Forward_No_Deletions       <- 0
    configTable$Forward_In_Frameshift      <- 0
    configTable$Forward_Out_Frameshift     <- 0
    configTable$Forward_Multiple_Deletions <- 0
    
    configTable$Reverse_No_Deletions       <- 0
    configTable$Reverse_In_Frameshift      <- 0
    configTable$Reverse_Out_Frameshift     <- 0
    configTable$Reverse_Multiple_Deletions <- 0
    
    # Deletions, ONLY valid
    configTable$Forward_Read_No_Deletions       <- 0 #This should always be zero, will keep this for debugging only
    configTable$Forward_Read_In_Frameshift      <- 0
    configTable$Forward_Read_Out_Frameshift     <- 0
    configTable$Forward_Read_Multiple_Deletions <- 0
    
    configTable$Reverse_Read_No_Deletions       <- 0
    configTable$Reverse_Read_In_Frameshift      <- 0
    configTable$Reverse_Read_Out_Frameshift     <- 0
    configTable$Reverse_Read_Multiple_Deletions <- 0
    
    # Only valid deletions
    configTable$Sum_Forward_Cuts <- 0
    configTable$Sum_Reverse_Cuts <- 0
    
    configTable$Sum_Is_Read <- 0               
               
    #-------------------INDEPENDENT-----------------------------------------
    configTable$Sum_Unique_Cuts_Passed_ALLIN        <- 0 # Independent cuts
    configTable$Sum_Unique_Cuts_Passed_SAMECUT      <- 0
    configTable$Sum_Unique_Cuts_Passed_PRIMERDIMER  <- 0
    
    configTable$Sum_Unique_Reads_Passed_ALLIN       <- 0 # Independent reads
    configTable$Sum_Unique_Reads_Passed_SAMECUT     <- 0
    configTable$Sum_Unique_Reads_Passed_PRIMERDIMER <- 0
    
    #-------------------SUCCESSIVELY-----------------------------------------
    configTable$Sum_Unique_Cuts_After_ALLIN        <- 0 # Successively cuts
    configTable$Sum_Unique_Cuts_After_SAMECUT      <- 0
    configTable$Sum_Unique_Cuts_After_PRIMERDIMER  <- 0
    
    configTable$Sum_Unique_Reads_After_ALLIN       <- 0 # Successively reads
    configTable$Sum_Unique_Reads_After_SAMECUT     <- 0
    configTable$Sum_Unique_Reads_After_PRIMERDIMER <- 0
  
  }

  # For each line of the config table
  for (j in 1:nrow(configTable)){
  
    # Set the times to zero
    {
      allTime          <- 0           
      readBasicData    <- 0
      processUnique    <- 0    
      recordStats      <- 0
      
      expandUniques        <- 0      
      initializeVariables  <- 0
      extractAlignmentInfo <- 0
      allInTiming          <- 0
      sameCutTiming        <- 0
      primerDimerTiming    <- 0
      frameShiftTiming     <- 0      
      applyCriterias       <- 0
      fillUniqueTable      <- 0      
    }
  
    # Start counting time before processing that config line, find the difference when is finish and write it in
    # the statistics.
    t3 <- Sys.time()

    # Initialize some variables
    currentLineTotalPrimers <- 0
    totalUniqueSequences <- 0

    #------------------------
    # CHECK FOR BASIC INFO
    #------------------------
    {
      t1 <- Sys.time()
    
      # Get the ID and barcode
      currentID           <- configTable[j,"ID"]
      currentBarcode      <- configTable[j,"Barcode"]
      forwardPrimerLength <- nchar(as.character(configTable[j,"Forward_Primer"]))
      reversePrimerLength <- nchar(as.character(configTable[j,"Reverse_Primer"]))
      ampliconLength      <- nchar(as.character(configTable[j,"Genome"]))
    
      # Progress feedback for the user
      print(paste("The processor",processID,"is doing line",j,"of",nrow(configTable), "Progress of: ",round(j/nrow(configTable)*100,2),"%"))
    
      # Make a folder for the results
      currentIDCodeResultsFolderName <- paste(analysisPath,"/",currentID,"_",currentBarcode,sep = '')
      dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
    
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      readBasicData <- td
    }
    
    # Get the file path for the unique table
    uniqueFilePath           <- paste(resultsFolder,'/',currentID,'_',currentBarcode,"/",currentID,'_',currentBarcode,
                                      "_uniques.txt", collapse='', sep='')

    # Check that the files exist and can be read
    uniqueReady           <- checkFileAccess(uniqueFilePath)

    if(uniqueReady == FALSE){

      if(uniqueReady == FALSE){
        print("-----------------------------------------")
        print("---|             WARNING!            |---")
        print("-----------------------------------------")
        print("WARNING!: Couldn't read this unique file.")
        print(uniqueFilePath)
        print("Does it exist? Do you have read access?")
        print("All warning and errors are register in the log")
        
        writeLines(paste("The processor ",processID), logFileConn)  
        writeLines(paste("Couldn't find the unique file ",uniqueFilePath), logFileConn)  
        writeLines("", logFileConn)  
        
      }
    
    }
    else{
      
      # Get the unique table into a dataframe
      uniqueTable <- read.table(uniqueFilePath, header=TRUE, sep="\t")
      
      # Keep the ID row for later, we will need it to compare with some info from the deletion table
      # Keep in mind that the the uniqueTable_i match the uniqueID_i
      uniqueIDs <- rownames(uniqueTable)
      
      # Redo the frequency. Now we don't care with respect the whole barcode file. Only with respect
      # the experiment.
      # How many sequences in total we have
      totalSequences <- sum(uniqueTable$Total)

      #CANT MAKE THIS WORK WITH MUTATE/TRANSFORM!! WHY!????
      if(nrow(uniqueTable)>0){
      for (k in 1:nrow(uniqueTable)){
        uniqueTable$Frequency[k] <- uniqueTable$Total[k] / totalSequences  
      }}

      cutsFilePath <- paste(currentIDCodeResultsFolderName,"/",currentID,"_",currentBarcode,
                                  "_cutsTable.txt",sep = '')

      cutsFD      <-  file(cutsFilePath, open="at")
        
      writeLines( paste("Gene_ID", "Hash_ID", "Start","End","Freq","Type","Location", "Wide",sep = '\t') 
                    , cutsFD)
      
      # If we have any candidates:
      t5 <- Sys.time()
      if(nrow(uniqueTable) > 0){
        
        #------------------------------------------------------
        # EXPAND THE UNIQUE TABLE
        #
        # Add empty colums where the rest of the info will go
        # This include the timing columns, deletions positions, and many other several statistics
        #------------------------------------------------------    
        {
          t1 <- Sys.time()
          
          uniqueTable$Forward_Frameshift <- 0
          uniqueTable$Reverse_Frameshift <- 0
          
          uniqueTable$Forward_Cuts <- 0
          uniqueTable$Reverse_Cuts <- 0
          
          uniqueTable$Is_Read <- FALSE             
          
          #-------------------INDEPENDENT-----------------------------------------
          uniqueTable$Unique_Cuts_Passed_ALLIN <- 0
          uniqueTable$Unique_Cuts_Passed_SAMECUT <- 0
          uniqueTable$Unique_Cuts_Passed_PRIMERDIMER <- 0
          
          #-------------------SUCCESSIVELY-----------------------------------------
          uniqueTable$Unique_Cuts_After_ALLIN <- 0
          uniqueTable$Unique_Cuts_After_SAMECUT <- 0
          uniqueTable$Unique_Cuts_After_PRIMERDIMER <- 0
          
          t2 <- Sys.time()
          td <- as.numeric(t2-t1, units = "secs")
          expandUniques <- td
        }  
        
        
        
        # For each of the candidates
        for (k in 1:nrow(uniqueTable)){
          
            t <- as.numeric(uniqueIDs[k]) # Crappy patchout from old code
#             print("K")
#             print(k)
#             print("T")
#             print(t)
#             print("F")
            #------------------------------------------------------
            # Initialize the variables for this unique combination
            #------------------------------------------------------
            {
              t1 <- Sys.time()  
              
              forwardLiteString <- as.character(uniqueTable$Forward_Deletions_String[k])
              reverseLiteString <- as.character(uniqueTable$Reverse_Deletions_String[k])
              
              forwardGenomeRelativeLiteString <- as.character(uniqueTable$Forward_Deletions_Genome_Coordinates_String[k])
              reverseGenomeRelativeLiteString <- as.character(uniqueTable$Reverse_Deletions_Genome_Coordinates_String[k])
              
              forwardFrameshiftClass <- 0
              reverseFrameshiftClass <- 0
              
              # All of these variables are for the cuts
              forwardCutsALLIN <- 0           #|
              reverseCutsALLIN <- 0           #| All of these are for the ALLIN criteria CUTS
              cutsALLIN        <- 0           #|  
              
              cutsSAMECUT <- 0                #|  Only this one is for the INTERVAL criteria, since its deals
                                              #|  with the forward and reverse at the same time. Keep in mind that
                                              #|  the total is still F + R; which are the same amount (ei: F x 2 )
              
              forwardCutsPRIMERDIMER <- 0     #|
              reverseCutsPRIMERDIMER <- 0     #| All of these are for the PRIMERDIMER criteria
              cutsPRIMERDIMER <- 0            #|  
              
              # These variables keep track of each data after each criteria is applied.
              cutsAfterALLIN       <- 0
              cutsAfterSAMECUT     <- 0
              cutsAfterPRIMERDIMER <- 0
              
              # All of these variables are for the reads
              readALLIN        <- FALSE       #| If the READ passed the ALLIN criteria , this is cutsALLIN > 0
              readSAMECUT      <- FALSE       #| If the READ passed the SAMECUT criteria , this is cutsSAMECUT > 0
              readPRIMERDIMER  <- FALSE       #| If the READ passed the PRIMERDIMER criteria , this is cutsPRIMERDIMER > 0
              
              readAfterALLIN        <- FALSE  #| If the READ passed the ALLIN criteria after the previous
              readAfterSAMECUT      <- FALSE  #| If the READ passed the SAMECUT criteria after the previous
              readAfterPRIMERDIMER  <- FALSE  #| If the READ passed the PRIMERDIMER criteria after the previous
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              initializeVariables <- initializeVariables + td
              
            }
            
            # -----------------------------------------------------
            # Extract the info from the alignments
            # -----------------------------------------------------
            {
              t1 <- Sys.time()
              
              # Get the alignment position, wide, and the length from forward and reverse
              forwardAlignmentPositions <- as.numeric(uniqueTable$Forward_Alignment_Positions[k])
              reverseAlignmentPositions <- as.numeric(uniqueTable$Reverse_Alignment_Positions[k])
              forwardWidePosition       <- as.numeric(uniqueTable$Forward_Wide_Position[k])
              reverseWidePosition       <- as.numeric(uniqueTable$Reverse_Wide_Position[k])
              forwardAlignmentLength    <- as.numeric(uniqueTable$Forward_Alignment_Length[k])
              reverseAlignmentLength    <- as.numeric(uniqueTable$Reverse_Alignment_Length[k])
              
              forwardData             <- getEventInfo(forwardLiteString)
              reverseData             <- getEventInfo(reverseLiteString)

              forwardDeletionsTotal   <- forwardData[[4]]
              forwardDeletionsStarts  <- forwardData[[5]]
              forwardDeletionsEnds    <- forwardData[[6]]
    
              reverseDeletionsTotal   <- reverseData[[4]]
              reverseDeletionsStarts  <- reverseData[[5]]
              reverseDeletionsEnds    <- reverseData[[6]]
              
              forwardGenomeData             <- getEventInfo(forwardGenomeRelativeLiteString)
              reverseGenomeData             <- getEventInfo(reverseGenomeRelativeLiteString)

              forwardGenomeDeletionsTotal   <- forwardGenomeData[[4]]
              forwardGenomeDeletionsStarts  <- forwardGenomeData[[5]]
              forwardGenomeDeletionsEnds    <- forwardGenomeData[[6]]

              reverseGenomeDeletionsTotal   <- reverseGenomeData[[4]]
              reverseGenomeDeletionsStarts  <- reverseGenomeData[[5]]
              reverseGenomeDeletionsEnds    <- reverseGenomeData[[6]]

              forwardRelativeDeletionsStarts <- forwardGenomeDeletionsStarts
              reverseRelativeDeletionsStarts <- reverseGenomeDeletionsStarts
              forwardRelativeDeletionsEnds   <- forwardGenomeDeletionsEnds
              reverseRelativeDeletionsEnds   <- reverseGenomeDeletionsEnds

              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              extractAlignmentInfo <- extractAlignmentInfo + td
              
            }

            # Now we are going to find out the total of valid deletion for each cut criteria which is activated.
            # It could happen that the firsts criterias don't eliminate anything and the last one kill them all. For that
            # cases we keep track of two different statistics. First we have a group of variable that keep track of the
            # result for one particular criteria. The other group of variable are updated after each criteria.
            
            # ------------------------------------------------------------------------------------------------        
            # Check that all events are within alignment position reach; or check that there is at least one.
            # ------------------------------------------------------------------------------------------------
            {
              t1 <- Sys.time()
  
              validForwardDeletionsALLIN  <- allInsideRange(forwardAlignmentPositions, forwardDeletionsStarts, forwardDeletionsEnds, ALL_IN_MODE, ALIGNMENT_DISTANCE, forwardWidePosition, FALSE, forwardAlignmentLength)
              validReverseDeletionsALLIN  <- allInsideRange(reverseAlignmentPositions, reverseDeletionsStarts, reverseDeletionsEnds, ALL_IN_MODE, ALIGNMENT_DISTANCE, reverseWidePosition, FALSE, reverseAlignmentLength)        
              
              # Track the stats for this criteria independently
              forwardCutsALLIN <- sum(validForwardDeletionsALLIN)      
              reverseCutsALLIN <- sum(validReverseDeletionsALLIN)      
              
              # Get the read results
              readALLIN        <- ((forwardCutsALLIN + reverseCutsALLIN)> 0)
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              allInTiming <- allInTiming + td
            }
            
            # --------------------------------------------------------        
            # Lets check the events have the same coordinates
            # --------------------------------------------------------
            {
              t1 <- Sys.time()
              
              intervalsResults <- sameCutHybrid(forwardRelativeDeletionsStarts, forwardRelativeDeletionsEnds,
                                                reverseRelativeDeletionsStarts, reverseRelativeDeletionsEnds,
                                                ampliconLength)
              
              # Extract the information from the result list
              validForwardDeletionsSAMECUTSINTERVALS <- intervalsResults[[1]]
              validReverseDeletionsSAMECUTSINTERVALS <- intervalsResults[[2]]
              
              # Track the stats for this criteria independently
              cutsSAMECUT <- sum(validForwardDeletionsSAMECUTSINTERVALS) + sum(validReverseDeletionsSAMECUTSINTERVALS)
              
              # Get the read results
              readSAMECUT              <- (cutsSAMECUT > 0)
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              sameCutTiming <- sameCutTiming + td
            }
       
            # --------------------------------------------------------        
            # Finding the primer dimer problem
            # -------------------------------------------------------- 
            {
              t1 <- Sys.time()          
              
              validForwardDeletionPRIMERDIMER <- primerDimerTest(forwardDeletionsStarts, forwardDeletionsEnds, forwardPrimerLength, reversePrimerLength, ampliconLength, DIMERBASE)
              validReverseDeletionPRIMERDIMER <- primerDimerTest(reverseDeletionsStarts, reverseDeletionsEnds, forwardPrimerLength, reversePrimerLength, ampliconLength, DIMERBASE)
              
              # Track the stats for this criteria independently
              forwardCutsPRIMERDIMER <- sum(validForwardDeletionPRIMERDIMER)        
              reverseCutsPRIMERDIMER <- sum(validReverseDeletionPRIMERDIMER)        
              
              # Get the read results
              readPRIMERDIMER              <- (cutsPRIMERDIMER > 0)
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              primerDimerTiming <- primerDimerTiming + td
            }
            
            # --------------------------------------------------------        
            # Test the deletion and see if they are in frameshift
            # --------------------------------------------------------
            {
              t1 <- Sys.time()
              
              forwardFrameshiftClass <- frameShift(forwardDeletionsStarts, forwardDeletionsEnds, forwardAlignmentLength)
              reverseFrameshiftClass <- frameShift(reverseDeletionsStarts, reverseDeletionsEnds, reverseAlignmentLength)
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              frameShiftTiming <- frameShiftTiming + td
            }
      
            #--------------------------------------------------
            # Apply all the criterias selected by the user
            #--------------------------------------------------
            {
              t1 <- Sys.time()
        
              # This arrays represent all the deletions in the forward and the revese alignment.
              # By default, all of them are set to TRUE because all of them are valid until a
              # criteria says otherwise.
              forwardCuts <- rep(TRUE,length(forwardDeletionsStarts))
              reverseCuts <- rep(TRUE,length(reverseDeletionsStarts))  
        
              # -- Apply that everything is inside rage. This criteria is applied always; you can choose to ignore it
              #    by changing the MODE to 5, or by giving a reach from AP or +500 (for example)
              forwardCuts <- forwardCuts * validForwardDeletionsALLIN
              reverseCuts <- reverseCuts * validReverseDeletionsALLIN
        
              cutsAfterALLIN       <- sum(forwardCuts) + sum(reverseCuts)
        
              readAfterALLIN       <- (cutsAfterALLIN > 0)
        
              # -- Now apply the criteria that says that the cuts must start and end at the same amplicon position
              if(SAME_CUTS == TRUE){
                forwardCuts <- forwardCuts * validForwardDeletionsSAMECUTSINTERVALS
                reverseCuts <- reverseCuts * validReverseDeletionsSAMECUTSINTERVALS
                cutsAfterSAMECUT <- sum(forwardCuts) + sum(reverseCuts)
              }
              else{
                cutsAfterSAMECUT <- cutsAfterSAMETOTAL
              }
              
              readAfterSAMECUT       <- (cutsAfterSAMECUT > 0)
              
              # -- Now apply the criteria that deal with the primer dimer problem.
              if(PRIMERDIMER == TRUE){
                forwardCuts <- forwardCuts * validForwardDeletionPRIMERDIMER
                reverseCuts <- reverseCuts * validReverseDeletionPRIMERDIMER
                cutsAfterPRIMERDIMER <- sum(forwardCuts)+sum(reverseCuts)
              }
              else{
                cutsAfterPRIMERDIMER <- cutsAfterSAMECUT
              }
        
              readAfterPRIMERDIMER       <- (cutsAfterPRIMERDIMER > 0)
        
              # Fill the deletion table with the cut info; now they can be a deletion or a cut
              if(length(forwardCuts) > 0){
              for(m in 1:length(forwardCuts)){
                type <- "deletion"
                if(forwardCuts[m] == 1){
                  type <- "cut"
                }
                  
                # TODO:
                # VERY VERY IMPORTANT
                # The Gotoh has a bug where it returns an interval of [start,end-1] instead of [start,ends] if ends is the
                # end of the amplicon. This piece of code correct that BUT THAT MUST BE CORRECTED IN GOTOH.CPP!!!!
                addjustedEndForward <- forwardRelativeDeletionsEnds[m]
                if(forwardRelativeDeletionsEnds[m] == ampliconLength-1){
                  addjustedEndForward <- addjustedEndForward + 1
                }
                 
                writeLines( paste(currentID,
                                  t,
                                  as.numeric(forwardRelativeDeletionsStarts[m]),
                                  as.numeric(addjustedEndForward),
                                  as.numeric(uniqueTable$Frequency[k]),
                                  type,
                                  "Forward",
                                  addjustedEndForward - forwardRelativeDeletionsStarts[m],
                                  sep = '\t') , cutsFD)
                  
              }}
        
              if(length(reverseCuts) > 0){
              for(m in 1:length(reverseCuts)){
                type <- "deletion"
                if(reverseCuts[m] == 1){
                  type <- "cut"
                }
            
                # TODO:
                # VERY VERY IMPORTANT
                # The Gotoh has a bug where it returns an interval of [start,end-1] instead of [start,ends] if ends is the
                # end of the amplicon. This piece of code correct that BUT THAT MUST BE CORRECTED IN GOTOH.CPP!!!!
                addjustedEndReverse <- reverseRelativeDeletionsEnds[m]
                if(reverseRelativeDeletionsEnds[m] == ampliconLength-1){
                  addjustedEndReverse <- addjustedEndReverse + 1
                }
            
                writeLines( paste(currentID,
                                  t,
                                  as.numeric(reverseRelativeDeletionsStarts[m]),
                                  as.numeric(addjustedEndReverse),
                                  as.numeric(-uniqueTable$Frequency[k]),
                                  type,
                                  "Reverse",
                                  addjustedEndReverse - reverseRelativeDeletionsStarts[m],
                                  sep = '\t') , cutsFD)
            
              }}
        
              # Get the cuts statistics after the last criteria is applied.
              forwardCuts <- sum(forwardCuts)
              reverseCuts <- sum(reverseCuts)
              cuts <- forwardCuts + reverseCuts
              total_cuts <- cuts * uniqueTable$Total[k]
        
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              applyCriterias <- applyCriterias + td        
        
            }
  
            #------------------------------------------------------
            # Fill the unique table with the results
            #------------------------------------------------------
            {
              t1 <- Sys.time()
  
              # Get the frameshift
              uniqueTable$Forward_Frameshift[k] <- forwardFrameshiftClass
              uniqueTable$Reverse_Frameshift[k] <- reverseFrameshiftClass
  
              uniqueTable$Forward_Cuts[k] <- forwardCuts
              uniqueTable$Reverse_Cuts[k] <- reverseCuts
  
              uniqueTable$Is_Read[k]        <- ((forwardCuts + reverseCuts)> 0)
  
              #-------------------INDEPENDENT-----------------------------------------
              uniqueTable$Unique_Cuts_Passed_ALLIN[k]       <- forwardCutsALLIN + reverseCutsALLIN
              uniqueTable$Unique_Cuts_Passed_SAMECUT[k]     <- cutsSAMECUT
              uniqueTable$Unique_Cuts_Passed_PRIMERDIMER[k] <- forwardCutsPRIMERDIMER + reverseCutsPRIMERDIMER
  
              #-------------------SUCCESSIVELY-----------------------------------------
              uniqueTable$Unique_Cuts_After_ALLIN[k]       <- cutsAfterALLIN
              uniqueTable$Unique_Cuts_After_SAMECUT[k]     <- cutsAfterSAMECUT
              uniqueTable$Unique_Cuts_After_PRIMERDIMER[k] <- cutsAfterPRIMERDIMER
              
              t2 <- Sys.time()
              td <- as.numeric(t2-t1, units = "secs")
              fillUniqueTable <- fillUniqueTable + td
              
            }
            
          } #For each of the unique rows
        
      }
      t6 <- Sys.time()
      td <- as.numeric(t6-t5, units = "secs")
      processUnique <- td

      #------------------------------------------------------
      # RECORD STATISTICS
      # Several statistics about deletions, cuts and reads 
      #------------------------------------------------------
      {
        t1 <- Sys.time()

        #print(uniqueTable)
        # Read Rate, what it is commonly call Cut Rate
        if(configTable$Sum_Is_Sequence[j]>0){

          configTable$Cut_Rate[j] <- sum(subset(uniqueTable, Is_Read==TRUE)$Total) / configTable$Sum_Is_Sequence[j]
          
          # Deletions, valid or not
          configTable$Forward_No_Deletions[j]       <- sum(subset(uniqueTable, Forward_Frameshift==1)$Total)
          configTable$Forward_In_Frameshift[j]      <- sum(subset(uniqueTable, Forward_Frameshift==2)$Total)
          configTable$Forward_Out_Frameshift[j]     <- sum(subset(uniqueTable, Forward_Frameshift==3)$Total)
          configTable$Forward_Multiple_Deletions[j] <- sum(subset(uniqueTable, Forward_Frameshift==4)$Total)
          
          configTable$Reverse_No_Deletions[j]       <- sum(subset(uniqueTable, Reverse_Frameshift==1)$Total)
          configTable$Reverse_In_Frameshift[j]      <- sum(subset(uniqueTable, Reverse_Frameshift==2)$Total)
          configTable$Reverse_Out_Frameshift[j]     <- sum(subset(uniqueTable, Reverse_Frameshift==3)$Total)
          configTable$Reverse_Multiple_Deletions[j] <- sum(subset(uniqueTable, Reverse_Frameshift==4)$Total)
          
          # Subset and get only those which are reads
          uniqueReadsTable <- subset(uniqueTable, Is_Read == TRUE)
          
          # If we have at least one read, find the frameshift for those
          if(nrow(uniqueReadsTable)>0){
          
            # Deletions, ONLY valid reads
            configTable$Forward_Read_No_Deletions[j]       <- sum(subset(uniqueReadsTable, Forward_Frameshift==1)$Total)
            configTable$Forward_Read_In_Frameshift[j]      <- sum(subset(uniqueReadsTable, Forward_Frameshift==2)$Total)
            configTable$Forward_Read_Out_Frameshift[j]     <- sum(subset(uniqueReadsTable, Forward_Frameshift==3)$Total)
            configTable$Forward_Read_Multiple_Deletions[j] <- sum(subset(uniqueReadsTable, Forward_Frameshift==4)$Total)
            
            configTable$Reverse_Read_No_Deletions[j]       <- sum(subset(uniqueReadsTable, Reverse_Frameshift==1)$Total)
            configTable$Reverse_Read_In_Frameshift[j]      <- sum(subset(uniqueReadsTable, Reverse_Frameshift==2)$Total)
            configTable$Reverse_Read_Out_Frameshift[j]     <- sum(subset(uniqueReadsTable, Reverse_Frameshift==3)$Total)
            configTable$Reverse_Read_Multiple_Deletions[j] <- sum(subset(uniqueReadsTable, Reverse_Frameshift==4)$Total)
            
          }
          
          # Only valid deletions
          configTable$Sum_Forward_Cuts[j] <- sum( uniqueTable$Forward_Cuts * uniqueTable$Total )
          configTable$Sum_Reverse_Cuts[j] <- sum( uniqueTable$Reverse_Cuts * uniqueTable$Total )
          
          # Only reads
          configTable$Sum_Is_Read[j] <- sum(subset(uniqueTable, Is_Read==TRUE)$Total)
          
          #-------------------INDEPENDENT-----------------------------------------
          configTable$Sum_Unique_Cuts_Passed_ALLIN[j]        <- sum( uniqueTable$Unique_Cuts_Passed_ALLIN * uniqueTable$Total ) # Independent cuts
          configTable$Sum_Unique_Cuts_Passed_SAMECUT[j]      <- sum( uniqueTable$Unique_Cuts_Passed_SAMECUT * uniqueTable$Total )
          configTable$Sum_Unique_Cuts_Passed_PRIMERDIMER[j]  <- sum( uniqueTable$Unique_Cuts_Passed_PRIMERDIMER * uniqueTable$Total )
          
          configTable$Sum_Unique_Reads_Passed_ALLIN[j]       <- sum(subset(uniqueTable, Unique_Cuts_Passed_ALLIN > 0)$Total) # Independent reads
          configTable$Sum_Unique_Reads_Passed_SAMECUT[j]     <- sum(subset(uniqueTable, Unique_Cuts_Passed_SAMECUT > 0)$Total)
          configTable$Sum_Unique_Reads_Passed_PRIMERDIMER[j] <- sum(subset(uniqueTable, Unique_Cuts_Passed_PRIMERDIMER > 0)$Total)
          
          #-------------------SUCCESSIVELY-----------------------------------------
          configTable$Sum_Unique_Cuts_After_ALLIN[j]        <- sum( uniqueTable$Unique_Cuts_After_ALLIN * uniqueTable$Total ) # Successively cuts
          configTable$Sum_Unique_Cuts_After_SAMECUT[j]      <- sum( uniqueTable$Unique_Cuts_After_SAMECUT * uniqueTable$Total )
          configTable$Sum_Unique_Cuts_After_PRIMERDIMER[j]  <- sum( uniqueTable$Unique_Cuts_After_PRIMERDIMER * uniqueTable$Total )
          
          configTable$Sum_Unique_Reads_After_ALLIN[j]       <- sum(subset(uniqueTable, Unique_Cuts_After_ALLIN > 0)$Total) # Successively reads
          configTable$Sum_Unique_Reads_After_SAMECUT[j]     <- sum(subset(uniqueTable, Unique_Cuts_After_SAMECUT > 0)$Total)
          configTable$Sum_Unique_Reads_After_PRIMERDIMER[j] <- sum(subset(uniqueTable, Unique_Cuts_After_PRIMERDIMER > 0)$Total)
          
          
        }
        else{

          configTable$Cut_Rate[j] <- NA
        
        }

        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        recordStats <- td        
        
        # Write the timings
        if(TIMING == TRUE){
          
          configTiming$Time_ReadBasic[j]            <- readBasicData
          configTiming$Time_ProcessRow[j]           <- processUnique
          configTiming$Time_RecordStat[j]           <- recordStats
          configTiming$Time_ExpandUniques[j]        <- expandUniques
          configTiming$Time_InitVariables[j]        <- initializeVariables
          configTiming$Time_ExtractAlignmentInfo[j] <- extractAlignmentInfo
          configTiming$Time_ALLIN[j]                <- allInTiming
          configTiming$Time_SAMECUT[j]              <- sameCutTiming
          configTiming$Time_PRIMERDIMER[j]          <- primerDimerTiming
          configTiming$Time_Frameshift[j]           <- frameShiftTiming
          configTiming$Time_ApplyCriterias[j]       <- applyCriterias
          configTiming$Time_FillUnique[j]           <- fillUniqueTable
          
#           configTable$Time_ProcessRow[j]    <- processUnique
#           
#           configTable$Time_ExpandUniques[j]        <- expandUniques
#           configTable$Time_InitVariables[j]        <- initializeVariables
#           configTable$Time_ExtractAlignmentInfo[j] <- extractAlignmentInfo
#           configTable$Time_ALLIN[j]                <- allInTiming
#           configTable$Time_SAMECUT[j]              <- sameCutTiming
#           configTable$Time_PRIMERDIMER[j]          <- primerDimerTiming
#           configTable$Time_Frameshift[j]           <- frameShiftTiming
#           configTable$Time_ApplyCriterias[j]       <- applyCriterias
#           configTable$Time_FillUnique[j]           <- fillUniqueTable
          

          
        }

        
      }
      
      # Write the unique table into disk
      write.table(uniqueTable, paste(currentIDCodeResultsFolderName,"/",currentID,"_",currentBarcode,
                                     "_uniques.txt",sep = ''), sep="\t")

      # Close the cuts table FDs
      close(cutsFD)

    }
    
    t4 <- Sys.time()
    td2 <- as.numeric(t4-t3, units = "secs")
    allTime <- td2
    
    configTiming$Total_Time[j] <- allTime
    
    # Now that the whole config is finish, write it to disk
    write.table(configTable, paste(analysisPath,"/",processID,"_configFile_results.txt",sep = ''), sep="\t")

    if(TIMING == TRUE){
      write.table(configTiming, paste(analysisPath,"/",processID,"_configFile_Timing",sep = ''), sep="\t")

    }
    
    
  }
  
  print(paste(processID, "finish!"))
  
  # Close the log file
  close(logFileConn)
  
}


  

