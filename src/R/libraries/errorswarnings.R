{
# This function check if the target RNA is in the amplicon.
# 
# The function takes the following paramenters:
#   
#   string targetPrimer: A sequence of nucleotides in a string format representing the target
#          amplicon:     A sequence of nucleotides in a string format representing the amplicon
#          ID:           The ID from where this target and amplicon came.
#          barcode:      The barcode from where this target and amplicon came.
#          configFilePath: The location on disk where you can find the config file with this ID and barcode.
# 
#   file descriptor logFileConn: A file descriptor in R which is pointing to the log file.
# 
# The function returns the following variables:
#   
#   bool , TRUE:  Everything went good.
#          FALSE: If the test fails and a warning is annotated.
# 
# Invariant:
#   
#   The file descriptor must be set to append.
}
checkTarget <- function (targetPrimer, amplicon, ID, barcode, configFilePath, logFileConn){

  toReturn <- TRUE
  
  targetPositions <- regexpr(targetPrimer, amplicon, ignore.case=TRUE, fixed=FALSE)[1]
  if(targetPositions <= 0 ){
    
    toReturn <- FALSE
    print ("WARNING!: Target not found in the amplicon. This event will be register in the log")
    writeLines("Couldn't find the target primer:", logFileConn)
    writeLines(toString(targetPrimer), logFileConn)
    writeLines("In amplicon:", logFileConn)
    writeLines(toString(amplicon), logFileConn)
    writeLines("For ID and barcode:", logFileConn)
    writeLines(toString(paste(ID, barcode, sep=" ")), logFileConn)
    writeLines("Which bellongs to the config file located in:", logFileConn)
    writeLines(toString(configFilePath), logFileConn)
    writeLines("\n", logFileConn)
    
  }
  
  return (toReturn)
  
}

{
# This function check if the forward anr reverse primer is in the amplicon.
#
# The forward primer should be the one coming from the configuration file, while the reverse primer should be
# reverse and complement. This function DOES NOT reverse complement the input. So if you look up for something
# that is suppose to containt a reverse complement it won't find it. However you can use it to find any pair of
# primers regarless if they are reverse complemented or not.
# 
# The function takes the following paramenters:
#   
#   string forwardPrimer: A sequence of nucleotides in a string format representing the forward primer
#          reversePrimer: A sequence of nucleotides in a string format representing the reverse primer. Usually
#                         you want to reverse complement this BEFORE giving it to the function. The function WILL NOT
#                         reverse complement it for you. (see: /libraries/tools.R , you have there the reverse
#                         complement function)
#          amplicon:     A sequence of nucleotides in a string format representing the amplicon
#          ID:           The ID from where this target and amplicon came.
#          barcode:      The barcode from where this target and amplicon came.
#          configFilePath: The location on disk where you can find the config file with this ID and barcode.
# 
#   file descriptor logFileConn: A file descriptor in R which is pointing to the log file.
# 
# The function returns the following variables:
#   
#   bool , TRUE:  Everything went good.
#          FALSE: If the test fails and a warning is annotated.
# 
# Invariant:
#   
#   The file descriptor must be set to append.
}
checkPrimers <- function (forwardPrimer, reversePrimerReverseComplementary, amplicon, ID, barcode,
                          configFilePath, logFileConn){

  toReturn <- TRUE
  
  forwardPrimerPosition <- regexpr(forwardPrimer, amplicon, ignore.case=TRUE, fixed=FALSE)[1]
  reversePrimerPosition <- regexpr(reversePrimerReverseComplementary, amplicon, ignore.case=TRUE, fixed=FALSE)[1]
  if(forwardPrimerPosition <= 0 || reversePrimerPosition <=0){
    toReturn <- FALSE
    print ("WARNING!: One of primer was not found in the amplicon. This event will be register in the log")
    writeLines("Couldn't find the forward primer or reverse primer reversed and complemented):", logFileConn)
    writeLines(toString(forwardPrimer), logFileConn)
    writeLines(toString(reversePrimerReverseComplementary), logFileConn)
    writeLines("In amplicon:", logFileConn)
    writeLines(toString(amplicon), logFileConn)
    writeLines("For ID and barcode:", logFileConn)
    writeLines(toString(paste(ID, barcode, sep=" ")), logFileConn)
    writeLines("Which bellongs to the config file located in:", logFileConn)
    writeLines(toString(configFilePath), logFileConn)
    writeLines("\n", logFileConn)
  }
  
  return (toReturn)
}

{
# This function check if the given alignment positions are valid
#
# The function takes the following paramenters:
#   
#   array<int> alignmentPositions: an array of integers representing the aligning positions inside the amplicon.
#
#   string   amplicon:     A sequence of nucleotides in a string format representing the amplicon
#            ID:           The ID from where this target and amplicon came.
#            barcode:      The barcode from where this target and amplicon came.
#            configFilePath: The location on disk where you can find the config file with this ID and barcode.
# 
#   file descriptor logFileConn: A file descriptor in R which is pointing to the log file.
# 
# The function returns the following variables:
#   
#   bool , TRUE:  Everything went good.
#          FALSE: If the test fails and a warning is annotated.
# 
# Invariant:
#   
#   The file descriptor must be set to append.
#   The alignmentPositions array can also be a single value. For R, every single variable is an array of size one.
#   The alignmentPositions must be 0 or less to indicate that is not found. The first character in the amplicon is
#   the character one, not zero.
}
checkPositions <- function (alignmentPositions, amplicon, ID, barcode, configFilePath, logFileConn){

  toReturn <- TRUE
  
  if(alignmentPositions[1] <= 0){
    warningTrigger <- FALSE
    print ("WARNING!: Aligment position not found in the amplicon. This event will be register in the log")
    writeLines("Couldn't find alignment position for amplicon:", logFileConn)
    writeLines(toString(amplicon), logFileConn)
    writeLines("For ID and barcode:", logFileConn)
    writeLines(toString(paste(ID, barcode, sep=" ")), logFileConn)
    writeLines("Which bellongs to the config file located in:", logFileConn)
    writeLines(toString(configFilePath), logFileConn)
    writeLines("\n", logFileConn)
  }
  
  return (toReturn)
  
}

{
# This function pre-process a config file and check that everything is in order
# 
# Its takes care of the following:
#   No IDs are duplicated
#   Every combination of barcode, forward primer and reverse primer is unique.
#   Each barcode has unique forward reads file and reverse read files.
#   Check that the read files exist and with reading access
# 
# If anything is wrong it will annotated in a log file.
# 
# The function takes the following paramenters:
# The function returns the following variables:
#   
#    bool , TRUE:  Everything went good.
#           FALSE: If the test fails and a warning is annotated.  
#   
# Invariant:
#
#    The dataframe must have the following columns:
#      "ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Experiment_Type", "Target_Primer",
#      "Forward_Primer","Reverse_Primer","Strand","Genome".
#    The config table is not empty.
}  
checkConfigFile <- function (configTable, configFilePath, logFileConn){

  
  
  toReturn <- TRUE
  
  # Get the total of rows in the config table
  totalRows <- nrow(configTable)
  
  #-------------------------------------
  #   NO NAs VALUES IN THE CONFIG FILE
  #-------------------------------------
  {
    goodRows      <- complete.cases(configTable)
    totalGoodRows <- sum(goodRows)
    
    if(totalGoodRows != nrow(configTable)){
    
      toReturn <- FALSE
      print ("WARNING!: Skipping config file, some rows have NA/NULL data!")    
      writeLines("This config file has bad rows due to NA/NULL values:", logFileConn)
      writeLines(toString(configFilePath), logFileConn)
      writeLines("The rows with bad data are:", logFileConn)
      writeLines( toString(which(goodRows==FALSE)), logFileConn)
      
    }
    
    
  }
  
  #----------------------------------
  #   NO IDS ARE DUPLICATED
  #----------------------------------
  {
    # Find out all the uniques IDs
    IDSubset <- subset(configTable, select = c("ID"))
    IDSubset <- ddply(IDSubset,.(IDSubset$"ID"),nrow)
    
    # If all IDs are unique, then the number of rows haven't change
    # If now we have less, we have duplicate IDs
    if(totalRows != nrow(IDSubset)){
    
      toReturn <- FALSE
      print ("WARNING!: Skipping config file, not all IDs are unique. This event will be register in the log")    
      writeLines("This config file has duplicates IDs:", logFileConn)
      writeLines(toString(configFilePath), logFileConn)
      
    }
  }
  
  #-----------------------------------------------------
  #   BARCODE, FORWARD, AND REVERSE PRIMERS ARE UNIQUE
  #-----------------------------------------------------  
  {
#     # Find the unique combinations of barcodes, forwards, and reverse primers.
#     barcodeAndPrimersSubset <- subset(configTable, select = c("Barcode","Forward_Primer","Reverse_Primer"))
#     barcodeAndPrimersSubset <- ddply(barcodeAndPrimersSubset,.(barcodeAndPrimersSubset$"Barcode",
#                                                                barcodeAndPrimersSubset$"Forward_Primer",
#                                                                barcodeAndPrimersSubset$"Reverse_Primer"),nrow)
#     colnames(barcodeAndPrimersSubset) <- c("Barcode","FWD","REV","Total")
#     
#     # If the sum of the totals doesn't match number of rows then we have duplicates.
#     if(sum(barcodeAndPrimersSubset$"Total") != nrow(barcodeAndPrimersSubset)){
#       
#       toReturn <- FALSE #TODO: Delete this part?
#       print ("WARNING!: Skipping config file, not all combination of Barcode, Forward primer and reverse primer is unique. This event will be register in the log")
#       writeLines("This config file has duplicates barcode + forward + reverse:", logFileConn)
#       writeLines(toString(configFilePath), logFileConn)
#       writeLines("These are wrong:", logFileConn)
#       for(i in 1:nrow(barcodeAndPrimersSubset)){
#         
#         if(barcodeAndPrimersSubset$"Total"[i] > 1){
#           writeLines(paste(toString(barcodeAndPrimersSubset$"Barcode"[i]),
#                            toString(barcodeAndPrimersSubset$"FWD"[i]),
#                            toString(barcodeAndPrimersSubset$"REV"[i]),
#                            toString(barcodeAndPrimersSubset$"Total"[i]),sep=" "), logFileConn)
#           
#         }
#         
#       }
#       
#     }  
  }
  
  #----------------------------------------------------------
  #   EACH BARCODE HAS UNIQUE FORWARD AND REVERSE READ FILES
  #----------------------------------------------------------     
  {
  # Find the unique combinations of barcodes, forwards and reverse files, and forward and reverse primers
  barcodeAndFilesSubset <- subset(configTable, select = c("Barcode","Forward_Reads_File","Reverse_Reads_File"))
  barcodeAndFilesSubset <- ddply(barcodeAndFilesSubset,.(barcodeAndFilesSubset$"Barcode",
                                                         barcodeAndFilesSubset$"Forward_Reads_File",
                                                         barcodeAndFilesSubset$"Reverse_Reads_File"
  ),nrow)
  
  colnames(barcodeAndFilesSubset) <- c("Barcode","FWD","REV", "Total")
  # In here we have the list of barcodes, forward and reverse reads which are uniques. We can have many of each
  # and that is not a problem yet.
  
  # Now we are going to get the unique barcodes.  
  barcodeTable           <- ddply(barcodeAndFilesSubset,.(barcodeAndFilesSubset$"Barcode"),nrow)
  colnames(barcodeTable) <- c("Barcode","Total")

  badCandidates <- subset(barcodeTable, Total > 1, select=c("Barcode", "Total"))
  
#   print(barcodeTable)
#   print(badCandidates)
#   print(barcodeAndFilesSubset)
  
  # We can have several barcodes, repeat on the config with the same files, thats is OK. But not with different
  # files. If that happen, we will have more rows in the barcodeAndFilesSubset than in the barcodeTable
  if(nrow(barcodeTable) != nrow(barcodeAndFilesSubset)){
    toReturn <- FALSE
    print ("WARNING!: Skipping config file, There is a non unique combination of barcode, forward file and reverse file.")
  
    writeLines("This config file has duplicates barcode + forward read file + reverse read file:", logFileConn)
    writeLines(toString(configFilePath), logFileConn)
    writeLines("These barcodes have a non unique read files", logFileConn)
    
    for(i in 1:nrow(badCandidates)){
      
      badRows <- subset(barcodeAndFilesSubset, Barcode == badCandidates$Barcode[i], select=c("Barcode","FWD","REV"))
    
      print(badRows)
      print(toString(badRows$Barcode[1]))
      
      for(j in 1:nrow(badRows)){
        
        writeLines(
          paste(
            toString(badRows$Barcode[j]),
            toString(badRows$FWD[j]),
            toString(badRows$REV[j]),
            sep=" "), logFileConn)
      }
      
    }
    
  }

  
  }
#   TODO: Fix this
#   #----------------------------------------------------------
#   #   EACH BARCODE HAS UNIQUE FORWARD AND REVERSE READ FILES, FORWARD AND REVERSE PRIMERS
#   #----------------------------------------------------------     
#   {
#     # Find the unique combinations of barcodes, forwards and reverse files, and forward and reverse primers
#     barcodeAndFilesSubset <- subset(configTable, select = c("Barcode","Forward_Reads_File","Reverse_Reads_File","Forward_Primer","Reverse_Primer"))
#     barcodeAndFilesSubset <- ddply(barcodeAndFilesSubset,.(barcodeAndFilesSubset$"Barcode",
#                                                            barcodeAndFilesSubset$"Forward_Reads_File",
#                                                            barcodeAndFilesSubset$"Reverse_Reads_File",
#                                                            barcodeAndFilesSubset$"Forward_Primer",
#                                                            barcodeAndFilesSubset$"Reverse_Primer"
#                                                            ),nrow)
#     
#     colnames(barcodeAndFilesSubset) <- c("Barcode","FWD","REV", "FWD_P", "RVS_P", "Total")
#     
#     # In here we have the list of barcodes, forward and reverse reads which are uniques. We can have many of each
#     # and that is not a problem yet.
#     # Now; the forward column must not have duplicates, and the same is true for the reverse column. If that happen
#     # we have he same file for different barcodes.
#     
#     # We count how many repetions we have en each row. If one of them is greater than one, then we have a problem
#     badCandidates <- subset(barcodeAndFilesSubset, Total > 1, select=c("Barcode", "Total"))
# 
#     if(nrow(badCandidates) > 1){
#     
#       toReturn <- FALSE
#       print ("WARNING!: Skipping config file, There is a non unique combination of barcode, forward file, reverse file, forward primer and reverse primer. This event will be register in the log")
#       writeLines("This config file has duplicates barcode + forward read file + reverse read file + forward primer + reverse primer:", logFileConn)
#       writeLines(toString(configFilePath), logFileConn)
#       writeLines("One of these is not like the others:", logFileConn)
#       for(i in 1:nrow(barcodeAndFilesSubset)){
#             
#         writeLines(paste(toString(barcodeAndFilesSubset$"Barcode"[i]),
#                          toString(barcodeAndFilesSubset$"FWD"[i]),
#                          toString(barcodeAndFilesSubset$"REV"[i]),
#                          toString(barcodeAndFilesSubset$"FWD_P"[i]),
#                          toString(barcodeAndFilesSubset$"RVS_P"[i]),
#                          toString(barcodeAndFilesSubset$"Total"[i]),sep=" "), logFileConn)
#               
#       }
#             
#     }
# 
#   }

  #----------------------------------------------------------
  #   WE HAVE ACCESS TO ALL READ FILES
  #---------------------------------------------------------- 
  { 
  
    uniquesForwards <- subset(barcodeAndFilesSubset, select = c("FWD"))
    uniquesReverses <- subset(barcodeAndFilesSubset, select = c("REV"))
        
    uniquesForwards <- ddply(uniquesForwards,.(uniquesForwards$"FWD"),nrow)
    uniquesReverses <- ddply(uniquesReverses,.(uniquesReverses$"REV"),nrow)
        
    colnames(uniquesForwards) <- c("FWD","Total")
    colnames(uniquesReverses) <- c("REV","Total")
    
    # Check the forward files
    for(i in 1:nrow(uniquesForwards)){
  
      # For each name, check if we have an absolute or relative path. Is a relative path if we have any /
      # If we have a relative path we need to find the reads files inside the same folder of the config file
      # Otherwise, just take the path given in the string
      
      # Check that we are in an absolute path
      absolutePath <- ( length(grep("/", uniquesForwards$"FWD"[i], ignore.case = TRUE)) >= 1 )
      
      
      if(absolutePath == TRUE){
        
        forwardReadsFilePath = uniquesForwards$"FWD"[i]
        
      }
      else{
        
        dataFolderPath <- getConfigFolder(configFilePath)
        forwardReadsFilePath = paste(dataFolderPath, "/",uniquesForwards$"FWD"[i],sep = '')
        
      }
      
      
      # Check that we have read access, NOTE!: If we do, the function returns 0, and -1 if we don't.
      if(file.access(toString(forwardReadsFilePath), mode = 4) == -1 ){
      
        toReturn <- FALSE
        print ("WARNING!: Skipping config file. Can't read file. This event will be register in the log")
        print (toString(forwardReadsFilePath))
        writeLines("This config file has a read file with no read access (maybe doesn't exist?):", logFileConn)
        writeLines(toString(configFilePath), logFileConn)
        writeLines("The file path is:", logFileConn)
        writeLines(toString(forwardReadsFilePath), logFileConn)
        
      }
    
    }
    # Check the reverse files
    for(i in 1:nrow(uniquesReverses)){
      
      
      
      # For each name, check if we have an absolute or relative path. Is a relative path if we have any /
      # If we have a relative path we need to find the reads files inside the same folder of the config file
      # Otherwise, just take the path given in the string
      
      # Check that we are in an absolute path
      absolutePath <- ( length(grep("/", uniquesForwards$"REV"[i], ignore.case = TRUE)) >= 1 )
      
      
      if(absolutePath == TRUE){
        
        reverseReadsFilePath = uniquesForwards$"REV"[i]
        
      }
      else{
        
        dataFolderPath <- getConfigFolder(configFilePath)
        reverseReadsFilePath = paste(dataFolderPath, "/",uniquesForwards$"REV"[i],sep = '')
        
      }
      
      # Check that we have read access, NOTE!: If we do, the function returns 0, and -1 if we don't.
      if(file.access(toString(reverseReadsFilePath), mode = 4) == -1 ){
        
        toReturn <- FALSE
        print ("WARNING!: Skipping config file. Can't read file. This event will be register in the log")
        print (toString(forwardReadsFilePath))
        writeLines("This config file has a read file with no read access (maybe doesn't exist?):", logFileConn)
        writeLines(toString(configFilePath), logFileConn)
        writeLines("The file path is:", logFileConn)
        writeLines(toString(reverseReadsFilePath), logFileConn)
        
      }
      
    }
  }
  
  return (toReturn)

}
  
{
# This function write info into a file that is human readable. It also tells
# the users via console that there was an error, and where he can find further
# verbose information.
#
# The function takes the following paramenters:
#
#     (string) currentConfigPath: The place where is the original config file.
#     (string) logFileName:       The path where the log file is placed.
#                                 Usually, the name goes like this:
#                                     resultsFolder + "/MasterLog.txt"
#     (file)   logFileConn:       The file object that refers to "logFileName".
#
# The function returns the following variables:
#     (void)
}
writeConfigError <- function(currentConfigPath, logFileName, logFileConn){

  print("-----------------------------------------")
  print("---|             ERROR!              |---")
  print("-----------------------------------------")
  print("ERROR!: Invalid config file.")
  print(currentConfigPath)
  print("   Please check the following: ")
  print("   -- No ID is duplicated")
  print("   -- Each combination of barcode, forward and reverse primer have different reads files")
  print("   -- All read files exist and you have read permission")
  print("   -- No NA/NULL values")
  print(" ")
  print("All errors and warnings are recorded in the log")
  print("Please check the Master Log for more detailed information about what caused this error")
  print(paste("Check '", logFileName, "' for detailed information"))
  
  writeLines("ERROR: Invalid config file: ", logFileConn)
  writeLines(currentConfigPath , logFileConn)
  writeLines("\n" , logFileConn)

}  
