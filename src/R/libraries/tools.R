{
# This function return the reverse complementary dna string of a given dna
# string.
# 
# The function respect uppercases. Any other character that is not
# a,c,g,t,A,C,G,T is ignored.
# 
# The function takes the following parameters:
# 
#   dna: (String) A string with the nucleotide sequence.
# 
# The function returns:
#
#   (String) A string with the nucleotides reversed and complemented
# 
# Invariant:
#   
#   The returned string must have the same length as the input string.
}
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

{
  # This function takes a lite string and gives back the events coordinates
  #
  # For the given example:
  # 0         1         2         3         4         5
  # 012345678901234567890123456789012345678901234567890123456789
  # NNNNNNNNNNN-----NNNNNNN-------NNNNNNNN--------NNNNNNNN------    Subject
  # ------NNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNN---NNNNNNNNNN    Pattern
  # 
  # The RCPP library returns the information regarding the start and ends of the events in the following way:
  # (Please notice that there is a +1 shift in the index)
  #
  # @4@ 12,16 * 24,30 * 39,46 * 55,60 * ! @3@ 1,6 * 31,33  * 48,50 *
  #
  # The function returns:
  #     - A list with 6 elements
  #     - The element 1 and element 4 is the total of insertions and deletions (in that order)
  #     - The element 2 and element 5 is where the insertions and deletions starts (in that order)
  #     - The element 3 and element 6 is where the insertions and deletion ends in (that order)
  #     - Elements 2,3,4,5 are integer arrays.
  #     - If the total of insertions/deletions is 0, the integer arrays will be equal to NULL
}
getEventInfo <- function(liteString){
  
  if(is.na(liteString) == FALSE){
  
	  # Divide the string in two pieces; one for the insertions and the other one for the deletions
	  alignmentData <- unlist(strsplit(liteString, "!", fixed = TRUE))
	  alignmentInsertionsData <- alignmentData[1]
	  alignmentDeletionsData  <- alignmentData[2]
	  alignmentMissmatchData  <- alignmentData[3]
	  
	  # The totals of events are after the first @ and before the second @
	  alignmentInsertionsTotal  <- as.numeric(unlist(strsplit(alignmentInsertionsData,"@",fixed = TRUE))[2])
	  alignmentDeletionsTotal   <- as.numeric(unlist(strsplit(alignmentDeletionsData ,"@",fixed = TRUE))[2])
	  alignmentMissmatchesTotal <- as.numeric(unlist(strsplit(alignmentMissmatchData ,"@",fixed = TRUE))[2])
	  
	  # We need 7 arrays, 2 for the insertions, 2 for the deletions, 3 for the missmatches
	  #
	  # For deletions and insertions:
	  # Each of those two arrays, will represent the starts and the ends of each event.
	  # So for event 1, we need to look where does it start in the start_array[1] and where does ends in end_array[1]
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
	  alignmentInsertionsBoundariesData  <- unlist(strsplit(alignmentInsertionsData,"@",fixed = TRUE))[3]
	  alignmentInsertionsBoundariesArray <- unlist(strsplit(alignmentInsertionsBoundariesData,"*",fixed = TRUE))
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
	  alignmentDeletionsBoundariesData <- unlist(strsplit(alignmentDeletionsData,"@",fixed = TRUE))[3]
	  alignmentDeletionsBoundariesArray <- unlist(strsplit(alignmentDeletionsBoundariesData,"*",fixed = TRUE))
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
	  alignmentMissmatchBoundariesData <- unlist(strsplit(alignmentMissmatchData,"@",fixed = TRUE))[3]
	  alignmentMissmatchBoundariesArray <- unlist(strsplit(alignmentMissmatchBoundariesData,"*",fixed = TRUE))
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
  
  }
  
  else{

	  return (list(0, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL))
    
  }
  
}

{
  # For a given string, detect how many groups of uppercases are there, where are
  # they, and how long are they.
  # 
  # For example:
  #
  #    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
  #
  # Has 4 groups of uppercases, and they have length 7,4,1, and 3 respectevily.
  # The function returns a list with the length of each one. So the length of the
  # list represent.
  # The function returns a list with 3 elements:
  # -- How many groups we found. If will be 1 if we found nothing.
  # -- The starting position for each AP. It would be -1 if we found nothing.
  # -- The wide for each AP. It will be 1 if we found nothing.
  # 
  # The function takes the following parameters:
  # 
  #   candidate: (String) A string with the nucleotide sequence.
  # 
  # The function returns:
  #
  #   (List) A list with three element. The first element is an integer with the
  #          total of upper case group; lets say is X. The second element is a
  #          list of integer of size X with the start position of each group. The
  #          third element is a list of integer with the wide of each group.
  # 
}
countUppercaseGroups <- function(candidate){
  
  # Add a tiny lowercase letter at the end of candidate. That will trick gregexpr, and avoid making weird cases later.
  candidate <- paste(candidate,"a",sep='')
  
  alignmentPositionsStarts <- gregexpr("(?<![A-Z])[A-Z]",  candidate, ignore.case = FALSE, perl = TRUE, fixed = FALSE)[[1]]
  alignmentPositionsEnds   <- gregexpr("([A-Z])([a-z,-])", candidate, ignore.case = FALSE, perl = TRUE, fixed = FALSE)[[1]]
  
  total <- length(alignmentPositionsStarts)
  
  length <- (alignmentPositionsEnds - alignmentPositionsStarts) + 1
  
  return(list(total,alignmentPositionsStarts,length))
  
}

{
# This function takes a dataframe and a number of processors and gives back the
# same dataframe divided into several subsets. Every subset crafted so each
# processors takes aproximatly the same amount of rows. Note that that doesn't
# mean that they will take the same workload. This is just a arbitrary way of
# dividing the work.
#
# See:
#
#   divideWorkBySize()
#   divideWorkByBarcode()
# 
# The function takes the following parameters:
# 
#   totalProcessors: (Int)       An integer with the number of processors.
#                                This is also the amount of groups that will
#                                result from this function. The number must be
#                                1 or bigger.
# 
#   configDataframe: (Dataframe) A dataframe you want to divide. It doesn't
#                                require any special column. But it must not be
#                                empty. It can be full of NA rows thought.
# 
# The function returns:
#
#   (List) A list with several dataframes. Each dataframe is a subgroup of 
#          configDataframe. If the number of rows of configDataframe is bigger
#          than totalProcessors, the list have length totalProcessors. If it is
#          smaller, the length is equal to the number of rows of
#          configDataframe.
#
}
divideWork <- function(totalProcessors, configDataframe){

  # Make some initial calculations
  totalConfigLines <- nrow(configDataframe)
  processorsLines <- rep(as.integer(totalConfigLines/totalProcessors),totalProcessors)
  
  if((totalConfigLines%%totalProcessors)>0){
  for(i in 1:((totalConfigLines%%totalProcessors)) ){    
    processorsLines[i] <- processorsLines[i] + 1
  }}
  
  processorsSubDataframes <- list()
  
  start <- 1
  end <- 0
  offset <- 0

  for(i in 1:totalProcessors){

    # Find out the ending line
    end <- offset + processorsLines[i]
    offset <- end
    
#     print(paste("For processor",i))
#     print(start)
#     print(end)
    
    # Make a subset based on that
    subData <- NULL
    
    if(start <= end){
      subData <- configDataframe[start:end,]
    }
    
    processorsSubDataframes[[i]] <- subData
    
    # Update the next start to the current end + 1
    start <- end + 1
    
  }
  
  return(processorsSubDataframes)
  
}

{
# This function divide the dataframe into similar dataframes based on the
# colName values. The function try to divide the work so each processor gets
# aproximaly the same amount of colName summations. This is commonly known as
# the PARTITION PROBLEM. This implementation uses a greedy approach with order
# O(nlogn).
#
# See:
#
#   divideWork()
#   divideWorkByBarcode()
# 
# The function takes the following parameters:
# 
#   totalProcessors: (Int)       An integer with the number of processors.
#                                This is also the amount of groups that will
#                                result from this function. The number must be
#                                1 or bigger.
# 
#   configDataframe: (Dataframe) A dataframe you want to divide. It must not be
#                                empty.
#
#   colName:         (String)    The name of a column in configDataframe. This
#                                column should be full of of values that can be
#                                sorted and can be added (typically either int
#                                or float). It can be negative values.
# 
# The function returns:
#
#   (List) A list with several dataframes. Each dataframe is a subgroup of 
#          configDataframe. If the number of rows of configDataframe is bigger
#          than totalProcessors, the list have length totalProcessors. If it is
#          smaller, the length is equal to the number of rows of
#          configDataframe.
#
}
divideWorkBySize <- function(totalProcessors, configDataframe, colName){

  # In here we are going to store how many lines goes for each processor
  # at the beggining is a list of empty lists (slow)
  processorsLines <- rep( list(list()), totalProcessors ) 
  
  # In here we store the sum for each processor
  processorsSums <- rep(0, totalProcessors)
  
  # Get the data sorted by value in decreased order
  rowsSorted <- order(configDataframe[,c(colName)], decreasing = TRUE)
  
  # Get the total of lines
  totalConfigLines <- nrow(configDataframe)
  
  # Initialize the lists. Each one get the top rows
  for(i in 1:totalProcessors){
  
    processorsLines[i] <- c(rowsSorted[i])
    processorsSums[i]  <- configDataframe[,c(colName)][rowsSorted[i]]
    
  }
  
  # Now; for the rest of the non-top rows, added to set that will render the smaller sum
  if(totalConfigLines > 1){
  for(i in (totalProcessors+1):totalConfigLines){
    
    #Initialize the variables for this iteration
    analyzeLine  <- rowsSorted[i]
    currentValue <- configDataframe[,c(colName)][analyzeLine]
    
    # Minimum and reference
    minimum <- processorsSums[1] + currentValue
    reference <- 1
    
    # Try all combinations (except the first one which is minimum by default)
    if(totalProcessors>=2){
    for(j in 2:totalProcessors){
    
      if(( processorsSums[j] + currentValue) < minimum ){
        
        minimun <- processorsSums[j] + currentValue
        reference <- j
      
      }
    
    }}
    
    # Add the minimum to the list and update everything
    processorsLines[[reference]] <- c(processorsLines[[reference]],analyzeLine)
    processorsSums[reference] <- processorsSums[reference] + currentValue
    
  }}
  
  # Finally, we have all the lines divided. Create each dataframe and return them
  processorsSubDataframes <- list()
  
  for(i in 1:totalProcessors){
    
#     print(paste("For processor",i))
#     print(processorsLines[[i]])
#     print(length(processorsLines[[i]]))
#     print(processorsSums[i])
    
    processorsSubDataframes[[i]] <- configDataframe[processorsLines[[i]],]
    
  }
  
  return(processorsSubDataframes)

}

{
# This function takes in a config dataframe and divide the work into barcodes.
# This means that instead of dividing the work by the number of sequences per
# line as in divideWorkBySize, it will consider the barcode group as a whole,
# and give packages based on barcodes.
#
# That way, we don't need to ddply afterwards many times, only once. Each
# barcode have a unique combination of files from where it gets its reads.
#
# See:
#
#   divideWork()
#   divideWorkBySize()
# 
# The function takes the following parameters:
# 
#   totalProcessors: (Int)       An integer with the number of processors.
#                                This is also the amount of groups that will
#                                result from this function. The number must be
#                                1 or bigger.
# 
#   configDataframe: (Dataframe) A dataframe you want to divide. It must not be
#                                empty.
#
#   configFilePath:  (String)    The config file from where the configDataframe
#                                was generated. The dataframe contain the
#                                forward and reverse files paths. If it is not
#                                the absolute path, we need to find where are
#                                the files using this path.
# 
# The function returns:
#
#   (List) A list with several dataframes. Each dataframe is a subgroup of 
#          configDataframe. If the number of rows of configDataframe is bigger
#          than totalProcessors, the list have length totalProcessors. If it is
#          smaller, the length is equal to the number of rows of
#          configDataframe.
#
}
divideWorkByBarcode <- function(totalProcessors, configDataframe, configFilePath){

  # Get the barcodes
  uniqueBarcodeDF    <- configDataframe[c("Barcode")]
  uniqueBarcodeDF    <- unique(uniqueBarcodeDF)
  
  # Now lets sort the config and the uniques by barcodes
  sortedConfigDataframe <- configDataframe[order(configDataframe$Barcode),]
  uniqueBarcodeDF <- data.frame(uniqueBarcodeDF[order(uniqueBarcodeDF$Barcode),])
  colnames(uniqueBarcodeDF) <- c("Barcode")
  
  # ---------------------------------------------------------------
  # ADD THE SIZE OF THE FORWARD FILE TO EACH BARCODE
  # ---------------------------------------------------------------
  {
    filePath <- ""
    
    uniqueBarcodeDF$Size <- 0
    
    currentBarcode <- 1
    
    # We always have at least one barcode. So we can fill the first one already
    
    # Check that we are in an absolute path or relative path
    # (Not ./folder/folder/file.txt path, but /folfer/file.txt or file.txt , something not unixy)
    absolutePathForward <- ( length(grep("/", sortedConfigDataframe$"Forward_Reads_File"[1], ignore.case = TRUE)) >= 1 )
    
    # Get the names of the forward read and reverse read files
    if(absolutePathForward == TRUE){
      filePath <- as.character(sortedConfigDataframe$"Forward_Reads_File"[1])
    }
    else{
      dataFolderPath <- getConfigFolder(configFilePath)
      filePath <- paste(dataFolderPath, "/",sortedConfigDataframe$"Forward_Reads_File"[1],sep = '')
    }
    
    #filePath <- as.string(configDataframe$Forward_Reads_File[1])
    uniqueBarcodeDF$Size[1] <- round(file.info(filePath)[["size"]]/1048576,2)  
    
    # Now check all the barcodes in the config
    for(i in 1:nrow(sortedConfigDataframe)){
    
      if(as.character(sortedConfigDataframe$Barcode[i]) != as.character(uniqueBarcodeDF$Barcode[currentBarcode])){
    
        currentBarcode <- currentBarcode + 1
        
        # Get the names of the forward read and reverse read files
        if(absolutePathForward == TRUE){
          filePath <- as.character(sortedConfigDataframe$"Forward_Reads_File"[i])
        }
        else{
          dataFolderPath <- getConfigFolder(configFilePath)
          filePath <- paste(dataFolderPath, "/",sortedConfigDataframe$"Forward_Reads_File"[i],sep = '')
        }
        
        #filePath <- as.string(configDataframe$Forward_Reads_File[1])
        uniqueBarcodeDF$Size[currentBarcode] <- round(file.info(filePath)[["size"]]/1048576,2)  
      
      }
      
    }
  
  }
  
  # ---------------------------------------------------------------
  # DIVIDE THE BARCODES FOR EACH PROCESSOR
  # ---------------------------------------------------------------
  {

    # In here we are going to store how many barcodes goes for each processor
    # at the beggining is a list of empty lists (slow)
    processorsBarcodes <- rep( list(list()), totalProcessors ) 
    
    # In here we store the sum for each processor
    processorsSums <- rep(0.0, totalProcessors)
    
    # Get the data sorted by value in decreased order
    rowsSorted <- order(uniqueBarcodeDF[,c("Size")], decreasing = TRUE)
    
    # Get the total of barcodes
    totalBarcodes <- nrow(uniqueBarcodeDF)
    
    # Initialize the lists. Each one get the top rows
    for(i in 1:totalProcessors){
      
      processorsBarcodes[i] <- c(rowsSorted[i])
      processorsSums[i]  <- uniqueBarcodeDF$Size[rowsSorted[i]]
      
    }
    
    # Now; for the rest of the non-top rows, added to set that will render the smaller sum
    if(totalProcessors < totalBarcodes){
    for(i in (totalProcessors+1):totalBarcodes){
      
      #Initialize the variables for this iteration
      analyzeLine <- rowsSorted[i]
      currentValue <- uniqueBarcodeDF$Size[rowsSorted[i]]
      
      # Minimum and reference
      minimum <- processorsSums[1] + currentValue
      reference <- 1
      
      # Try all combinations (except the first one which is minimum by default)
      if(totalProcessors>=2){
        for(j in 2:totalProcessors){
          
          if(( processorsSums[j] + currentValue) < minimum ){
            
            minimun <- processorsSums[j] + currentValue
            reference <- j
            
          }
          
        }}
      
      # Add the minimum to the list and update everything
      processorsBarcodes[[reference]] <- c(processorsBarcodes[[reference]],analyzeLine)
      processorsSums[reference] <- processorsSums[reference] + currentValue
      
    }}
    
  }

  # ---------------------------------------------------------------
  # CREATE THE DATAFRAME FOR EACH PROCESSOR
  # ---------------------------------------------------------------
  {
  # Now we have the Barcode ID for each processor, we just need to divide the original config dataframe
  # base upon their barcode
    
  processorsSubDataframes <- rep( list(c()), totalProcessors ) 
  
  currentCandidate <- 1
  
  
  # For each of the groups in the barcodes group
  for(i in 1:length(processorsBarcodes)){
    
    # For each of the member of that group, which are Barcodes IDs
    for(j in 1:length(processorsBarcodes[[i]])){
      
      # Get the Barcode from the list of barcodes
      candidateBarcode <- as.character(uniqueBarcodeDF$Barcode[processorsBarcodes[[i]][j]])
      
      # Find out all the lines from the original config file that are suppose to come in here
      for(k in 1:nrow(configDataframe)){
      
        # If the Barcodes match
        if(as.character(configDataframe$Barcode[k]) == candidateBarcode){
  
          processorsSubDataframes[[i]] <- c(processorsSubDataframes[[i]],k)
          
        }
      
      }
  
    }
    
  }
  
  # We have discover which rows should go into each dataframe. Now, convert those integers into actuall dataframe
  for(i in 1:length(processorsSubDataframes)){
  
    print(paste("For processor",i))
#     print(as.character(uniqueBarcodeDF[processorsBarcodes[[i]],]$Barcode))
    print("This many barcodes")
    print(length(processorsBarcodes[[i]]))
    print("Total forward reads filesize of (MB)")
    print(processorsSums[i])
    
    processorsSubDataframes[[i]] <- configDataframe[processorsSubDataframes[[i]],]
    
  }
  
  
  }
  
  
  return (processorsSubDataframes)

}

{
# This function takes all the files in a given folder that match a regex
# pattern and put them toguether into a single file.
#
# BE CAREFULL: This function __WILL_DELETE__ all the original little files
# afterward if one of the parameters it set to TRUE.
#
# The function takes the following parameters:
# 
#   targetFolder:  (String) The folder where the files are suppose to be.
#
#   regex:         (String) The files to be combined, must comply with this
#                           regex.
#
#   finalFileName: (String) The name of the final file.
#
#   header:        (Bool)   If header is true, it will take the header of the
#                           first file only and ignore the rest. Otherwise it
#                           will repeat the header many times.
#
#   delete:        (Bool)   If true, it will delete the original files
#
#   isrecursive:   (Bool)   If recursive is true, it will look inside all the
#                           folder of the given folder for more files which 
#                           name match the regex.
# 
# The function returns:
#
#   (String) A string with the nucleotides reversed and complemented
# 
# Invariant:
#   
#   The returned string must have the same length as the input string.
}
unifyFiles <- function (targetFolder, regex, finalFileName, header, delete, isrecursive){
  
  # Get all the files
  allFiles <- list.files(targetFolder, recursive = isrecursive)
  
  # Now get those only with the regex expresion
  candidates <- grep(regex, allFiles, fixed=FALSE, ignore.case=FALSE)  
  
  if(length(candidates)>0){
    
    # Inside this folder it will be a file named as the user desired
    finalFilePath <- paste(targetFolder, "/", finalFileName, sep = '')
    finalFileFD   <- file(finalFilePath, open="w")
    
    for(i in 1:length(candidates)){

      candidateName <- paste(targetFolder,"/",allFiles[candidates[i]],sep="")
      
      # Read the file into a string variable
      textFile <- readLines(candidateName,encoding="UTF-8")  
      
#       print("From:")
#       print(candidateName)
#       print("Get:")
#       print(textFile)
#       
      if(header == TRUE && i!=1){
        textFile <- textFile[2:length(textFile)]
      }
      
      # Write the text into the unified folder
      writeLines(textFile, finalFileFD)
      
      if(delete == TRUE){
        # Remove the old file
        file.remove(candidateName)
      }
      
    }

    # Close the file descriptor for the final file
    close(finalFileFD)

  }

  return (length(candidates))

}

{
# This funcion takes a txt file, with a file path per row, and delete all those files.
#
# The function takes the following parameters:
#
#    listOfFilesFileName: (string) Path to the file that containts 0 or more files paths that will be deleted
#
#    deleteSource:        (bool)   If you want that the file listOfFilesFileName to be ALSO deleted, set this
#                                  to TRUE. Default is FALSE
# The function returns:
#  
#    (int) How many files where deleted. This could be less than the orinal files, as it will only delete
#          files that you have write privilege.
}
deleteFiles <- function (listOfFilesFileName, deleteSource = FALSE){

  # Keep track of how many files we delete
  totalDeleted <- 0

  # Read all the files candidates
  textFile <- readLines(toString(listOfFilesFileName),encoding="UTF-8")
  
  # For each of the candidates
  for (k in 1:length(textFile)){
  
    # Check if we have write permissions for that file
    access <- checkFileWriteAccess(textFile[k])
    
    # If we do, delete it.
    if(access == TRUE){
      totalDeleted <- totalDeleted + 1
      file.remove(textFile[k])
    }
    
  }
  
  return (totalDeleted)

}

{
# This function check if the given file exist and can be read.
# 
# The function takes the following parameters:
# 
#   filePath: (String) A string the path to the file.
# 
# The function returns:
#
#   (Bool) TRUE if we have file access. This means:
#          -- If the file path is not NULL.
#          -- If the file path is not NA.
#          -- If the file path is not "Null", "null", "na", ...
#          -- If the file path have read permissions.
#                 
#          FALSE otherwise.
# 
}
checkFileAccess <- function (filePath){
  
  result <- FALSE
  
  # Check that the config file exist and we have read access
  if(!is.null(filePath[1])){ #First we need to check if it is null or we will get a warning
    
    if(!is.na(filePath[1])){ #Same with NA
      
      # Check that we don't have weird NULL strings
      if(nchar(filePath[1]) != 0    && filePath[1] !="NULL" && filePath[1]!="null" &&
               filePath[1]  != "na" && filePath[1] != "NA"){
        
        # Check that we have read access, NOTE!: If we do, the function returns 0, and -1 if we don't.
        if(file.access(filePath[1], mode = 4) == 0 ){
          
          result <- TRUE
          
        }
      }  
    }
  }
  
  return (result)
  
}

{
  # This function check if the given file exist and can be write.
  # 
  # The function takes the following parameters:
  # 
  #   filePath: (String) A string the path to the file.
  # 
  # The function returns:
  #
  #   (Bool) TRUE if we have file write access. This means:
  #          -- If the file path is not NULL.
  #          -- If the file path is not NA.
  #          -- If the file path is not "Null", "null", "na", ...
  #          -- If the file path have write permissions.
  #                 
  #          FALSE otherwise.
  # 
}
checkFileWriteAccess <- function (filePath){
  
  result <- FALSE
  
  # Check that the config file exist and we have read access
  if(!is.null(filePath[1])){ #First we need to check if it is null or we will get a warning
    
    if(!is.na(filePath[1])){ #Same with NA
      
      # Check that we don't have weird NULL strings
      if(nchar(filePath[1]) != 0    && filePath[1] !="NULL" && filePath[1]!="null" &&
           filePath[1]  != "na" && filePath[1] != "NA"){
        
        # Check that we have write access, NOTE!: If we do, the function returns 0, and -1 if we don't.
        if(file.access(filePath[1], mode = 2) == 0 ){
          
          result <- TRUE
          
        }
      }  
    }
  }
  
  return (result)
  
}

{
  # This function get a reads file and put into a dataframe.
  # 
  # The files are usually stored in a .gz file. If the .gz file has been unzip
  # before the function won't try to unzip it, it will read the .fastq instead.
  # If the file is a .gz, and there is no .fastq in the same folder with the same
  # name, the function will unzip the .gz file and leave the .fastq in the same
  # folder.
  #
  # Prerrequisites:
  #
  #   The file must exist, and you must have reads rights.
  #
  #   If the file is zipped, you must have write rights for the writing folder.
  #
  # The function takes the following parameters:
  #   
  #   fileName:   (String) A String with the absolute path to the file.
  #
  #   TEMPFOLDER: (String) Optional, a String with the absolute path to a folder.
  #                        In this folder is where the files are going to be unziped.
  #                        If nothing is specify, the files will be unziped in the same
  #                        folder where the zipped files are.
  #
  #   tempFileConn: (File) Optional, a file where we write the path of the temporal file.
  #                        If a file is unzipped, we might want to consider it a temoporal file
  #                        and delete it later. In here we keep track of it.
  # 
  # The function returns:
  #
  #   If there is no errors:
  #   (Dataframe) A dataframe with the following columns:
  # 
  #       [01] Sequence line
  #       [02] Quality
  #  
  #   If the FASTQ file has a number of files that is not divisible by 4:
  #   NULL
  # 
  # Invariant:
  #   
  #   The unzipped file contain a plain text file. This file is formatted with a
  #   FASTQ format, meaning that is repeating the following pattern over an over
  #   again:
  #   - Line with the ID and other metadata
  #   - Line with the sequence [01]
  #   - '+ line'
  #   - Line with the quality of the nucleotides [02]
  # 
  #   This means that the number of lines in that file must be divisible by 4.
  # 
  #   Worth mentioning that this function is suppose to be use for the forward
  #   and reverse reads. The number of lines in each of those FASTQ files must
  #   be the same.
  
}
getReadsFile <- function(fileName, TEMPFOLDER = NULL, tempFileConn = NULL, cropLeft = 0, cropRight = 0) {
  
  unzipReadsFileName <- NULL
  table.df           <- NULL
  
  # Check that the fileName has .gz in it.
  # If it does we need to unzip it. Otherwise we can start reading it directly.
  
  # If it is compressed
  if(length(grep(".gz", fileName, ignore.case = TRUE))==1){

    # Get the path of the original path without the .gz
    pathUnzipName <- sub(".gz", "", fileName, fixed=TRUE)
    
    # Get the name of the original file without the gz
    fileUnzipName <- basename(pathUnzipName)
    
    # Create the path of the final file
    unzipReadsFileName <- pathUnzipName # this is the path in the same folder as the .gz
    
    if(!is.null(TEMPFOLDER)){
      unzipReadsFileName <- paste(TEMPFOLDER,"/",fileUnzipName, sep='')
    }
    
    # Unzip the file
    gunzip(fileName, destname = unzipReadsFileName, skip = TRUE, overwrite = FALSE, remove = FALSE) #If it is already unzipped, use that
#     print("final name")
#     print(unzipReadsFileName)
    # If the users wants, keep track of the filepath for the unzipped file.
    if(!is.null(tempFileConn)){
    
      writeLines(unzipReadsFileName , tempFileConn)
      
    }
    
  }
  # If it is not compressed
  else{
    
    unzipReadsFileName <- fileName
    
  }
  
  # Read the file into a string variable
  textFile <- readLines(toString(unzipReadsFileName),encoding="UTF-8")
  
  if(length(textFile) %% 4 == 0){
    
    totalRowsTable = length(textFile)/4 # This MUST be always multiple of 4
    
    # Creates the matrix where the forward info is going to be stored
    # It has the following columns:
    # [01] Instrument      |
    # [02] Run             |
    # [03] Flowcell ID     |  This is the ID for that row
    # [04] Flowcell Line   |
    # [05] Tile            |
    # [06] X               |
    # [07] Y               |
    # [08] Member
    # [09] Filtered
    # [10] Control
    # [11] Index
    # [12] Sequence line
    # [13] + line
    # [14] Quality
    matrixTable <- matrix(ncol=2, nrow=totalRowsTable)
    
    # Fill the matrix with forward reads  
    z <- 0
    #pBar <- txtProgressBar(style=3)
    for (k in 1:length(textFile)){
      
      # Move forward the progress bar
      #setTxtProgressBar(pBar, k/length(textFile))
      
      # Check if we are in the METADATA, SEQUENCE, + , QUALITY line
      z <- k %% 4
      
      # Sequence line
      if(z==2){
        matrixTable[ceiling(k/4),1] <-  substr(textFile[k], cropLeft+1, nchar(textFile[k])-cropRight)
      }
      
      # Quality line
      if(z==0){
        matrixTable[ceiling(k/4),2] <- substr(textFile[k], cropLeft+1, nchar(textFile[k])-cropRight)
      }
      
    }
    
    # Transform the matrix into a dataframe
    table.df <- data.frame(matrixTable)  
    colnames(table.df) <- c("Sequence","Quality") 
    
  }
  else{
    
    print(fileName)
    print("INCORRECT FILE: This file has a number of lines that is not multiple of 4!")
    
  }
  
  return (table.df)
  
}

{
  # Given a config file path, get the folder containing that config file
  #
  # /home/myUser/myGenesData/experiment.txt
  #
  # The function returns:
  #
  # myGenesData
  #
  # The function takes the following paramenters:
  #   string configFilePath: A string with the absolute path pointing to the config file
  #
  # The function returns the following variables:
  #   
  #   string :  A string with the name of the folder containing the config file
  #
  # Invariant:
  #
  #   The file path must be Unix style
}
getConfigFolder <- function(configFilePath){
  
  # Separate the name by points.
  auxiliaryNames <- unlist(strsplit(configFilePath, ".", fixed = TRUE))
  
  # The path can have an IP address, or not. In any case, we need to add all the list toguether again
  # and skip that last member, while putting a dot in between each member.
  finalPath = auxiliaryNames[1]
  
  if(length(auxiliaryNames) > 2){
    for (i in 2:(length(auxiliaryNames)-1)){
      
      finalPath = paste(finalPath,".",auxiliaryNames[i])
      
    }}
  
  # In here; we have the complete path including the config file name minus the extension
  # We just need to take out the name of the config file away
  configName <- getConfigName(configFilePath)
  
  #   print(configFilePath)
  #   print(configName)
  #   print(finalPath)
  
  
  finalPath <- substr(finalPath, 1, (nchar(finalPath) - nchar(configName)) )  
  
  #return (finalPath) Isn't the following much more efficient?
  
  # Divide by /
  foldersNames <- unlist(strsplit(configFilePath, "/", fixed = TRUE))
  
  # Put toguether everything again except the last element of the list (which is the name of the config file)
  finalPath <- paste(foldersNames[1:length(foldersNames)-1], collapse = '/')
  
  # We need to add a / at the start of the finalPath
  #finalPath <- paste("/",finalPath,sep='') #NOPE! Because the first / in the config path was splitted, and give us a "" (empty string) at position [1], so collapsing with / will give us the initial /
  
  return (finalPath)
  
}

{
  # Given a config file path, get the name of the config file. For example, for the path:
  #
  # /home/myUser/myGenesData/experiment.txt
  #
  # The function returns:
  #
  # experiment.txt
  #
  # The function takes the following paramenters:
  #   string configFilePath: A string with the absolute path pointing to the config file
  #
  # The function returns the following variables:
  #   
  #   string :  A string with the name of the config file
  #
  # Invariant:
  #
  #   The file path must be Unix style
}
getConfigName <- function(configFilePath){
  
  # Separate the name by points.
  auxiliaryNames <- unlist(strsplit(configFilePath, ".", fixed = TRUE))
  
  # The path can have an IP address, or not. In any case, the second name from the end of that list
  # contain the name of the parent folder. Inside that name, we split by the folder separator "/"
  auxiliaryFolders <- unlist(strsplit(auxiliaryNames[length(auxiliaryNames) - 1], "/", fixed = TRUE))
  
  # The last element is the name of the config file
  configName <- auxiliaryFolders[length(auxiliaryFolders)]
  
  return (configName)
  
}

{
  # Gives a config file path and check if exist, is readable, etc...
  #
  # The function takes the following paramenters:
  #   string configFilePath: A string with the absolute path pointing to the config file
  #
  # The function returns the following variables:
  #   
  #   bool , TRUE:  Everything went good.
  #          FALSE: If the config file doesn't exist or if we don't have read access
  #                 If the string provided was 'NULL', 'NA', '0', 'null' or 'na'
  #
  # Invariant:
  #
  #   The file path must be Unix style
}
checkReady <- function(configFilePath){
  
  result <- FALSE
  
  # Check that the config file exist and we have read access
  if(!is.null(configFilePath[1])){ #First we need to check if it is null or we will get a warning
    
    if(!is.na(configFilePath[1])){ #Same with NA
      
      # Check that we don't have weird NULL strings
      if(nchar(configFilePath[1]) != 0 && configFilePath[1]!="NULL" && configFilePath[1]!="null" &&
           configFilePath[1] != "na" && configFilePath[1] != "NA"){
        
        # Check that we have read access, NOTE!: If we do, the function returns 0, and -1 if we don't.
        if(file.access(configFilePath[1], mode = 4) == 0 ){
          
          result <- TRUE
          
        }
      }  
    }
  }
  
  return (result)
  
}

checkFolderAcess <- function(folderPath){}
