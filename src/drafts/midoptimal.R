#Final script

decipher_installed <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

parallel_installed <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

utils_installed    <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

biostring_installed<- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

ggplot2_installed  <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

reshape2_installed <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

dplyr_installed    <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

total_processors <- 1       # Set this to the number of processors you want to use.
                            # You can set it to NULL if you want all that are available.

skip_N_Nucleotides <- TRUE  # Some sequences have faulty nucleotides labels with N.
                            # If we find a sequence like that in either forwards or reverse,
                            # we skip that aligment

average_Quality <- 0        # The FASTQ file have a quality for each nucleotide, being ! the lower and ~ the highest
                            # In ASCII :
                            # !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
                            # This quality number goes from 0 to 100, being 0 the lowest quality and 100 the highest
                            # You can write whatever number in between and the program will find out the apropiate
                            # character encoding. The filter work by converting each character to a number, and then
                            # find the average. If the average fall above this threshold then we take the sequence.
  
min_Quality <- 0            # Same as before, but this is the minimum quality for ALL nucleotide. If one of them
                            # has quality BELLOW this threshold, then the sequence is skipped.

writeAlignments <- FALSE     # Write the aligments results into disk 

writeRAW <- FALSE            # Write the aligments indels and mismatches into disk

############################################
#Install the packages if necessary
############################################
{
if(!decipher_installed){
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("DECIPHER")
  
}
if(!parallel_installed){
  
  install.packages("doParallel")
  
}
if(!utils_installed){
  
  install.packages("R.utils")
}
if(!biostring_installed){
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  
}
if(!ggplot2_installed){
  
  install.packages("ggplot2")
}
if(!reshape2_installed){
  
  install.packages("reshape2")
}
if(!dplyr_installed){
  install.packages("dplyr")
}

}
############################################
#Summons the libraries and get to use functions from there.
############################################
{
library(DECIPHER)
library(doParallel)
library(Biostrings)
library(R.utils) #Unzipping the .gz files
library(reshape2) #Rearrenge data for ggplot
library(ggplot2)
library(dplyr) #Makes code look nicer; implemented hash for duplicated sequences with this.
cl <- makeCluster(total_processors, outfile="") #outfile tells the workers to print on the terminal
registerDoParallel(cl)
}
############################################
#Algorith starts here
############################################


#Look into a folder call data, everything inside there, that is a folder, is interesting data (ie: Cas9_data)
#Everything else that is not a folder is just ignored.
parentFolderFiles <- list.files(getwd())
dataFolder <- grep("^data$", parentFolderFiles, fixed=FALSE, ignore.case=FALSE)
if(length(dataFolder)!=1){
  stop("No data folder found!")
}
allFilesNames <- list.files(paste(getwd(), "/data",sep = ''))
foldersIndexes  <- grep("^[^.]+$", allFilesNames, fixed=FALSE) #Find strings that DONT have a dot

#Gives an error if no directories are found
if(length(foldersIndexes)<0){
  stop("/data , No folders found!")
}

#Initialize array with the folders name. These are the folders that we will analyse
allFoldersNames <- rep(NA, length(foldersIndexes))

#Get the folders names from the files names array
for (i in 1:length(foldersIndexes)){
  allFoldersNames[i] <- allFilesNames[foldersIndexes[i]]
}

print(paste("Found a total of",length(allFoldersNames),"folders"))

#Check if there is a folder call /result in cwd, if not create a new one
dir.create(file.path(getwd(), "/results"), showWarnings = FALSE)

#Make a folder call /run_YYYYmmddHHMMSS where we are going to write the results for this run.
currentTime <- Sys.time()
timeStamp <-  strftime(currentTime,"%Y%m%d%H%M%S")
resultsFolder <- paste(getwd(), "/results/run_", timeStamp, sep = '')
dir.create(file.path(resultsFolder), showWarnings = TRUE)

#TODO: Filter the config files and take those you are interested.

#For each of the data folder/config files (there is only one config file per folder, ie: Cas9_data = run11.txt)
#--------------------------------------------------------------------------------------------------------------
for (i in 1:length(allFoldersNames)){
  
  #Find out the folder in where we are now
  currentDataFolderName <- paste(getwd(), "/data/",allFoldersNames[i],sep = '')
  print(paste("Working on",currentDataFolderName))
  
  #Make a folder in the results directory for this one
  currentResultsFolderName <- paste(resultsFolder,"/",allFoldersNames[i],"_results",sep = '')
  dir.create(file.path(currentResultsFolderName), showWarnings = TRUE)
  
  #Get all the files in that data folder
  folderFiles <- list.files(currentDataFolderName)
  
  #Get the config file name
  configFileName  <- folderFiles[grep(".conf", folderFiles, fixed=TRUE)] #Find a strings that contain .conf
  
  #If there is not config file, or if thre is more than one halt.
  if(length(configFileName)!=1){
    stop("NO CONFIG FOUND IN THIS FOLDER!!!")
  }
  
  configFileName  <- paste(currentDataFolderName,"/",configFileName,sep = '')
  
  #Read the config file into a dataset
  print("Creating config file dataframe")
  configTable <- read.table(configFileName, header=FALSE)
  colnames(configTable) <- c("ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Experiment_Type","Target_Primer","Forward_Primer","Reverse_Primer","Strand","Genome")
  
  #Prepare the Config_file_i_results.txt matrix
  #That have the following columns [ID,Barcode,FR_file,RR_file,F_Insertions,F_Deletions,F_Mismatches,R_Insertions,R_Deletions,R_Mismatches]
  matrixConfigFileResults<- matrix(ncol=10, nrow=nrow(configTable))

  #Creates a matrix where all the timings will go
  #It has these columns [concept, accumulated_time, times_annotated, 100 vector]
  #Recopilate the following concepts
  # [01] - Printing: This is when the user is presented with some sensible info.
  # [02] - Create_Datastructure: When we create empty dataset or matrices.
  # [03] - Unzipping_File: Extract the txt file from the gz file.
  # [04] - Reading_FASTQ_File: After unzipping, put the FASTQ txt into a string variable.
  # [05] - Fill_Datastructure: Transform that string into a datastructure with some meaning.
  # [06] - Transform_into_Dataframe: A matrix get transformed into a dataframe.
  # [07] - Get_Info: Put into a variable some relevant information.
  # [08] - Making_folders: Create an appropiate folder tree to write results inside and such.
  # [09] - Recycle_data: Uses some old datastructure to annotate some new stuff so we don't need to allocate memory again.
  # [10] - Initialize_variables: Reset variables to the default state, typically in a loop.
  # [11] - Search_for_primers: Time taken to looks for the forward and reverse primers, also look for the target.
  # [12] - Pairwise_FandR: Pairwise the forward + genome, also pairwise reverse + genome
  # [13] - Write Pairwise: The results of the pairwise goes into disk
  # [14] - Write RAW: Indels and mismatch information goes into disk
  # [15] - Write STATS: Writes statistics about the RAW into a datastructure ???
  # [16] - Config_line_complete: How long does it take to process a config line.
  # [17] - Write_variables: Get some info and modify the variables value (only in N-nucleotides, throw away?)
    
  # The way this works is like this:
  # We use the accumulated time to measure the total time and get the average. We also get the first 100 samples
  # for each group in orther to get a nice aproximation of the percentiles.
  
  # When the line is finish we get the average, sd, and quartile information for each concept, and then
  # we add it to the absolute matrix with rbind. That would take around 0.5s for each config line.
  # At the end we have a lot of info in the absolute matrix from where we take the data to draw the boxplots
  
  # The absolute matrix have these columns: [concept, average, sd, 0Q, 1Q, 2Q, 3Q, 4Q(MAX), Absolute_Frequency]
  matrixTimingAbsolute <- NULL

  #For each line in the config row
  for (j in 1:nrow(configTable)){

    #Initialize the timing variables
    matrixTiming <- matrix(0, ncol=103, nrow=17)
    t1 <- Sys.time()
    t2 <- Sys.time()
    t3 <- Sys.time()
    t4 <- Sys.time()
    td <- difftime(t2,t1)
    td2 <- difftime(t4,t3)
    
    matrixTiming[1,1]  <- "Printing"
    matrixTiming[2,1]  <- "Create_Datastructure"
    matrixTiming[3,1]  <- "Unzipping_File"
    matrixTiming[4,1]  <- "Reading_FASTQ_File"
    matrixTiming[5,1]  <- "Fill_Datastructure"
    matrixTiming[6,1]  <- "Transform_into_Dataframe"
    matrixTiming[7,1]  <- "Get_info"
    matrixTiming[8,1]  <- "Making_folders"
    matrixTiming[9,1]  <- "Recycle_data"
    matrixTiming[10,1] <- "Initialize_variables"
    matrixTiming[11,1] <- "Search_for_primers"
    matrixTiming[12,1] <- "Pairwise_FandR"
    matrixTiming[13,1] <- "Write_Pairwise"
    matrixTiming[14,1] <- "Write_RAW"
    matrixTiming[15,1] <- "Write_STATS"
    matrixTiming[16,1] <- "Config_Line_Complete"
    matrixTiming[17,1] <- "Writing_variables"
    
    t3 <- Sys.time()
    
    {
    
    #Get the ID and barcode
    t1 <- Sys.time()
    currentID <- configTable[j,"ID"]
    currentBarcode <- configTable[j,"Barcode"]
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Get info
    matrixTiming[7,2] <- as.numeric(matrixTiming[7,2]) + td
    matrixTiming[7,3] <- as.numeric(matrixTiming[7,3]) + 1
    if(as.numeric(matrixTiming[7,3]<101)){
      matrixTiming[7,as.numeric(matrixTiming[7,3])+3] <- td    
    }
        
    #Progress feedback for the user
    t1 <- Sys.time()
    print("________________________________________")
    print(paste("Working on",currentID, currentBarcode))
    print(paste("Line",j,"of",nrow(configTable),round(j/nrow(configTable)*100,2),"% for",currentDataFolderName))
    print("----------------------------------------")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3]<101)){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }
    
    #Make a folder for the results of the forwards/reverse reads that are going to come from here
    #The folder is inside currentResultsFolderName  (ie: run_??/Cas9_data_results/) and will be named
    #after the ID and Barcode for these forwards and reverse reads (ie:)
    # run_??/Cas9_data_results/phf8_1_1
    # run_??/Cas9_data_results/phf8_2_2
    # run_??/Cas9_data_results/phf8_3_3 ...
    t1 <- Sys.time()
    currentIDCodeResultsFolderName <- paste(currentResultsFolderName,"/",currentID,"_",currentBarcode,sep = '')
    dir.create(file.path(currentIDCodeResultsFolderName), showWarnings = TRUE)
        
    #Inside that folder, it will be another folder where we write the RAW alignmnets results
    currentAlignmentsResultsFolderName <- paste(currentIDCodeResultsFolderName,"/alignments",sep='')
    dir.create(file.path(currentAlignmentsResultsFolderName), showWarnings = TRUE)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Making folders
    matrixTiming[8,2] <- as.numeric(matrixTiming[8,2]) + td
    matrixTiming[8,3] <- as.numeric(matrixTiming[8,3]) + 1
    if(as.numeric(matrixTiming[8,3])<101){
      matrixTiming[8,as.numeric(matrixTiming[8,3])+3] <- td
    }
    
    #TODO: Filter by the line options, if it is interesting continue, otherwise skip to the next one
    #TODO: Check for the reverse amplicon, What was that again?
    
    #Get the names of the forward read and reverse read files
    t1 <- Sys.time()
    forwardReadsFileName = paste(currentDataFolderName, "/",configTable[j,"Forward_Reads_File"],sep = '')
    reverseReadsFileName = paste(currentDataFolderName, "/",configTable[j,"Reverse_Reads_File"],sep = '')

    #Get the target, forward, and reverse primer
    forwardPrimer <- configTable[j,"Forward_Primer"]
    reversePrimer <- configTable[j,"Reverse_Primer"]
    targetPrimer  <- configTable[j,"Target_Primer"]
      
    #Make the complementary reverse of the reverse primer
    myReverseComplementPrimer <- reverseComplement(DNAStringSet(reversePrimer))
      
    #Get the genome sequence
    genome <- configTable[j,"Genome"]
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Get info
    matrixTiming[7,2] <- as.numeric(matrixTiming[7,2]) + td
    matrixTiming[7,3] <- as.numeric(matrixTiming[7,3]) + 1
    if(as.numeric(matrixTiming[7,3])<101){
      matrixTiming[7,as.numeric(matrixTiming[7,3])+3] <- td
    }
        
    #-----------------------------------------------------------------
    #Get the actual FORWARD read file, unzip it, read it, an put it into a dataframe.
    #-----------------------------------------------------------------
    {
    t1 <- Sys.time()  
    print(paste("Unzipping",forwardReadsFileName))
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }
    
    #Unzip the file
    t1 <- Sys.time()
    gunzip(forwardReadsFileName, overwrite = TRUE, remove = FALSE)
    unzipForwardReadsFileName <- sub(".gz", "", forwardReadsFileName, fixed=TRUE)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Uzipping
    matrixTiming[3,2] <- as.numeric(matrixTiming[3,2]) + td
    matrixTiming[3,3] <- as.numeric(matrixTiming[3,3]) + 1
    if(as.numeric(matrixTiming[3,3])<101){
      matrixTiming[3,as.numeric(matrixTiming[3,3])+3] <- td
      
    }
    
    #Print feedback of the reading for the user, very large files can take a while to read
    t1 <- Sys.time()
    print(paste("Reading file of size",round(file.info(unzipForwardReadsFileName)[["size"]]/1048576,2),"MiBs"))
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }
    
    #Read the file into a string variable
    t1 <- Sys.time()
    textForwardFile <- readLines(unzipForwardReadsFileName,encoding="UTF-8")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Read FASTQ file
    matrixTiming[4,2] <- as.numeric(matrixTiming[4,2]) + td
    matrixTiming[4,3] <- as.numeric(matrixTiming[4,3]) + 1
    if(as.numeric(matrixTiming[4,3])<101){
      matrixTiming[4,as.numeric(matrixTiming[4,3])+3] <- td
    }
    
    #The number of lines must be multiple of 4 ALWAYS!
    if(length(textForwardFile) %% 4 != 0){
      print(forwardReadsFileName)
      stop("INCORRECT FORWARD FILE: The forward file has a number of lines that is not multiple of 4!")
    }
    
    t1 <- Sys.time()
    totalRowsForwardsTable = length(textForwardFile)/4 #This MUST be always multiple of 4
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Get info
    matrixTiming[7,2] <- as.numeric(matrixTiming[7,2]) + td
    matrixTiming[7,3] <- as.numeric(matrixTiming[7,3]) + 1
    if(as.numeric(matrixTiming[7,3])<101){
      matrixTiming[7,as.numeric(matrixTiming[7,3])+3] <- td
      
    }
    
    #Creates the matrix where the forward info is going to be stored    
    t1 <- Sys.time()
    matrixForwardsTable <- matrix(ncol=14, nrow=totalRowsForwardsTable)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Create_datastructure
    matrixTiming[2,2] <- as.numeric(matrixTiming[2,2]) + td
    matrixTiming[2,3] <- as.numeric(matrixTiming[2,3]) + 1
    if(as.numeric(matrixTiming[2,3])<101){
      matrixTiming[2,as.numeric(matrixTiming[2,3])+3] <- td
    }
    
    t1 <- Sys.time()
    print("Filling matrix for forward reads")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
      
    }
    
    #Fill the matrix with forward reads
    t1 <- Sys.time()
    z <- 0
    pBar <- txtProgressBar(style=3)
    for (k in 1:length(textForwardFile)){
      
      setTxtProgressBar(pBar, k/length(textForwardFile))
      z <- k %% 4
        
      #ID line
      if(z==1){
        
        #Separate fist by space, the left side is the Instrument until Y-coordinate (7 values), the right side
        #is from member until index (4 values); everything separated by ":"
        leftRight <- unlist(strsplit(textForwardFile[k], " ", fixed = TRUE))
              
        leftSide <- unlist(strsplit(leftRight[1], ":",fixed = TRUE))
        rightSide <- unlist(strsplit(leftRight[2], ":",fixed = TRUE))
              
        instrument <- leftSide[1]
        run <- leftSide[2]
        flowcellID <- leftSide[3]
        flowcellLine <- leftSide[4]
        tile <- leftSide[5]
        x <- leftSide[6]
        y <- leftSide[7]
        member <- rightSide[1]
        filtered <- rightSide[2]
        control <- rightSide[3]
        index <- rightSide[4]
  
        matrixForwardsTable[ceiling(k/4),1] <- instrument
        matrixForwardsTable[ceiling(k/4),2] <- run
        matrixForwardsTable[ceiling(k/4),3] <- flowcellID
        matrixForwardsTable[ceiling(k/4),4] <- flowcellLine
        matrixForwardsTable[ceiling(k/4),5] <- tile
        matrixForwardsTable[ceiling(k/4),6] <- x
        matrixForwardsTable[ceiling(k/4),7] <- y
        matrixForwardsTable[ceiling(k/4),8] <- member
        matrixForwardsTable[ceiling(k/4),9] <- filtered
        matrixForwardsTable[ceiling(k/4),10] <- control
        matrixForwardsTable[ceiling(k/4),11] <- index

      }
          
      #Sequence line
      if(z==2){
        matrixForwardsTable[ceiling(k/4),12] <- textForwardFile[k]
      }
        
      #+ line
      if(z==3){
        matrixForwardsTable[ceiling(k/4),13] <- textForwardFile[k]
      }
      
      #Quality line
      if(z==0){
        matrixForwardsTable[ceiling(k/4),14] <- textForwardFile[k]
      }
        
    }
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Fill_datastructure
    matrixTiming[5,2] <- as.numeric(matrixTiming[5,2]) + td
    matrixTiming[5,3] <- as.numeric(matrixTiming[5,3]) + 1
    if(as.numeric(matrixTiming[5,3])<101){
      matrixTiming[5,as.numeric(matrixTiming[5,3])+3] <- td
      
    }
    
    #Transform the matrix into a dataframe
    #TODO: Check by making the dataframe from the begging
    t1 <- Sys.time()
    forwardsTable <- data.frame(matrixForwardsTable)
    colnames(forwardsTable) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Member","Filtered","Control","Index","Sequence","Extra","Quality")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Transform_into_datafame
    matrixTiming[6,2] <- as.numeric(matrixTiming[6,2]) + td
    matrixTiming[6,3] <- as.numeric(matrixTiming[6,3]) + 1
    if(as.numeric(matrixTiming[6,3])<101){
      matrixTiming[6,as.numeric(matrixTiming[6,3])+3] <- td
      
    }
    #-------- END FORWARDS --------------
    
    }

    #-----------------------------------------------------------------
    #Get the actual REVERSE read file, unzip it, read it, an put it into a dataframe.
    #-----------------------------------------------------------------
    {
    t1 <- Sys.time()
    print(paste("Unzipping",reverseReadsFileName))
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
      
    }
    
    t1 <- Sys.time()
    gunzip(reverseReadsFileName, overwrite = TRUE, remove = FALSE)
    unzipReverseReadsFileName <- sub(".gz", "", reverseReadsFileName, fixed=TRUE)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Uzipping
    matrixTiming[3,2] <- as.numeric(matrixTiming[3,2]) + td
    matrixTiming[3,3] <- as.numeric(matrixTiming[3,3]) + 1
    if(as.numeric(matrixTiming[3,3])<101){
      matrixTiming[3,as.numeric(matrixTiming[3,3])+3] <- td
      
    }
    
    t1 <- Sys.time()
    print(paste("Reading file of size",round(file.info(unzipReverseReadsFileName)[["size"]]/1048576,2),"MiBs"))
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
      
    }
    
    t1 <- Sys.time()
    textReverseFile <- readLines(unzipReverseReadsFileName,encoding="UTF-8")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Read FASTQ file
    matrixTiming[4,2] <- as.numeric(matrixTiming[4,2]) + td
    matrixTiming[4,3] <- as.numeric(matrixTiming[4,3]) + 1
    if(as.numeric(matrixTiming[4,3])<101){
      matrixTiming[4,as.numeric(matrixTiming[4,3])+3] <- td
    }

    #The number of lines must be multiple of 4 ALWAYS!
    if(length(textReverseFile) %% 4 != 0){
      print(reverseReadsFileName)
      stop("INCORRECT FORWARD FILE: The forward file has a number of lines that is not multiple of 4!")
    }

    t1 <- Sys.time()
    totalRowsReversesTable = length(textReverseFile)/4 #This MUST be always multiple of 4
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Get info
    matrixTiming[7,2] <- as.numeric(matrixTiming[7,2]) + td
    matrixTiming[7,3] <- as.numeric(matrixTiming[7,3]) + 1
    if(as.numeric(matrixTiming[7,3])<101){
      matrixTiming[7,as.numeric(matrixTiming[7,3])+3] <- td
      
    }
#     print(td)
#     print(matrixTiming[7,2])
#     print(matrixTiming[7,3])
#     print(matrixTiming[7,4])
#     print(matrixTiming[7,5])
#     print(matrixTiming[7,6])
#     n <- readline("Stop4: ")
    
    t1 <- Sys.time()
    print("Creating empty matrix for reverse reads")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }
    
    t1 <- Sys.time()
    matrixReversesTable <- matrix(ncol=14, nrow=totalRowsReversesTable)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Create_datastructure
    matrixTiming[2,2] <- as.numeric(matrixTiming[2,2]) + td
    matrixTiming[2,3] <- as.numeric(matrixTiming[2,3]) + 1
    if(as.numeric(matrixTiming[2,3])<101){
      matrixTiming[2,as.numeric(matrixTiming[2,3])+3] <- td
    }
    
    t1 <- Sys.time()
    print("Filling matrix for reverse reads")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }
    
    
    t1 <- Sys.time()
    z <- 0
    pBar <- txtProgressBar(style=3)
    for (k in 1:length(textReverseFile)){
  
      setTxtProgressBar(pBar, k/length(textReverseFile))
      z <- k %% 4
  
      #ID line
      if(z==1){            
        
        #Separate fist by space, the left side is the Instrument until Y-coordinate (7 values), the right side
        #is from member until index (4 values); everything separated by ":"
        leftRight <- unlist(strsplit(textReverseFile[k], " ", fixed = TRUE))
          
        leftSide <- unlist(strsplit(leftRight[1], ":",fixed = TRUE))
        rightSide <- unlist(strsplit(leftRight[2], ":",fixed = TRUE))
          
        instrument <- leftSide[1]
        run <- leftSide[2]
        flowcellID <- leftSide[3]
        flowcellLine <- leftSide[4]
        tile <- leftSide[5]
        x <- leftSide[6]
        y <- leftSide[7]
        member <- rightSide[1]
        filtered <- rightSide[2]
        control <- rightSide[3]
        index <- rightSide[4]
          
        matrixReversesTable[ceiling(k/4),1] <- instrument
        matrixReversesTable[ceiling(k/4),2] <- run
        matrixReversesTable[ceiling(k/4),3] <- flowcellID
        matrixReversesTable[ceiling(k/4),4] <- flowcellLine
        matrixReversesTable[ceiling(k/4),5] <- tile
        matrixReversesTable[ceiling(k/4),6] <- x
        matrixReversesTable[ceiling(k/4),7] <- y
        matrixReversesTable[ceiling(k/4),8] <- member
        matrixReversesTable[ceiling(k/4),9] <- filtered
        matrixReversesTable[ceiling(k/4),10] <- control
        matrixReversesTable[ceiling(k/4),11] <- index
    
      }
      
      #Sequence line
      if(z==2){
        matrixReversesTable[ceiling(k/4),12] <- textReverseFile[k]
      }
        
      #+ line
      if(z==3){
        matrixReversesTable[ceiling(k/4),13] <- textReverseFile[k]
      }
        
      #Quality line
      if(z==0){
        matrixReversesTable[ceiling(k/4),14] <- textReverseFile[k]
      }
  
    }
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Fill_datastructure
    matrixTiming[5,2] <- as.numeric(matrixTiming[5,2]) + td
    matrixTiming[5,3] <- as.numeric(matrixTiming[5,3]) + 1
    if(as.numeric(matrixTiming[5,3])<101){
      matrixTiming[5,as.numeric(matrixTiming[5,3])+3] <- td
    }
    
    t1 <- Sys.time()
    reversesTable <- data.frame(matrixReversesTable)
    colnames(reversesTable) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Member","Filtered","Control","Index","Sequence","Extra","Quality")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Transform_into_datafame
    matrixTiming[6,2] <- as.numeric(matrixTiming[6,2]) + td
    matrixTiming[6,3] <- as.numeric(matrixTiming[6,3]) + 1
    if(as.numeric(matrixTiming[6,3])<101){
      matrixTiming[6,as.numeric(matrixTiming[6,3])+3] <- td
    }
    }
      
    #Make the matrix for the forward, and reverse alignments statistics.
    #The matrix would have this columns
    #Instrument ID, ... , X, Y, Forward Found, Reverse Found, Target Found, Number of insertions, Number of deletions, Number of mistchmatches, Alignment txt file
    #We can use the matrices from before since they are the same dimension, and half of the data is already there
    t1 <- Sys.time()
    matrixStatsForwardsTable <- matrixForwardsTable
    matrixStatsReversesTable <- matrixReversesTable
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Recycling_data
    matrixTiming[9,2] <- as.numeric(matrixTiming[9,2]) + td
    matrixTiming[9,3] <- as.numeric(matrixTiming[9,3]) + 1
    if(as.numeric(matrixTiming[9,3])<101){
      matrixTiming[9,as.numeric(matrixTiming[9,3])+3] <- td
    }

    {
    #Now we are going to make alignments and write the results into disk.
    #Note that the Alignment txt file of row i in both forward and reverse must point to the same file!
    #The folder structure is like this:
    #Data_folder_1 (ie: Cas9_data)
    # ...
    #Data_folder_X
    #Data_folder_1_results (ie: Cas9_data_results)
    #----Config_file_1_results (ie: run11_results)
    #-------- Results 1 (ie: 1_S1_L001_R1_001,1_S1_L001_R2_001)
    #------------ Alignments
    #---------------- Alignment_1.txt
    #---------------- ...
    #---------------- Alignment_Z.txt
    #------------ StatsForwardsTable.txt [Instrument ID,..., X, Y, # IN, # DEL, # MISS, Alignment txt file]
    #------------ StatsReverseTable.txt  [Instrument ID,..., X, Y, # IN, # DEL, # MISS, Alignment txt file]
    #-------- Results Y
    #-------- Config_file_1_results.txt (ie: run11_results.txt) [ID,Barcode,FR_file,RR_file,Insertions,Deletions,Mismatches]
    }    
      
    
    
    #Print feedback of the status of the program for the user
    t1 <- Sys.time()
    print("Looking for candidates to align")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Printing
    matrixTiming[1,2] <- as.numeric(matrixTiming[1,2]) + td
    matrixTiming[1,3] <- as.numeric(matrixTiming[1,3]) + 1
    if(as.numeric(matrixTiming[1,3])<101){
      matrixTiming[1,as.numeric(matrixTiming[1,3])+3] <- td
    }

    #Make a new dataframe with unique sequences, for forwards and reverses
    #This dataframe will have the columns [Sequence, Absolute_Frequency, Result_file]
    #Result_file is set to NA at the beggining.
    #The fist sequence that arrive there will find NA, do the alignment and write the result somewhere in disk
    #That will give us the disk path file and we set it into the Result_file column.
    #The next equal sequence that arrives, will not find NA, and just skip the alignment and copypaste the path
    #On its own summary row.
    #The unique sequences statistics goes into each config line folder later on

    #Group the forward and reverse tables by sequence, add nrows, and a blank string field

    #For each of the candidates
    pBar <- txtProgressBar(style=3)
    for (k in 1:nrow(forwardsTable)){
      
      setTxtProgressBar(pBar, k/nrow(forwardsTable))
      
      #Set flags to nothing found
      t1 <- Sys.time()
      forwardFound <- FALSE
      reverseFound <- FALSE
      targetFound <- FALSE
      
      totalInsertionsForward <- 0
      totalInsertionsReverse <- 0
      
      totalDeletionsForward <- 0
      totalDeletionsReverse <- 0
      
      totalMismatchForward <- 0
      totalMismatchReverse <- 0
      
      fileID <- "none"
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      #Timing Initialize_variables
      matrixTiming[10,2] <- as.numeric(matrixTiming[10,2]) + td
      matrixTiming[10,3] <- as.numeric(matrixTiming[10,3]) + 1
      if(as.numeric(matrixTiming[10,3])<101){
        matrixTiming[10,as.numeric(matrixTiming[10,3])+3] <- td    
      }
      
      #Get the sequences
      t1 <- Sys.time()
      candidateForwardSequence <- forwardsTable[k,"Sequence"]
      candidateReverseSequence <- reversesTable[k,"Sequence"]
      candidateForwardQuality  <- forwardsTable[k,"Quality"]
      candidateReverseQuality  <- reversesTable[k,"Quality"]
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      #Timing Get info
      matrixTiming[7,2] <- as.numeric(matrixTiming[7,2]) + td
      matrixTiming[7,3] <- as.numeric(matrixTiming[7,3]) + 1
      if(as.numeric(matrixTiming[7,3])<101){
        matrixTiming[7,as.numeric(matrixTiming[7,3])+3] <- td
      }
      
      #Get how many N Nuclotides are each sequence. Importat for filtering and statistic porpuses
      totalNReverse <- length(grep("N", strsplit(toString(candidateReverseSequence),"")[[1]], fixed=TRUE))
      totalNForward <- length(grep("N", strsplit(toString(candidateForwardSequence),"")[[1]], fixed=TRUE))
      
      #Filter by the sequence options, if it is interesting continue, otherwise skip to the next one
      filter_pass <- TRUE
      
      skip_N_Nucleotides_Filter <- TRUE
      average_Quality_Filter <- TRUE
      min_Quality_Filter <- TRUE
      
      #Check if we have N nucleotides
      if(skip_N_Nucleotides == TRUE && (totalNReverse>0 || totalNForward >0)){
        skip_N_Nucleotides_Filter <- FALSE
      }
      
      #Get the average quality
      #Transform the quality string into an array of char
      qualityArrayForward <- strsplit(as.character(candidateForwardQuality),"")
      qualityArrayReverse <- strsplit(as.character(candidateReverseQuality),"")
      #For each of those, accumulate and find the average
      accumulateF <- 0
      accumulateR <- 0
      for(l in range(length(qualityArrayForward))){
      
        accumulateF <- strtoi(charToRaw(qualityArrayForward[[1]][l]),16L) + accumulateF
        
      }
      averageF <- accumulateF/length(qualityArrayForward)
      
      for(l in range(length(qualityArrayReverse))){
        
        accumulateR <- strtoi(charToRaw(qualityArrayReverse[[1]][l]),16L) + accumulateR
        
      }
      averageR <- accumulateR/length(qualityArrayReverse)
      
      #Find the minimums
      minimumF <- strtoi(charToRaw(min(qualityArrayForward[[1]])),16L) - 33
      minimumR <- strtoi(charToRaw(min(qualityArrayReverse[[1]])),16L) - 33
        
      #TODO: Transform into a 0,100 scale
      #Check againts the filters    
      if(average_Quality > averageF || average_Quality > averageR){
        average_Quality_Filter <- FALSE
      }
        
      if(min_Quality > minimumF || min_Quality > minimumR){
        min_Quality_Filter <- FALSE
      }
      
      filter_pass = skip_N_Nucleotides_Filter && average_Quality_Filter && min_Quality_Filter
      
      
      #TODO: Have we use this sequence already?
      
      #Copypaste the file with the alignment
      
      
    
      #Search for the forward, reverse and target
      {
      
      t1 <- Sys.time()
      #Search for the forward primer at the begginning of the forward read sequence
      forwardPrimerPositions <- regexpr(forwardPrimer, candidateForwardSequence, fixed=TRUE)[1]
      #for(z in 1:length(forwardPrimerPositions){}
      if(forwardPrimerPositions > 0){
        forwardFound <- TRUE
      }
          
      #Search for the reverse primer at the beggining of the reverse read sequence
      reversePrimerPositions <- regexpr(reversePrimer, candidateReverseSequence, fixed=TRUE)[1]
      if(reversePrimerPositions > 0){
        reverseFound <- TRUE
      }
    
      #Search for the target in the forward read sequence, if is not in the sequence prompt a warning
      targetPrimerPositions <- regexpr(targetPrimer, candidateForwardSequence, fixed=TRUE)[1]
      if(targetPrimerPositions > 0){
        targetFound <- TRUE
      }
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      #Timing Search_for_primers
      matrixTiming[11,2] <- as.numeric(matrixTiming[11,2]) + td
      matrixTiming[11,3] <- as.numeric(matrixTiming[11,3]) + 1
      if(as.numeric(matrixTiming[11,3])<101){
        matrixTiming[11,as.numeric(matrixTiming[11,3])+3] <- td    
      }
      
      }
          
      #Filter out the sequences that do not have a forward primer or reverse primer; warning for those
      #which don't have a target (by writing that on the summary file)
          
      if(filter_pass && (forwardFound && reverseFound)){
                        
        
        t1 <- Sys.time()
        #Find the reverse complementary of the reverse sequence for the multiple alignment
        candidateComplementarySequence <- reverseComplement(DNAStringSet(candidateReverseSequence))

        #Get the forward and reverse candidate and make a DNASet with them
        forrevSequences <- DNAStringSet(c(toString(candidateComplementarySequence)[1], toString(candidateForwardSequence)))
        
        #That DNA set is pairwise aligned with the genome
        myPairwise <- pairwiseAlignment(forrevSequences, genome, type = "global")
        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        #Timing Pairwise_FandR
        matrixTiming[12,2] <- as.numeric(matrixTiming[12,2]) + td
        matrixTiming[12,3] <- as.numeric(matrixTiming[12,3]) + 1
        if(as.numeric(matrixTiming[12,3])<101){
          matrixTiming[12,as.numeric(matrixTiming[12,3])+3] <- td    
        }
        
        #The pairwise object contain some information that we would like to extract
        #Number of insertions, Number of deletions, Number of mistchmatches, Alignment txt file

        #Write the Alignment txt file, put in the folder and add the info to each matrix
        #The folder would be in currentAlignmentsResultsFolderName (ie: run_??/Cas9_data_results/phf8_1_1/Alignments)
        #Ultimatly in this folder we won't have 100K txts, but only those whose primers matched
        
        t1 <- Sys.time()
        fileID <- paste(k,forwardsTable[k,"Instrument_ID"],"_",forwardsTable[k,"Run_ID"],"_",forwardsTable[k,"X"],"_",forwardsTable[k,"Y"],sep='')
        if(writeAlignments == TRUE){
          writePairwiseAlignments(myPairwise, file=paste(currentAlignmentsResultsFolderName,"/","Alignment_",fileID,".fasta",sep=''), block.width=80) # Writes sequences to file in FASTA format.
        }
        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        #Timing Write_Pairwise
        matrixTiming[13,2] <- as.numeric(matrixTiming[13,2]) + td
        matrixTiming[13,3] <- as.numeric(matrixTiming[13,3]) + 1
        if(as.numeric(matrixTiming[13,3])<101){
          matrixTiming[13,as.numeric(matrixTiming[13,3])+3] <- td    
        }
        
        t1 <- Sys.time()
        #Get the number of insertions for Forward and Reverse and put it in each matrix
        totalInsertionsForward <- length(insertion(myPairwise)[[2]])
        totalInsertionsReverse <- length(insertion(myPairwise)[[1]])          

        #Get the number of deletions for Forward and Reverse and put it in each matrix
        totalDeletionsForward <- length(deletion(myPairwise)[[2]])
        totalDeletionsReverse <- length(deletion(myPairwise)[[1]])

        #Get the number of mistchmatches for Forward and Reverse and put it in each matrix

        #Account for the N nucleotides, those shall not be count as a mismatch
        totalMismatchForward <- nmismatch(myPairwise)[2] - totalNForward
        totalMismatchReverse <- nmismatch(myPairwise)[1] - totalNReverse
        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        #Timing Writing_variables
        matrixTiming[17,2] <- as.numeric(matrixTiming[17,2]) + td
        matrixTiming[17,3] <- as.numeric(matrixTiming[17,3]) + 1
        if(as.numeric(matrixTiming[17,3])<101){
          matrixTiming[17,as.numeric(matrixTiming[17,3])+3] <- td    
        }
        

        #Get the RAW insertions, deletions, and mistchmatches information and plug it into another file near the alignment

        #Direct output to a file
        #S4 Strings appears to have a lot of print,cat,write issues
        #this is the only methods that I found that works
        t1 <- Sys.time()
        
        if(writeRAW == TRUE){
          sink(paste(currentAlignmentsResultsFolderName,"/","RAW_",fileID,".fasta",sep=''), append=FALSE, split=FALSE)
  
          print("Insertions (Reverse, Forward)")
          print(insertion(myPairwise))
          print("Deletions (Reverse, Forward)")
          print(deletion(myPairwise))
          print("Mismatches (Reverse, Forward)")
          print(mismatchTable(myPairwise))
          
          # return output to the terminal
          sink()
        }
        
        t2 <- Sys.time()
        td <- as.numeric(t2-t1, units = "secs")
        #Timing Write_RAW
        matrixTiming[14,2] <- as.numeric(matrixTiming[14,2]) + td
        matrixTiming[14,3] <- as.numeric(matrixTiming[14,3]) + 1
        if(as.numeric(matrixTiming[14,3])<101){
          matrixTiming[14,as.numeric(matrixTiming[14,3])+3] <- td    
        }
            
      }
          
      
      #Write the statistics
      {
      #Write the Forwards Stats
        
      t1 <- Sys.time()
        
      matrixStatsForwardsTable[k,8] = forwardFound
      matrixStatsForwardsTable[k,9] = reverseFound
      matrixStatsForwardsTable[k,10] = targetFound
      matrixStatsForwardsTable[k,11] = totalInsertionsForward
      matrixStatsForwardsTable[k,12] = totalDeletionsForward
      matrixStatsForwardsTable[k,13] = totalMismatchForward
      matrixStatsForwardsTable[k,14] = paste("Alignment_",fileID,sep='')

      #Write the Reverse Stats
      matrixStatsReversesTable[k,8] = forwardFound
      matrixStatsReversesTable[k,9] = reverseFound
      matrixStatsReversesTable[k,10] = targetFound
      matrixStatsReversesTable[k,11] = totalInsertionsReverse
      matrixStatsReversesTable[k,12] = totalDeletionsReverse
      matrixStatsReversesTable[k,13] = totalMismatchReverse
      matrixStatsReversesTable[k,14] = paste("Alignment_",fileID,sep='')
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      #Timing Fill_datastructure
      matrixTiming[5,2] <- as.numeric(matrixTiming[5,2]) + td
      matrixTiming[5,3] <- as.numeric(matrixTiming[5,3]) + 1
      if(as.numeric(matrixTiming[5,3])<101){
        matrixTiming[5,as.numeric(matrixTiming[5,3])+3] <- td
        
      }
      }
    }
   
    #Now write a txt file with the statistics tables for the forwards and reverse alignments
    #It goes into currentIDCodeResultsFolderName (ie: run_??/Cas9_data_results/phf8_1_1)
    {
    #Convert matrices into dataframes
    t1 <- Sys.time()
    statsForwardsTable <- data.frame(matrixStatsForwardsTable)
    colnames(statsForwardsTable) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Forward_Found","Reverse_Found","Target_Found","Total_Insertions","Total_Deletions","Total_Mismatch","Alignment_file")

    statsReversesTable <- data.frame(matrixStatsReversesTable)
    colnames(statsReversesTable) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Forward_Found","Reverse_Found","Target_Found","Total_Insertions","Total_Deletions","Total_Mismatch","Alignment_file")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Transform_into_datafame
    matrixTiming[6,2] <- as.numeric(matrixTiming[6,2]) + td
    matrixTiming[6,3] <- as.numeric(matrixTiming[6,3]) + 1
    if(as.numeric(matrixTiming[6,3])<101){
      matrixTiming[6,as.numeric(matrixTiming[6,3])+3] <- td
      
    }

    #Need to convert characters into proper numbers
    t1 <- Sys.time()
    statsForwardsTable[, 11] <- as.numeric(as.character( statsForwardsTable[, 11] ))
    statsForwardsTable[, 12] <- as.numeric(as.character( statsForwardsTable[, 12] ))
    statsForwardsTable[, 13] <- as.numeric(as.character( statsForwardsTable[, 13] ))
    
    statsReversesTable[, 11] <- as.numeric(as.character( statsReversesTable[, 11] ))
    statsReversesTable[, 12] <- as.numeric(as.character( statsReversesTable[, 12] ))
    statsReversesTable[, 13] <- as.numeric(as.character( statsReversesTable[, 13] ))
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Fill_datastructure
    matrixTiming[5,2] <- as.numeric(matrixTiming[5,2]) + td
    matrixTiming[5,3] <- as.numeric(matrixTiming[5,3]) + 1
    if(as.numeric(matrixTiming[5,3])<101){
      matrixTiming[5,as.numeric(matrixTiming[5,3])+3] <- td
      
    }
    
    t1 <- Sys.time()
    #Write the file
    write.table(statsForwardsTable, paste(currentIDCodeResultsFolderName,"/",currentID,"_",currentBarcode,"_forwards.txt",sep = ''), sep="\t") 
    write.table(statsReversesTable, paste(currentIDCodeResultsFolderName,"/",currentID,"_",currentBarcode,"_reverses.txt",sep = ''), sep="\t") 
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Write_STATS
    matrixTiming[15,2] <- as.numeric(matrixTiming[15,2]) + td
    matrixTiming[15,3] <- as.numeric(matrixTiming[15,3]) + 1
    if(as.numeric(matrixTiming[15,3])<101){
      matrixTiming[15,as.numeric(matrixTiming[15,3])+3] <- td    
    }
    }

    #Now we are finish with a line in the config file, write the total stats
    #That have the following columns [ID,Barcode,FR_file,RR_file,F_Insertions,F_Deletions,F_Mismatches,R_Insertions,R_Deletions,R_Mismatches]
    t1 <- Sys.time()
    matrixConfigFileResults[j,1] <- currentID
    matrixConfigFileResults[j,2] <- currentBarcode
    matrixConfigFileResults[j,3] <- forwardReadsFileName
    matrixConfigFileResults[j,4] <- reverseReadsFileName
    matrixConfigFileResults[j,5] <- sum(statsForwardsTable$"Total_Insertions")
    matrixConfigFileResults[j,6] <- sum(statsForwardsTable$"Total_Deletions")
    matrixConfigFileResults[j,7] <- sum(statsForwardsTable$"Total_Mismatch")
    matrixConfigFileResults[j,8] <- sum(statsReversesTable$"Total_Insertions")
    matrixConfigFileResults[j,9] <- sum(statsReversesTable$"Total_Deletions")
    matrixConfigFileResults[j,10] <- sum(statsReversesTable$"Total_Mismatch")
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    #Timing Fill_datastructure
    matrixTiming[5,2] <- as.numeric(matrixTiming[5,2]) + td
    matrixTiming[5,3] <- as.numeric(matrixTiming[5,3]) + 1
    if(as.numeric(matrixTiming[5,3])<101){
      matrixTiming[5,as.numeric(matrixTiming[5,3])+3] <- td
    }

    }
   
    t4 <- Sys.time()
    td2 <- as.numeric(t4-t3, units = "secs")
    #Timing Config_Line_Complete
    matrixTiming[16,2] <- as.numeric(matrixTiming[16,2]) + td2
    matrixTiming[16,3] <- as.numeric(matrixTiming[16,3]) + 1
    if(as.numeric(matrixTiming[16,3])<101){
      matrixTiming[16,as.numeric(matrixTiming[16,3])+3] <- td2
    }


    #Get the timing into a dataframe
    timing <- data.frame(matrixTiming, stringsAsFactors = FALSE)
    colnames(timing) <- c("Concept","Accumulated_Time_s","Total_Observations")
    for(k in 2:103){
      timing[, k] <- as.numeric(as.character( timing[, k] ))
    }
    
    #For each of the concept in the dataframe of timing
    #Write a line in matrixAbsolute with the statistics for that concept
    allConcepts <- levels(factor(timing$Concept))
    
    for (k in 1:length(allConcepts)){
      
      print(timing[k,1])
      
      #Find the first 0 in the row, that will be the last sample to take (excluded)
      index <- 103
      for (n in 5:103){ #We know that there is at least one, hence start from sample 2 (column 5)
        if(timing[k,n]==0){
          index <- (n-1)
          break
        }
      }
      
      #Get the vector with [concept, average, sd, 0Q, 1Q, 2Q, 3Q, 4Q(MAX),Absolute_Frequency]
      samplesSubset <- NULL
      if(index==4){
        samplesSubset <- timing[k,4]
      }
      else{
        samplesSubset <- subset(timing[k,c(4:index)])
      }
      
      samplesVector <- as.vector(t(samplesSubset))
      sResults <- summary(samplesVector)
      
      statVector <- c(timing[k,1],sResults[[4]],sd(samplesVector),sResults[[1]],sResults[[2]],sResults[[3]],sResults[[5]],sResults[[6]],timing[k,3])
      
      
      #Each of those vectors goes into the absolute matrix with rbind
      #If the average is 0 means that we didn't got any reading on that (for example, no candidates were pairwise)
      #In that case we bind nothing
      if(sResults[[4]]>0){
        matrixTimingAbsolute <- rbind(matrixTimingAbsolute, statVector)
      }
      
    }



  }
  
  #Now that everything is finish for this data folder
  
  #Write the final Config_file_i_results.txt (ie: run11_results.txt)
  #That have the following columns [ID,Barcode,FR_file,RR_file,F_Insertions,F_Deletions,F_Mismatches,R_Insertions,R_Deletions,R_Mismatches]
  
  #First create the dataset from the matrix
  configFileResults <- data.frame(matrixConfigFileResults)
  colnames(configFileResults) <- c("ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Forwards_Insertions","Forwards_Deletions","Forwards_Mismatches","Reverse_Insertions","Reverse_Deletions","Reverse_Mismatches")  

  #Now write the file into a .txt file in currentResultsFolderName (ie: run_??/Cas9_data_results/)
  write.table(configFileResults, paste(currentResultsFolderName,"/","configFile_results",sep = ''), sep="\t")   
  
  
  #We have now finish the timing. Now write some statistics and plot something relevant.
  #The graphics of these goes into the folder with the results.
  # The absolute timing columns look like this: [concept, average, sd, 0Q, 1Q, 2Q, 3Q, 4Q(MAX)]

  #Get the matrix into a dataframe
  timingAbsolute <- data.frame(matrixTimingAbsolute, stringsAsFactors = FALSE)
  colnames(timingAbsolute) <- c("Concept","average","sd","MIN","1Q","2Q","3Q","MAX","Absolute_Frequency")
  for(k in 2:9){
    timingAbsolute[, k] <- as.numeric(as.character( timingAbsolute[, k] ))
  }

  #Show the results to the users
  allConcepts <- levels(factor(timingAbsolute$Concept))
  for (j in 1:length(allConcepts)){
   
    print(allConcepts[j])
    print(summary(subset(timingAbsolute, Concept==allConcepts[j])$average))
    print(nrow(subset(timingAbsolute, Concept==allConcepts[j])))
    
  }

  #Write the raw timing results in disk
  write.table(timingAbsolute, paste(currentResultsFolderName,"/","timing_results.txt",sep = ''), sep="\t")   


  #Lets find out the proportion of each concept with respect the config lines

  #Make a matrix where to rbind the percentages
  matrixPercentage <- NULL
  
  #Find the average time taken to complete a config line
  configLineAverageTime = summary(subset(timingAbsolute, Concept=="Config_Line_Complete")$average)[[4]]

  #Get a subset of all the data that is not a config line time
  noConfigDataFrame <- subset(timingAbsolute, Concept!="Config_Line_Complete")

  #For each of those concept, find the proportion with respect the time of the config line
  allConcepts <- levels(factor(noConfigDataFrame$Concept))
  for (j in 1:length(allConcepts)){
    
    ssConcept <- subset(noConfigDataFrame, Concept==allConcepts[j])
    
    for(k in 1:nrow(ssConcept)){
      
      matrixPercentage <- rbind(matrixPercentage, c(allConcepts[j],(ssConcept[k,]$average * ssConcept[k,]$Absolute_Frequency)/configLineAverageTime))
      
    }
    
  }

  #When everything is finish show statistics of those concepts and %
  percentages <- data.frame(matrixPercentage)
  colnames(percentages) <- c("Concept","Time_s")
  percentages[, 2] <- as.numeric(as.character( percentages[, 2] ))
  allConcepts <- levels(factor(percentages$Concept))
  
  for (i in 1:length(allConcepts)){
    print(allConcepts[i])
    
    print(quantile(subset(percentages, Concept==allConcepts[i])$Time_s))
    print(summary(subset(percentages, Concept==allConcepts[i])$Time_s))
    print(nrow(subset(percentages, Concept==allConcepts[i])))
  }
  
  #Transform into a boxplot graphic

  print("Duplicates?")

  test.m <- melt(percentages)
  
  myBoxplot <- ggplot(test.m, aes(factor(variable), value)) + geom_boxplot(aes(fill = Concept))
  
  svg(paste(currentResultsFolderName,"/","plots.svg",sep = ''))
  myBoxplot
  dev.off()


}