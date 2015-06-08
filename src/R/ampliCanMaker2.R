############################################
# Constant declarations
############################################
{
TOTAL_PROCESSORS <- 1       # Set this to the number of processors you want
                             # to use. You can set it to NULL if you want all
                             # that are available.
  
SKIP_BAD_NUCLEOTIDES <- TRUE  # Some sequences have faulty nucleotides labels
                              # with N. If we find a sequence like that in
                              # either forwards or reverse, we skip that
                              # aligment
  
AVERAGE_QUALITY <- 0        # The FASTQ file have a quality for each nucleotide, being ! the lower and ~ the highest
                            # In ASCII :
                            # !"#$%&'()*+,-./0123456789:;<=>?@
                            # ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
                            # abcdefghijklmnopqrstuvwxyz{|}~
                            # This quality number goes from 0 to 100, being 0
                            # the lowest quality and 100 the highest
                            # You can write whatever number in between and the
                            # program will find out the apropiate character
                            # encoding. The filter work by converting each
                            # character to a number, and then find the average.
                            # If the average fall above this threshold then we
                            # take the sequence.
  
MIN_QUALITY <- 0            # Same as before, but this is the minimum quality
                            # for ALL nucleotide. If one of them
                            # has quality BELLOW this threshold, then
                            # the sequence is skipped.
  
WRITE_ALIGNMENTS            <- 2       # Write the aligments results into disk
                                       # 0 - Write nothing
                                       # 1 - Write only the summary file
                                       # 2 - Write also a verbose file with all the alignments in the same .txt
                                       # 3 - Write also every individual .txt file

  

SCORING_MATRIX <- "NUC44"  # The scoring matrix that you wish to use
GAP_OPENING    <- 50       # The opening gap score 
GAP_EXTENSION  <- 0        # The gap extension score
GAP_ENDING     <- FALSE    # If you want that the ending gap count for the scoring set this to TRUE
FAR_INDELS     <- TRUE     # If AA-- AAAA should be considered to be an indel from 3 to 4.

TIMING <- TRUE
    
#CONFIG_LIST <- c("/export/kjempetujafs/valen-group/projects/Cas9/Danielle/20141110_Run_1/20141112_F.conf") # Danielle Original Forward
#CONFIG_LIST <- c("/export/kjempetujafs/valen-group/projects/Cas9/Danielle/20141110_Run_1/20141112_R.conf") # Danielle Original Reverse
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/configs/Danielle/20141117_F2.txt") # Danielle Christmas Forward
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/configs/Danielle/20141117_R2.txt") # Danielle Christmas Reverse

CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/git/ampliCan/res/Cas9_toy/run11.conf") # Toy data

#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Jason/config_Jason2.txt") # Jason Original
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Jason2/150212-MiSeq-20246231/JasonAfterChristmas.txt") # Jason Chritmas
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Jason3/all/Jason3.txt") # Jason April
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Jason3/all/Jason4.txt") # Jason April Corrected


#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Built/wernerConfig.txt") # Werner
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/miseq3_DSP1.conf") # Werner SN0048153
#CONFIG_LIST <- c("/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/miseq3_DSP2.conf") # Werner SN0048152
  
}

############################################
#Summons the libraries and get to use functions from there.
############################################
{
  
# Utilities
library(R.utils) # Unzipping the .gz files
library(plyr)
library(dplyr) # Makes code look nicer; implemented hash for duplicated sequences with this.
library(Rcpp) # Alignments
library(gotoh)

# Making things in parallel
library(doParallel)
cl <- makeCluster(TOTAL_PROCESSORS, outfile="") #outfile tells the workers to print on the terminal
registerDoParallel(cl)
  
# Our own developed stuff.
#sourceCpp("libraries/gotoh.cpp") # Alignments.

source("libraries/filters2.R") # Filtering.
#source("libraries/mailing.R")
#source("libraries/inandout2.R") # How to read and write from the disk.
#source("libraries/cutscriterias.R") # How deletions are defined. TODO: This is not here use right?
source("libraries/tools2.R") # Minor stuff like the reverse complement of a DNA sequence.
source("libraries/errorswarnings2.R") # Handles pre-parsing, error messages, warning to the users, and so on.
source("libraries/configAlignments2.R") # Function that deals with the config files and get running everything.
#source("libraries/recovery.R") #TODO Do we need this anymore ?
}

############################################
# Algorith starts here
############################################

print("START")

currentConfigPath <- CONFIG_LIST[1]

# Check that the config file exist and is readable, etc...
readyConfig <- checkReady(currentConfigPath)

# If it is a bad config show the error and stop the program; otherwise, continue working on it
if(readyConfig == FALSE){
  
  print("-----------------------------------------")
  print("---|             ERROR!              |---")
  print("-----------------------------------------")
  print("ERROR!: Couldn't read this config file.")
  print(currentConfigPath)
  print("Does it exist? Do you have read access?")
  
}else{

  # Get the name of the config file
  configFileName <- getConfigName(currentConfigPath)
  
  # Get the path of the config folder
  configFolderName <- getConfigFolder(currentConfigPath)
  
  # Prepare the folder with the results
  # Check if there is a folder call /result in cwd, if not create a new one
  dir.create(file.path(getwd(), "/results"), showWarnings = FALSE)
  
  # Make a folder call /results/run_YYYYmmddHHMMSS_<config name>
  # where we are going to write the results for this run.
  {
    currentTime <- Sys.time()
    timeStamp <-  strftime(currentTime,"%Y%m%d%H%M%S")
    resultsFolder <- paste(getwd(), "/results/run_", timeStamp, "_", configFileName, sep = '')
    dir.create(file.path(resultsFolder), showWarnings = TRUE)
    
    # Inside this folder it will be a file named log where we will write relevant information
    logFileName <- paste(resultsFolder, "/MasterLog.txt", sep = '')
    logFileConn <- file(logFileName, open="at")
    
    # Write the constants into the master logs
    writeLines(paste("Total Processors:      " , TOTAL_PROCESSORS) , logFileConn)
    writeLines(paste("Skip Bad Nucleotides:  " , SKIP_BAD_NUCLEOTIDES) , logFileConn)
    writeLines(paste("Average Quality:       " , AVERAGE_QUALITY) , logFileConn)
    writeLines(paste("Minimum Quality:       " , MIN_QUALITY) , logFileConn)
    writeLines(paste("Write Alignments Mode: " , WRITE_ALIGNMENTS) , logFileConn)
    writeLines(paste("Scoring Matrix:        " , SCORING_MATRIX) , logFileConn)
    writeLines(paste("Gap Opening:           " , GAP_OPENING) , logFileConn)
    writeLines(paste("Gap Extension:         " , GAP_EXTENSION) , logFileConn)
    writeLines(paste("Gap Ending:            " , GAP_ENDING) , logFileConn)
    writeLines(paste("Far Indels:            " , FAR_INDELS) , logFileConn)
    
  }
  
  # Later on, we are going to generate analysis based on this alignments, plots, etc...
  # To avoid a mess in the folder, we are going to save the alignments a folder called
  # /run_XYZ/configFileName_alignments
  currentResultsFolderName <- paste(resultsFolder,"/",configFileName,"_alignments",sep = '')
  dir.create(file.path(currentResultsFolderName), showWarnings = TRUE)
  
  # In particular, for the plots we are going to need information of all the sequences that
  # are not assigned to anything. This information is at barcode level, so we are going
  # to make an special folder for this information
  currentUnassignedFolderName <- paste(currentResultsFolderName,"/unassigned_sequences",sep = '')
  dir.create(file.path(currentUnassignedFolderName), showWarnings = TRUE)
  
  # Read the config file into a dataset
  print("Creating config file dataframe")
  configTable <- read.table(currentConfigPath, header=FALSE, sep="\t")
  colnames(configTable) <- c("ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Experiment_Type",
                             "Target_Primer","Forward_Primer","Reverse_Primer","Strand","Genome")
  
  # Check that the config file is correct and we have access to everything
  # Pre-process the config file, if everything is right, continue. Otherwise stop the program
  configValid <- checkConfigFile(configTable, currentConfigPath, logFileConn)
  
  if(configValid == FALSE){
    
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
  
  else{
  
    # Divide the config dataframe into several subset of roughly equal size
    processorsSubDataframes <- divideWorkByBarcode(TOTAL_PROCESSORS, configTable, currentConfigPath)
    
    # Analyze the config file in parallel 
    parallelPackages = c("Rcpp","R.utils","plyr","dplyr","gotoh")
    
    foreach(j=1:TOTAL_PROCESSORS , .packages=parallelPackages) %dopar% {
      
      makeAlignment5(SKIP_BAD_NUCLEOTIDES, AVERAGE_QUALITY, MIN_QUALITY, WRITE_ALIGNMENTS,
                    SCORING_MATRIX, GAP_OPENING, GAP_EXTENSION, GAP_ENDING,FAR_INDELS,
                    processorsSubDataframes[[j]], resultsFolder, currentResultsFolderName, j,
                    currentConfigPath, TIMING)
      
    }
    
    print("Finish make alignments")
    
    # Put all the logs and all the configs toguether
    totalLogs    <- unifyFiles(resultsFolder,            "SUBLOG"             , "alignmentLog.txt"        ,FALSE,TRUE,FALSE)
    totalConfigs <- unifyFiles(currentResultsFolderName, "configFile_results" , "config_results.txt"      ,TRUE, TRUE,FALSE)
    totalConTime <- unifyFiles(currentResultsFolderName, "configFile_Timing"  , "configFile_Timing.txt"   ,TRUE, TRUE,FALSE)
    totalBarcode <- unifyFiles(currentResultsFolderName, "barcodeFile_results", "barcodeFile_results.txt" ,TRUE, TRUE,FALSE)
    totalBarTime <- unifyFiles(currentResultsFolderName, "barcodeFile_Timing" , "barcodeFile_Timing.txt"  ,TRUE, TRUE,FALSE)
    
    
    
  }
  
  # Clsoe the master log file descriptor
  close(logFileConn)
  
  # Stop the cluster and go home
  stopCluster(cl)
  
  print("Finish!")
  print(paste("Look for your results in",resultsFolder))
  
  t2 <- Sys.time()
  td <- as.numeric(t2-currentTime, units = "secs")
  
  print(paste("Total time: ",td * TOTAL_PROCESSORS ))

}
