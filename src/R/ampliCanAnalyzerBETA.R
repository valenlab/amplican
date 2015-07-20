############################################
# Constant declarations
############################################
{
  TOTAL_PROCESSORS <- 1       # Set this to the number of processors you want to use.
                               # You can set it to NULL if you want all that are available.
    
  ALL_IN_MODE <- 3            # This variable tells which critiria do we follow to decide if only one indel, or many
                              # indels, or all indel, or no indel is within range of one, all, or none alignment
                              # position and many combinations; please check the documentation in cutscriterias.R
                              # Mode 1 means that all cuts must be within reach
                              # Mode 3 means that at least one cut must be within reach
                              # Mode 5 means we don't care about the reach, so effectively turning this criteria off
  
  ALIGNMENT_DISTANCE <- 5     # When looking for valid indels, the indel must be with in the range of the alignment
                              # position. The range is -ALIGNMENT_DISTANCE from the left , and
                              # + (ALIGNMENT_WIDE + this distance) to the right
  
  SAME_CUTS  <- TRUE        # This variable tells if the valid events must be in the same absolute position in the
                            # Forward alignment and the Reverse alignment
  
  PRIMERDIMER <- TRUE       # This variable tells if we should filter primer dimers occurences.
  
  DIMERBASE <- 10           # The legnth of the constant bases that we use to adjust the primer dimer situation
  
  TIMING <- TRUE            # Incluse the timing results inside the final folder. This is intended for debugging.
  
  # Test
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150424215630_run11" # Toy data
  
  # Danielle Old?
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150217101758_20141112_F"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150219083907_20141112_R"
  
  # Jason
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150410175748_config_Jason2"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150410180646_JasonAfterChristmas"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150506233015_Jason3"
  TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150508172357_Jason4" # Jason last corrected
  
  
  # Danielle
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150424220252_20141112_F"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150321043346_20141112_R"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150321043622_20141117_F2"
  #TARGET_FOLDER <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150321043549_20141117_R2"
  
  
}

############################################
#Summons the libraries and get to use functions from there.
############################################
{  

  # Making things in parallel
  library(doParallel)
  cl <- makeCluster(TOTAL_PROCESSORS, outfile="") #outfile tells the workers to print on the terminal
  registerDoParallel(cl)
  
  # Our own developed stuff.
  source("libraries/cutscriterias2.R") # How deletions are defined.
  source("libraries/tools2.R") # Minor stuff like the reverse complement of a DNA sequence.
  source("libraries/configAnalyze2.R") # Function that deals with the config files and get running everything.
  
}

############################################
# Algorith starts here
############################################

print("START")

# We are going to write the results in the same folder, inside another folder called:
# Analysis_<analysis variables>_YYYYmmddHHMMSS
{
  currentTime <- Sys.time()
  timeStamp <-  strftime(currentTime,"%Y%m%d%H%M%S")
  analysisPath <- paste(TARGET_FOLDER, "/Analysis_",ALL_IN_MODE,"_",ALIGNMENT_DISTANCE,"_",
                        SAME_CUTS,"_", PRIMERDIMER,"_", DIMERBASE,"_",timeStamp, sep = '')
  dir.create(file.path(analysisPath), showWarnings = TRUE)
  
  # Inside this folder it will be a file named log where we will write relevant information
  logFileName <- paste(analysisPath, "/AnalysisLog.txt", sep = '')
  logFileConn <- file(logFileName, open="at")
  
  # Write the analysis variables into the log
  writeLines("The analysis was made with the following variables: ", logFileConn)
  writeLines(paste("ALL IN MODE           : ", ALL_IN_MODE), logFileConn)
  writeLines(paste("ALIGNMENT DISTANCE    : ", ALIGNMENT_DISTANCE), logFileConn)
  writeLines(paste("SAME CUTS             : ", SAME_CUTS), logFileConn)
  writeLines(paste("PRIMER DIMER          : ", PRIMERDIMER), logFileConn)
  writeLines(paste("PRIMER DIMER CONSTANT : ", DIMERBASE), logFileConn)
  writeLines("", logFileConn) # Empty line
}

# TODO:
# Check out that we have the following:
# - A folder called *_results
# - Inside that folder there is a file called config_results.txt
# - The number of folders inside *_results is the same as the number
#   of lines in the .txt file (not counting the header)
{

  # Get the result folder
  resultsPath     <- list.files(path = TARGET_FOLDER, pattern="_alignments", full.names = TRUE)
  allPaths        <- NULL
  allFoldersNames <- NULL
  configPath      <- NULL
  
  if(length(resultsPath) > 0){
    # Get the list of folders
    allPaths        <- list.dirs(path = resultsPath, full.names = TRUE,  recursive=FALSE)  
    allFoldersNames <- list.dirs(path = resultsPath, full.names = FALSE, recursive=FALSE)    
    configPath      <- list.files(path = resultsPath, pattern="config_results.txt", full.names = TRUE)
  
  }
  
  # If something fails, display on screen, write it on the log, and exit
  else{
    
    print("ERROR!: No result path found")
    print("ERROR!: No config .txt found?")
  
  }

}

# If we have everything, keep going
if(!is.null(allPaths) && !is.null(allFoldersNames) && !is.null(configPath)){

  # Load the config file into a dataframe
  configDF <-  read.table(configPath, header=TRUE, sep="\t")
  
#   # Prepare the barcode dataframe and write it into disk
#   # We need to forward this information to the plots
#   keeps <- c("Barcode", "Total_Pre_N_Filter", "Total_Post_N_Filter", "Experiment_Unique_Sequences", "Total_Experiment_Sequences")
#   barcodeDF <- configDF[keeps]
#   barcodeDFOriginal <- configDF[keeps]
#   barcodeDF <- unique(barcodeDF)
#   write.table(barcodeDF, paste(analysisPath,"/barcode_results.txt", sep=''), sep="\t")
  
  # Clean up the Dataframe, we don't need most of the columns
  keeps <- c("ID","Barcode", "Forward_Primer","Reverse_Primer", "Genome", "Sum_Is_Sequence")
  configDF <- configDF[keeps]
  
  #print(configDF)
  
  # Now we divide the work among the diffent processors.
  processorsSubDataframes <- divideWorkBySize(TOTAL_PROCESSORS, configDF,"Sum_Is_Sequence")

  # Analyze the config file in parallel 
  foreach(j=1:TOTAL_PROCESSORS , .packages=c()) %dopar% {
    
    analyzeConfig(ALL_IN_MODE, ALIGNMENT_DISTANCE, SAME_CUTS, PRIMERDIMER, DIMERBASE,
                  BOTHDIMERS, processorsSubDataframes[[j]], analysisPath, j, resultsPath, TIMING)
    
    
  }
  
  # Put all the logs and all the configs toguether
  totalLogs    <- unifyFiles(analysisPath,             "SUBLOG",             "myFinalLog.txt",          FALSE,TRUE,FALSE)
  totalConTime <- unifyFiles(currentResultsFolderName, "configFile_Timing",  "analysisFile_Timing.txt", TRUE, TRUE,FALSE)
  totalConfigs <- unifyFiles(analysisPath,             "configFile_results", "analysis_results.txt",    TRUE,TRUE,FALSE)
  
  print("Finish analysis")
  
  
}

# Close the master log file descriptor
close(logFileConn)

# Stop the cluster and go home
stopCluster(cl)

print("Finish!")
print(paste("Look for your results in", analysisPath))