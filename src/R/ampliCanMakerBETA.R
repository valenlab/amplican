############################################
# Constant declarations
############################################

# This function takes a configuration file and generate the alignments.

# Parameters:
#
# (string) CONFIG         The path to your configuration file. For example:
#
#                         /Home/johndoe/.../AmpliCan/res/Cas9_toy/run11.conf
#
#
# (int) TOTAL_PROCESSORS  Set this to the number of processors you want to use.
#                         You can set it to NULL if you want all that are
#                         available. Default = 1.
#
# (bool)SKIP_BAD_NUCLEOTIDES  Some sequences have faulty nucleotides labels
#                             with N. If we find a sequence like that in
#                             either forwards or reverse, we skip that
#                             aligment. Default = TRUE
#
# (int) AVERAGE_QUALITY The FASTQ file have a quality for each nucleotide,
#                       being ! the lower and ~ the highest. In ASCII :
#
#                              !"#$%&'()*+,-./0123456789:;<=>?@
#                              ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`
#                              abcdefghijklmnopqrstuvwxyz{|}~
#
#                        This quality variable goes from 0 to 100, being 0
#                        the lowest quality and 100 the highest. You can write
#                        whatever number in between, and the program will find
#                        out the apropiate character encoding. The filter works
#                        by converting each character to a number, and then
#                        finding the average. If the average fall above this
#                        threshold then we take the sequence. Default is 0.
#
# (int) MIN_QUALITY      Same as before, but this is the minimum quality for
#                        ALL nucleotide. If one of them has quality BELLOW this
#                        threshold, then the sequence is skipped. Default is 0.
#
# (int) WRITE_ALIGNMENTS  How to write the aligments results into disk:
#                         0 - Write nothing
#                         1 - Write only the summary file
#                         2 - Write also a verbose file with all the alignments
#                             in the same .txt file (Default option)
#                         3 - Write also every individual .txt file. Be aware
#                             that this option generated thousands of tiny
#                             files which can bottleneck the run.
#
# (string) SCORING_MATRIX The scoring matrix that you wish to use. The
#                         available options are listed here:
#                         - NUC44 (Default)
#
# (int) GAP_OPENING       The opening gap score. Default is 50.
#
# (int) GAP_EXTENSION     The gap extension score. Default is 0.
#
# (bool) GAP_ENDING       If you want that the ending gap count for the
#                         alignment score set this to TRUE. Default is FALSE.
#
# (bool) FAR_INDELS       If the ending/starting gap should be considered to be
#                         an indel. Default is TRUE.
#
#                         NOTE: If you want to filter these out from the plots,
#                               there is an option to do so in the plot
#                               function. You don't need to do it here.
#
# (string) RESULT_FOLDER  Where do you want your results to be stored. The
#                         program will make a folder inside RESULT_FOLDER. ie:
#
#                         RESULT_FOLDER:
#                             /Home/johndoe/myCRISPRData
#
#                         Run data for CONFIGFILE run11.conf
#                             /Home/johndoe/myCRISPRData/run_20150723103145_run11
#                         Run data for CONFIGFILE myconfig.txt
#                             /Home/johndoe/myCRISPRData/run_20150723103145_myconfig
#
#                         If nothing is specified, the folder will be call
#                         "results" and will be placed at the same folder as
#                         the CWD().
#
#                         The folder created for this run will have a timestamp
#                         with this format:
#                         run_ + YYYYMMDDhhmmss + _ + CONFIGFILE without extension
#
# (bool) DELETEFQ         Your FASTQ files can be compressed in a file. If this
#                         happens, the program needs to uncompress them first.
#                         When the program is finish, you can choose to delete
#                         the new uncompressed files. However, this may not be
#                         desirable. If you plan to run it again with different
#                         options, you will need to uncompress everything again
#                         instead of only one time.
#
#                         Set this to TRUE if you want to automatically delete
#                         all the temporal files at the end of the run. Set
#                         this to FALSE, if you don't want to.
#
#                         If your original files were not compressed, this
#                         option won't delete them, even if you set it to TRUE.
#                         You will need to manually delete them later. 
#                         Default is FALSE.
#
# (string) TEMPFOLDER     Your FASTQ files can be compressed in a file. If this
#                         happens, the program needs to uncompress them first.
#                         In order to do so, we need a folder where they will
#                         be placed. In here, you can specify the path where
#                         you want this files to be placed. If you don't
#                         specify a path, they will be placed in the same
#                         folder where the original files are. In any case you
#                         will need writing permissions in the folder in order
#                         to uncompress everything.
#
# (bool) TIMING           You can choose to generate a timing file for
#                         debugging and benchmarking. Default is FALSE.

ampliCanMaker <- function (CONFIG, TOTAL_PROCESSORS = 1,
                           SKIP_BAD_NUCLEOTIDES = TRUE, AVERAGE_QUALITY = 0, 
                           MIN_QUALITY = 0, WRITE_ALIGNMENTS = 2,
                           SCORING_MATRIX = "NUC44", GAP_OPENING = 50,
                           GAP_EXTENSION = 0, GAP_ENDING = FALSE, 
                           FAR_INDELS = TRUE, RESULT_FOLDER = NULL,
                           DELETEFQ = FALSE, TEMPFOLDER = NULL, 
                           TIMING = FALSE){

  ############################################
  #Summons the libraries and get to use functions from there.
  ############################################
  {
    
  # Utilities
  library(R.utils) # Unzipping the .gz files
  library(plyr)
  library(dplyr) # Makes code look nicer; implemented hash for duplicated sequences with this.
  library(Rcpp) # Alignments
  library(Gotoh)
  
  # Making things in parallel
  library(doParallel)
  cl <- makeCluster(TOTAL_PROCESSORS, outfile="") #outfile tells the workers to print on the terminal
  registerDoParallel(cl)
  
  # Our own developed stuff.
  source("libraries/filters.R") # Filtering.
  source("libraries/tools.R") # Minor stuff like the reverse complement of a DNA sequence.
  source("libraries/errorswarnings.R") # Handles pre-parsing, error messages, warning to the users, and so on.
  source("libraries/configAlignments.R") # Function that deals with the config files and get running everything.

  }

  ############################################
  # Algorith starts here
  ############################################

  print("START")
  
  currentConfigPath <- CONFIG

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
    
  }
  else{

    # Get the name of the config file
    configFileName <- getConfigName(currentConfigPath)
    
    # Get the path of the config folder
    configFolderName <- getConfigFolder(currentConfigPath)
    
    # Prepare the folder with the results
    resultsFolder <- ""
    currentTime <- Sys.time()
    timeStamp <-  strftime(currentTime,"%Y%m%d%H%M%S")
    
    # If nothing is specified, the folder will be call results" and will be
    # placed at the same folder as this script.
    if(is.null(RESULT_FOLDER)){
      # Check if there is a folder call /result in cwd, if not create a new one
      dir.create(file.path(getwd(), "/results"), showWarnings = FALSE) 
      resultsFolder <- paste(getwd(), "/results/run_", timeStamp, "_", configFileName, sep = '')
    
    }
    else{
      resultsFolder <- paste(RESULT_FOLDER, "/run_", timeStamp, "_", configFileName, sep = '')
    }
    
    # Make a folder call /RESULT_FOLDER/run_YYYYmmddHHMMSS_<config name>
    # where we are going to write the results for this run.
    {
      dir.create(file.path(resultsFolder), showWarnings = TRUE)
      
      # Inside this folder it will be a file named log where we will write relevant information
      logFileName <- paste(resultsFolder, "/MasterLog.txt", sep = '')
      logFileConn <- file(logFileName, open="at")
      
      # Inside that folder, it will be another file, that keep track of all temporal files
#       tempFileName <- paste(resultsFolder, "/TemporalFiles.txt", sep = '')
#       tempFileConn <- file(tempFileName, open="at")
      
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
      
      writeConfigError(currentConfigPath, logFileName, logFileConn)
      
    }
    else{
    
      # Divide the config dataframe into several subset of roughly equal size
      processorsSubDataframes <- divideWorkByBarcode(TOTAL_PROCESSORS, configTable, currentConfigPath)
      
      # Analyze the config file in parallel 
      parallelPackages = c("Rcpp","R.utils","plyr","dplyr","Gotoh")
      
      foreach(j=1:TOTAL_PROCESSORS , .packages=parallelPackages) %dopar% {
        
        source("libraries/filters.R") # Filtering.        
        source("libraries/tools.R") # Minor stuff like the reverse complement of a DNA sequence.        
        source("libraries/errorswarnings.R") # Handles pre-parsing, error messages, warning to the users, and so on.        
        source("libraries/configAlignments.R") # Function that deals with the config files and get running everything.
        
        #print(reverseComplement("aaaatcgat"))
        
        makeAlignment(SKIP_BAD_NUCLEOTIDES, AVERAGE_QUALITY, MIN_QUALITY, WRITE_ALIGNMENTS,
                      SCORING_MATRIX, GAP_OPENING, GAP_EXTENSION, GAP_ENDING,FAR_INDELS,
                      processorsSubDataframes[[j]], resultsFolder, currentResultsFolderName, j,
                      currentConfigPath, TIMING, TEMPFOLDER)
        
      }
      
      print("Finish make alignments")
      
      # Put all the logs and all the configs toguether
      totalLogs    <- unifyFiles(resultsFolder,            "SUBLOG"             , "alignmentLog.txt"        ,FALSE,TRUE,FALSE)
      totalTemp    <- unifyFiles(resultsFolder,            "TEMPORALFILES.txt"  , "temporal_files.txt"      ,FALSE,TRUE,FALSE)
      totalConfigs <- unifyFiles(currentResultsFolderName, "configFile_results" , "config_results.txt"      ,TRUE, TRUE,FALSE)
      totalConTime <- unifyFiles(currentResultsFolderName, "configFile_Timing"  , "configFile_Timing.txt"   ,TRUE, TRUE,FALSE)
      totalBarcode <- unifyFiles(currentResultsFolderName, "barcodeFile_results", "barcodeFile_results.txt" ,TRUE, TRUE,FALSE)
      totalBarTime <- unifyFiles(currentResultsFolderName, "barcodeFile_Timing" , "barcodeFile_Timing.txt"  ,TRUE, TRUE,FALSE)
      
    }
  
    # Close the master log file descriptor.
    close(logFileConn)
    
    # Close the temporal files list file descriptor
#     close(tempFileConn)
    
    # Stop the cluster.
    stopCluster(cl)
    
    # If the user want to delete the uncompressed results, do it now.
    if(DELETEFQ == TRUE){
      print("Total temporal files deleted")
      #print(paste(resultsFolder,"/temporal_files.txt",sep=''))
      print(deleteFiles(paste(resultsFolder,"/temporal_files.txt",sep=''), deleteSource = FALSE))
    }
    
    print("Finish!")
    print(paste("Look for your results in",resultsFolder))
  
    t2 <- Sys.time()
    td <- as.numeric(t2-currentTime, units = "secs")
    
    print(paste("Total time: ",td * TOTAL_PROCESSORS ))

  }
  
}
