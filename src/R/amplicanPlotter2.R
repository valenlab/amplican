############################################
# Constant declarations
############################################
{

  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/res/Alignment_Toy/Analysis_Test"
  #ANALYSIS_PATH <- "/scratch/ampliCan/toySet/Alignment_Toy/Analysis_Test"
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150424215630_run11/Analysis_3_5_TRUE_TRUE_10_20150429113211"
  ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150511142604_run11/Analysis_3_5_TRUE_TRUE_10_20150511142938"
  
  
  # Jason
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150410175748_config_Jason2/Analysis_3_5_TRUE_TRUE_10_20150410181017"
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150410180646_JasonAfterChristmas/Analysis_3_5_TRUE_TRUE_10_20150410181240"
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150506233015_Jason3/Analysis_3_5_TRUE_TRUE_10_20150507200426"
  #ANALYSIS_PATH  <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150508172357_Jason4/Analysis_3_5_TRUE_TRUE_10_20150511150106"
  
  
  # Danielle
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150320154239_20141112_F/Analysis_3_5_TRUE_TRUE_10_20150320222005"
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150321043346_20141112_R/Analysis_3_5_TRUE_TRUE_10_20150324104553"
  #ANALYSIS_PATH <- "/Home/ii/rafaelc/Desktop/git/ampliCan/src/R/results/run_20150321043622_20141117_F2/Analysis_3_5_TRUE_TRUE_10_20150324062211"
  
  TOTAL_PROCESSORS <- 1       # Set this to the number of processors you want
                               # to use. You can set it to NULL if you want all
                               # that are available.
  
  SELECTED_THEME <- 1 # Which Theme do you want
                      # 1: Default
                      # 2: Colorblind
                      # 3: Greyscale
                      # 4: User define (given by argument on the function call)
  
  # Colors for the cut rates plots
  CUT_RATE    <- c("#00bfc4","#66CC33", "#AAAAAA")                      # BLUE
  NO_CUT_RATE <- c("#f8766d","#000099", "#333333")  # 23DB19 FF331D ?   # RED
  
  # Colors for the frameshift plot
  NO_DELETION <- c("#999999","#000099", "#111111") # GREY    c("#666666","#666666") "#999999"
  IN_FRAME    <- c("#f5ff70","#000000", "#777777") # ORANGE  c("#f8766d","#66CC33") "#FF9966"
  FRAMESHIFT  <- c("#70ff9d","#66CC33", "#FFFFFF") # GREEN   c("#7cae00","#000099") "#7cae00"
  MULTIPLES   <- c("#ff70ca","#FFFFFF", "#CCCCCC") # MAGENTA c("#619cff","#FFFFFF") "#CC66FF"
  
  # Colors for the archplots
  DELETION_COLOR   <- c("#CC0000","#000099", "#333333") # RED
  CUT_COLOR        <- c("#66CC33","#66CC33", "#AAAAAA") # BLUE
  POSITION_WINDOW  <- c("#3399FF","#FFFFFF", "#000000")
  START            <- c("#00FF00","#9900CC", "#000000")
  END              <- c("#00FF00","#9900CC", "#000000")
  
  # The archplot can be all be normalized between -1 to 1, or be "zoom in" between their respective
  # maximum and minimum frequency (which will be more closer to 0 than -1 and 1). If you want it
  # to be between -1 and 1 always, set this variable to TRUE
  FREQONE <- TRUE
  
  # The arch can represent one sequence each, or we can accumulate all the events that start and end
  # at the same interval. If you want to accumulate the height of the archs, set this variable to TRUE.
  # But if you want each arch to represent each sequence, set this to FALSE.
  ACCUMULATE <- TRUE
  
  # You can filter start and end deletions. Usually you will get a lot of deletions that start at 1
  # and ends at the end of the amplicon, that are there because the sequences are just short and that part
  # don't alignt with anything. If you want to do that set this variable to TRUE. Otherwise set it to FALSE.
  # This will skip all deletions (not cuts) that start at 1 and all deletions that end at length(amplicon)
  FILTERENDS <- TRUE
  
  # Colors for the barcodes plots
  NOASSIGNED <- c("#999999","#000000", "#111111") # GREY
  SEQUENCE   <- c("#f8766d","#000099", "#333333") # RED  c("#527FDF","#66CC33") 
  READ       <- c("#00bfc4","#66CC33", "#AAAAAA") # BLUE c("#7FDF52","#FFFFFF")
  PASS       <- c("#7FDF52","#FFFFFF", "#FFFFFF") # GREEN
  FAILED     <- c("#000000","#000000", "#000000") # BLACK c("#DF527F","#000000")
  
  # Theme for the plots
  PLOTBACK <- "#FFFFFF"
  PLOTGRID <- "#888888"
  LINETYPE <- "dotted"
  
  DUPLICATE_INTO_FOLDERS <- FALSE # The plots are going to be group into a single folder. Beside that,
                                  # you can copy each of the archplots and the histogram plots into
                                  # each of the individual folders that compose the analysis
  FORMAT <- "png"
  
  GROUPBY <- 20 # If we have a lot of rows we will plot giagantic plots
                # instead, group by this number and make many tiny plots for cutRates, and so on
  
  # Archplots and mutations plots are images which width is desproportional to anything else.
  # The image size is quite complex as it involve several variables to draw nicelly text, bars, and many things
  # Please use one of this preselected sizes:
  # 4K   - 11811 x 2100 - 1000x178 mm²
  # 1080 - x 1080
  # 720  - x 720
  # VGA  - x 480
  SIZE   <- "4K" #2100 x 11811
}


############################################
#Summons the libraries and get to use functions from there.
############################################
{
  
  # Utilities
  
  # Plotting stuff
  library(reshape2) # Rearrenge data for ggplot
  library(ggplot2)
  library(GenomicRanges) # For the archplots
  library(ggbio)

  
  # Making things in parallel
  library(doParallel)
  cl <- makeCluster(TOTAL_PROCESSORS, outfile="") #outfile tells the workers to print on the terminal
  registerDoParallel(cl)
  
  # Our own developed stuff.
#   source("libraries/inandout2.R") # How to read and write from the disk.
#   source("libraries/cutscriterias2.R") # How deletions are defined.
  source("libraries/tools2.R") # Minor stuff like the reverse complement of a DNA sequence.
#   source("libraries/errorswarnings2.R") # Handles pre-parsing, error messages, warning to the users, and so on.
  source("libraries/configPlotting.R") # Function that deals with the config files and get running everything.
#   source("libraries/recovery2.R")
  
}


############################################
# Algorith starts here
############################################

print("START")

# We are going to write the results in the same folder, inside another folder called:
# Plot_<plot variables>_YYYYmmddHHMMSS
{
  
  # Creates the folder for the plots
  currentTime <- Sys.time()
  timeStamp   <- strftime(currentTime,"%Y%m%d%H%M%S")
  plotPath    <- paste(ANALYSIS_PATH, "/Plots_",SELECTED_THEME,"_",FORMAT,"_",timeStamp, sep = '')
  
  dir.create(file.path(plotPath), showWarnings = TRUE)
  
  # Creates a folder for all the data that is shown in the plots
  dataPath   <- paste(plotPath,"/Data",sep='',collapse='')
  dir.create(file.path(dataPath), showWarnings = TRUE)
  
  # Inside this folder it will be a file named log where we will write relevant information
  logFileName <- paste(plotPath, "/PlotsLog.txt", sep = '')
  logFileConn <- file(logFileName, open="w")
  
  # Write the analysis variables into the log
  writeLines("The plot was made with the following variables: ", logFileConn)
  writeLines(" -- Cut Rate -- ", logFileConn)
  writeLines(paste("CUT_RATE           : ", CUT_RATE[SELECTED_THEME]), logFileConn)
  writeLines(paste("NO_CUT_RATE        : ", NO_CUT_RATE[SELECTED_THEME]), logFileConn)
  writeLines(" -- Arch plot -- ", logFileConn)
  writeLines(paste("DELETION_COLOR     : ", DELETION_COLOR[SELECTED_THEME]), logFileConn)
  writeLines(paste("CUT_COLOR          : ", CUT_COLOR[SELECTED_THEME]), logFileConn)
  writeLines(paste("POSITION_WINDOW    : ", POSITION_WINDOW[SELECTED_THEME]), logFileConn)
  writeLines(paste("START              : ", START[SELECTED_THEME]), logFileConn)
  writeLines(paste("END                : ", END[SELECTED_THEME]), logFileConn)
  writeLines(" -- Others -- ", logFileConn)
  writeLines(paste("DUPLICATE_INTO_FOLDERS : ", DUPLICATE_INTO_FOLDERS), logFileConn)  
  writeLines(paste("GROUP BY               : ", GROUPBY), logFileConn)  
  writeLines("", logFileConn) # Empty line
  
  # Find out where is the alignments, we need to get some data from there.
  analysisTree        <- strsplit(ANALYSIS_PATH, "/", fixed = TRUE)[[1]]
  analysisTreeDeep    <- length(analysisTree)
  alignmentTree       <- analysisTree[1:analysisTreeDeep -1]
  ALIGNMENTFOLDER     <- paste(alignmentTree, collapse='/')
  ALIGNMENTFOLDER     <- list.files(path = ALIGNMENTFOLDER, pattern = "_alignments", full.names = TRUE)  
  UNNASIGNEDFILES     <- list.files(path = paste(ALIGNMENTFOLDER,"/unassigned_sequences",sep='',collapse=''), full.names = TRUE)
  
}

# Read the list of folders for this analysis. We are going to make several plots for each of them.

# This plots will be made on a single processor

# Read the config file from the analysis folder and get the ID+Barcode and Cut Rate
analysisDF <-  read.table(paste(ANALYSIS_PATH, "/analysis_results.txt", sep=''), header=TRUE, sep="\t")

# Read the config file from the alignment folder and get the barcode info

#----------------------------
# --- CUT RATE PLOT (Both normalize and absolute)
#----------------------------
{
 
  print("Plotting Cut Rates")
  
  keeps <- c("ID","Barcode", "Cut_Rate")
  cutRateDF <- analysisDF[keeps]
  cutRateDF <- cutRateDF[ order(as.numeric(row.names(cutRateDF))), ]
  cutRateDF$Complement <- 1 - cutRateDF$Cut_Rate
  
  cutRateDF$ID_Barcode <- paste(cutRateDF$ID,"_",cutRateDF$Barcode,sep='')
  
  subCutRatesDF   <- split(cutRateDF, factor(sort(rank(row.names(cutRateDF))%%ceiling(nrow(cutRateDF)/GROUPBY))))

  # We are plotting cut rates in groups of 20 (default)
  # If we happen to have a lot of rows we will create an enourmous image
  
  for (j in 1:length(subCutRatesDF)){
    
    cutRateDF <- subCutRatesDF[[j]][,c("ID_Barcode","Cut_Rate","Complement")]
    
    # If all the data for the complement and the Cut Rate is NA, the ggplot will fail.
    # We need to make sure that at least one is 0
    if(is.na(cutRateDF$Complement[1])){
      cutRateDF$Complement[1] <- 0  
    }
    
    colnames(cutRateDF) <- c("ID_Barcode", "Cut Reads", "Not Cut")
    
    cutRateDFMelted <- melt( cutRateDF , id.var="ID_Barcode")
    colnames(cutRateDFMelted) <- c("ID_Barcode","variable", "Relative")
    
    ggplot(cutRateDFMelted, aes(x = ID_Barcode, y = Relative, fill = variable)) +
      geom_bar(stat = "identity", colour="black") +
      scale_fill_manual(values=c(CUT_RATE[SELECTED_THEME], NO_CUT_RATE[SELECTED_THEME])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.background = element_rect(fill = PLOTBACK),
            panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
            panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    ggsave(paste(plotPath,"/Cut_Rates_",j,".", FORMAT ,sep = ''),width= nrow(subCutRatesDF[[j]]) * 20 , units="mm", limitsize=FALSE)
    
    
    fileName <- paste(dataPath,"/Cut_Rates_",j,".txt",sep = '', collapse = '')
    write.table(cutRateDFMelted, file = fileName, quote = FALSE, sep = "\t")
    
  }
  
  keeps <- c("ID", "Barcode", "Sum_Is_Sequence", "Sum_Is_Read")
  cutRateAbsoluteDF <- analysisDF[keeps]
  
  cutRateAbsoluteDF <- cutRateAbsoluteDF[ order(as.numeric(row.names(cutRateAbsoluteDF))), ]
  cutRateAbsoluteDF$Sum_Is_Sequence <- cutRateAbsoluteDF$Sum_Is_Sequence - cutRateAbsoluteDF$Sum_Is_Read
  cutRateAbsoluteDF$ID_Barcode <- paste(cutRateAbsoluteDF$ID,"_",cutRateAbsoluteDF$Barcode,sep='')
  subCutRatesAbsoluteDF   <- split(cutRateAbsoluteDF, factor(sort(rank(row.names(cutRateAbsoluteDF))%%ceiling(nrow(cutRateAbsoluteDF)/GROUPBY))))
  
    
  for (j in 1:length(subCutRatesAbsoluteDF)){
    
    cutRateAbsoluteDF <- subCutRatesAbsoluteDF[[j]][,c("ID_Barcode","Sum_Is_Read","Sum_Is_Sequence")]
    colnames(cutRateAbsoluteDF) <- c("ID_Barcode", "Cut Reads", "Not Cut")
    
    cutRateAbsoluteDFMelted <- melt( cutRateAbsoluteDF , id.var="ID_Barcode")
    colnames(cutRateAbsoluteDFMelted) <- c("ID_Barcode","variable", "Reads")
    
    ggplot(cutRateAbsoluteDFMelted, aes(x = ID_Barcode, y = Reads, fill = variable)) +
      geom_bar(stat = "identity", colour="black") +
      scale_fill_manual(values=c(CUT_RATE[SELECTED_THEME], NO_CUT_RATE[SELECTED_THEME])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.background = element_rect(fill = PLOTBACK),
            panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
            panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    ggsave(paste(plotPath,"/Cut_Rates_Absolute_",j,".", FORMAT ,sep = ''),width= nrow(subCutRatesDF[[j]]) * 20 , units="mm", limitsize=FALSE)
    
    fileName <- paste(dataPath,"/Cut_Rates_Absolute_",j,".txt",sep = '', collapse = '')
    write.table(cutRateAbsoluteDFMelted, file = fileName, quote = FALSE, sep = "\t")
    
  }
  
  
}

#----------------------------
# --- FRAMESHIFT PLOT
#----------------------------
{

  print("Plotting Frameshifts")
  
  # Now; lets make the in-frame and frameshift statistics.
  frameShiftTable <- analysisDF[,c("ID","Barcode",
                                   "Sum_Is_Sequence",
                                   "Forward_No_Deletions","Forward_In_Frameshift",
                                   "Forward_Out_Frameshift","Forward_Multiple_Deletions",
                                   "Sum_Is_Read",
                                   "Forward_Read_No_Deletions","Forward_Read_In_Frameshift",
                                   "Forward_Read_Out_Frameshift","Forward_Read_Multiple_Deletions")]
  
  # Order by row ID
  frameShiftTable <- frameShiftTable[ order(as.numeric(row.names(frameShiftTable))), ]
  
  # Get the percentage for the frameshift of ALL deletions
  frameShiftTable$Percentage_No_Deletion   <- frameShiftTable$Forward_No_Deletions / frameShiftTable$Sum_Is_Sequence
  frameShiftTable$Percentage_In_Frameshift <- frameShiftTable$Forward_In_Frameshift / frameShiftTable$Sum_Is_Sequence
  frameShiftTable$Percentage_Out_Frameshift <- frameShiftTable$Forward_Out_Frameshift / frameShiftTable$Sum_Is_Sequence
  frameShiftTable$Percentage_Multiple_Deletions <- frameShiftTable$Forward_Multiple_Deletions / frameShiftTable$Sum_Is_Sequence
  
  # Get the percentage for the frameshift of ONLY valid deletions
  frameShiftTable$Percentage_Reads_No_Deletion   <- frameShiftTable$Forward_Read_No_Deletions / frameShiftTable$Sum_Is_Read
  frameShiftTable$Percentage_Reads_In_Frameshift <- frameShiftTable$Forward_Read_In_Frameshift / frameShiftTable$Sum_Is_Read
  frameShiftTable$Percentage_Reads_Out_Frameshift <- frameShiftTable$Forward_Read_Out_Frameshift / frameShiftTable$Sum_Is_Read
  frameShiftTable$Percentage_Reads_Multiple_Deletions <- frameShiftTable$Forward_Read_Multiple_Deletions / frameShiftTable$Sum_Is_Read
  
  frameShiftTable$ID_Barcode <- paste(frameShiftTable$ID,"_",frameShiftTable$Barcode,sep='')
  
  subframeShiftDF <- split(frameShiftTable, factor(sort(rank(row.names(frameShiftTable))%%ceiling(nrow(frameShiftTable)/GROUPBY))))
  
  # Once more, we try to avoid very wide bar plots
  for (j in 1:length(subframeShiftDF)){
    
    # FRAMESHIFT FOR ALL THE SEQUENCES
    # ---------------------------------
    {
      
    frameshiftDF <- subframeShiftDF[[j]][,c("ID_Barcode","Percentage_No_Deletion","Percentage_In_Frameshift","Percentage_Out_Frameshift","Percentage_Multiple_Deletions")]
    
    # If all the data for the percentages is NA or NaN
    # We need to make sure that at least one is 0
    if(is.nan(frameshiftDF$Percentage_No_Deletion[1])){
      frameshiftDF$Percentage_No_Deletion[1]        <- 0
      frameshiftDF$Percentage_In_Frameshift[1]      <- 0
      frameshiftDF$Percentage_Out_Frameshift[1]     <- 0
      frameshiftDF$Percentage_Multiple_Deletions[1] <- 0
    }
    
    colnames(frameshiftDF) <- c("ID_Barcode","No Deletion","Preserve Frame","Frameshift","Multiple Deletions")
    
    # Reorder so the frameshift is at the botoom of the plot
    frameshiftDF <- frameshiftDF[c("ID_Barcode", "Frameshift", "Preserve Frame", "Multiple Deletions", "No Deletion" )]
    
    frameShiftTableMelted <- melt(frameshiftDF, id.var="ID_Barcode")
    colnames(frameShiftTableMelted) <- c("ID_Barcode","variable", "Frequency_Relative")
    
    #labels <- frameShiftTable$ID_Barcode
    ggplot(frameShiftTableMelted, aes(x = ID_Barcode, y = Frequency_Relative, fill = variable)) +
          geom_bar(stat = "identity", colour="black") +
          scale_fill_manual(values=c(FRAMESHIFT[SELECTED_THEME], IN_FRAME[SELECTED_THEME], MULTIPLES[SELECTED_THEME], NO_DELETION[SELECTED_THEME] )) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right", 
                panel.background = element_rect(fill = PLOTBACK),
                panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    ggsave(paste(plotPath,"/Frameshift_",j,".", FORMAT ,sep = ''),
           width= nrow(subframeShiftDF[[j]]) * 20 ,
           units="mm", limitsize=FALSE)
    
    fileName <- paste(dataPath,"/Frameshift_",j,".txt",sep = '', collapse = '')
    write.table(frameShiftTableMelted, file = fileName, quote = FALSE, sep = "\t")
    
    
    }
    
    # FRAMESHIFT ONLY FOR READS
    # ---------------------------------
    {
    
      frameshiftReadsDF <- subframeShiftDF[[j]][,c("ID_Barcode","Percentage_Reads_No_Deletion","Percentage_Reads_In_Frameshift","Percentage_Reads_Out_Frameshift","Percentage_Reads_Multiple_Deletions")]
      
      # If all the data for the percentages is NA or NaN
      # We need to make sure that at least one is 0
      if(is.nan(frameshiftReadsDF$Percentage_Reads_No_Deletion[1])){
        frameshiftReadsDF$Percentage_Reads_No_Deletion[1]        <- 0
        frameshiftReadsDF$Percentage_Reads_In_Frameshift[1]      <- 0
        frameshiftReadsDF$Percentage_Reads_Out_Frameshift[1]     <- 0
        frameshiftReadsDF$Percentage_Reads_Multiple_Deletions[1] <- 0
      }
      
      
      # Change columns names so they are more recognizable
      colnames(frameshiftReadsDF) <- c("ID_Barcode", "No Deletion", "Preserve Frame", "Frameshift", "Multiple Deletions")
      
      # Reorder so the frameshift is at the botoom of the plot
      frameshiftReadsDF <- frameshiftReadsDF[c("ID_Barcode", "Frameshift", "Preserve Frame", "Multiple Deletions", "No Deletion")]
      
      frameShiftTableMelted <- melt(frameshiftReadsDF, id.var="ID_Barcode")
      colnames(frameShiftTableMelted) <- c("ID_Barcode","variable", "Frequency_Relative")
      
      #labels <- frameShiftTable$ID_Barcode
      ggplot(frameShiftTableMelted, aes(x = ID_Barcode, y = Frequency_Relative, fill = variable)) +
        geom_bar(stat = "identity", colour="black") +
        scale_fill_manual(values=c(FRAMESHIFT[SELECTED_THEME], IN_FRAME[SELECTED_THEME], MULTIPLES[SELECTED_THEME], NO_DELETION[SELECTED_THEME])) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right",
              panel.background = element_rect(fill = PLOTBACK),
              panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
              panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
      
      ggsave(paste(plotPath,"/Frameshift_ReadsOnly_",j,".", FORMAT ,sep = ''),
             width= nrow(subframeShiftDF[[j]]) * 20 ,
             units="mm", limitsize=FALSE)
      
      fileName <- paste(dataPath,"/Frameshift_ReadsOnly",j,".txt",sep = '', collapse = '')
      write.table(frameShiftTableMelted, file = fileName, quote = FALSE, sep = "\t")
      
    }
      
    
    
  }

}

#----------------------------
# --- BARCODE PLOT
#----------------------------
{
  
  print("Plotting Barcodes Info")
  
  #Make the barcode dataframe
  barcodeDF <-  read.table(paste(ALIGNMENTFOLDER, "/barcodeFile_results.txt", sep=''), header=TRUE, sep="\t")
  uniqueBarcodes <- barcodeDF$Barcode
  
  # -- Prepare the final dataframe where everything goes
  barcodeTableFinal=data.frame(matrix(NA, nrow=length(uniqueBarcodes), ncol=9))
  names(barcodeTableFinal) <- c("Barcode","Total_Pre_N_Filter","Total_Post_N_Filter",
                                "Failed_N_Filter","Experiment_Unique_Sequences",
                                "Barcode_Total_Sequences",
                                "Barcode_Total_Reads",
                                "Barcode_Unnasigned","Barcode_non_read")
  
  # -- For each of those barcodes
  for(j in 1:length(uniqueBarcodes)){
            
    #Make a subset with that barcode (keep the barcode row)
    barcodeSubset <- subset(barcodeDF, Barcode == uniqueBarcodes[j])
    
    #Make a subset with that barcode ID from the analysis (might be several)
    analysisSubset <- subset(analysisDF, Barcode == uniqueBarcodes[j])
            
    #Update the rest of the columns in the new dataframe
    barcodeTableFinal$Barcode[j] <- toString(uniqueBarcodes[j])
    barcodeTableFinal$Total_Pre_N_Filter[j] <- barcodeSubset$Total_Pre_N_Filter[1]
    barcodeTableFinal$Total_Post_N_Filter[j] <- barcodeSubset$Total_Post_N_Filter[1]
    barcodeTableFinal$Experiment_Unique_Sequences[j] <- barcodeSubset$Experiment_Unique_Sequences[1]
    
    barcodeTableFinal$Barcode_Total_Sequences[j] <- sum(analysisSubset$Sum_Is_Sequence)
    barcodeTableFinal$Barcode_Total_Reads[j]     <- sum(analysisSubset$Sum_Is_Read)
            
    # These columns are just a combination of the previous ones, but we need it for melting the data later
    barcodeTableFinal$Failed_N_Filter[j]    <- barcodeSubset$Total_Pre_N_Filter[1]  - barcodeSubset$Total_Post_N_Filter[1]
    barcodeTableFinal$Barcode_Unnasigned[j] <- barcodeSubset$Total_Post_N_Filter[1] - barcodeTableFinal$Barcode_Total_Sequences[j]
    barcodeTableFinal$Barcode_non_read[j]   <- barcodeTableFinal$Barcode_Total_Sequences[j] - barcodeTableFinal$Barcode_Total_Reads[j]
            
            
  }
          
  # Get the data we need for each of the plot
  barcodeAssignedData <- subset(barcodeTableFinal, select=c(Barcode,Barcode_Unnasigned,Barcode_non_read,Barcode_Total_Reads))
  barcodeFilterData   <- subset(barcodeTableFinal, select=c(Barcode,Failed_N_Filter,Total_Post_N_Filter))

  barcodeAssignedData <- barcodeAssignedData[ order(as.numeric(row.names(barcodeAssignedData))), ]
  barcodeFilterData   <- barcodeFilterData[   order(as.numeric(row.names(barcodeFilterData)))  , ]
  
  # Change the names of the columns for human users
  colnames(barcodeAssignedData) <- c("Barcode", "Unnasigned", "No cut", "Cut Read")
  colnames(barcodeFilterData)   <- c("Barcode", "Failed"    , "Pass")
  
  # Reorder so the important things are at the bottom of the plot  
  barcodeAssignedData <- barcodeAssignedData[c("Barcode", "Cut Read", "No cut", "Unnasigned")]
  barcodeFilterData   <- barcodeFilterData[  c("Barcode", "Pass"    , "Failed")]
  
  # Divide for the assigned barcodes
  subBarcodeAssignedDF <- split(barcodeAssignedData, factor(sort(rank(as.numeric(row.names(barcodeAssignedData)))%%ceiling(nrow(barcodeAssignedData)/GROUPBY))))
  
  for (j in 1:length(subBarcodeAssignedDF)){
    
    barcodeAssignedMelted <- melt(subBarcodeAssignedDF[[j]], id.var="Barcode")
    
    ggplot(barcodeAssignedMelted, aes(x = Barcode, y = value, fill = variable)) +
      geom_bar(stat = "identity", colour="black") +
      scale_fill_manual(values=c(READ[SELECTED_THEME], SEQUENCE[SELECTED_THEME], NOASSIGNED[SELECTED_THEME])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.background = element_rect(fill = PLOTBACK),
            panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
            panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    ggsave(paste(plotPath,"/barcode_assigned_",j,".",FORMAT ,sep = ''), width= (nrow(subBarcodeAssignedDF[[j]]) * 20) + 50 , units="mm", limitsize=FALSE)
    
    fileName <- paste(dataPath,"/barcode_assigned_",j,".txt",sep = '', collapse = '')
    write.table(barcodeAssignedMelted, file = fileName, quote = FALSE, sep = "\t")
    
  }
  
  # Divide for the filters
  subBarcodeFiltersDF <- split(barcodeFilterData, factor(sort(rank(as.numeric(row.names(barcodeFilterData)))%%ceiling(nrow(barcodeFilterData)/GROUPBY))))
  
  for (j in 1:length(subBarcodeAssignedDF)){
    
    barcodeFilterMelted <- melt(subBarcodeFiltersDF[[j]], id.var="Barcode")
    
    ggplot(barcodeFilterMelted, aes(x = Barcode, y = value, fill = variable)) +
      geom_bar(stat = "identity", colour="black") +
      scale_fill_manual(values=c(PASS[SELECTED_THEME], FAILED[SELECTED_THEME])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),                panel.background = element_rect(fill = PLOTBACK),
            panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
            panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    ggsave(paste(plotPath,"/barcode_filters_",j,".",FORMAT ,sep = ''),
           width= (nrow(subBarcodeAssignedDF[[j]]) * 20) + 50 ,
           units="mm", limitsize=FALSE)
    
    fileName <- paste(dataPath,"/barcode_filters_",j,".txt",sep = '', collapse = '')
    write.table(barcodeFilterMelted, file = fileName, quote = FALSE, sep = "\t")
    
  }
  
  # Now make the unassigned histograms for each barcode
  for (j in 1:length(UNNASIGNEDFILES)){
    
    # Get the barcode name
    barcode <- strsplit(strsplit(UNNASIGNEDFILES[j], "/", fixed = TRUE)[[1]][length(strsplit(UNNASIGNEDFILES[j], "/", fixed = TRUE)[[1]])], "_", fixed = TRUE)[[1]][1]
    
    # Get the data
    unnasignedDF <- read.table(UNNASIGNEDFILES[j])
    
    # Sort by Total
    unnasignedDF <- unnasignedDF[with(unnasignedDF, order(-Total)), ]
    
    # Get the top 15
    unnasignedDF <- unnasignedDF[1:15,]
    
    # Find out the row ID and add it to the leyend of each column
    baseZero <- nchar(toString(nrow(unnasignedDF)))
    unnasignedDF$ID <- ""
    
    for(p in 1:nrow(unnasignedDF)){
      
      unnasignedDF$ID[p] <- rownames(unnasignedDF[p,])
      
      #Find how many zeros we need to fix the sort
      totalZeros  <- baseZero - nchar(toString(p))
      stringZeros <- paste(rep("0",totalZeros),collapse='')
      
      unnasignedDF$ID[p] <- paste(stringZeros,p,"_*_",toString(unnasignedDF$ID[p]),sep='')
      
    }
    
    # Make the plot and save it
    
    unnasignedPlot <- ggplot(data=unnasignedDF, aes(x=ID, y=Total), na.rm = TRUE)
    unnasignedPlot <- unnasignedPlot + geom_bar(stat="identity", na.rm = TRUE, fill = "black")          
    unnasignedPlot <- unnasignedPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                           panel.background = element_rect(fill = PLOTBACK),
                                           panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                                           panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                                           panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    
    # Save the file into the plot folder
    ggsave(paste(plotPath,"/",barcode,"_unnasigned.",FORMAT ,sep = ''))
    
    fileName <- paste(dataPath,"/",barcode,"_unnasigned.txt",sep = '', collapse = '')
    write.table(unnasignedDF, file = fileName, quote = FALSE, sep = "\t")
    
    
  }
  
  
}



#-------------------------------------
# --- HISTOGRAM AND ARCHPLOT FOR EACH
#-------------------------------------
{
  
  print("Plotting Archplots")
  
  # Divide the config dataframe into several subset of roughly equal size
  processorsSubDataframes <- divideWork(TOTAL_PROCESSORS, analysisDF)
  
  # Analyze the config file in parallel 
  parallelPackages = c("ggplot2","GenomicRanges", "ggbio", "plyr", "reshape2")
  
  foreach(j=1:TOTAL_PROCESSORS , .packages=parallelPackages) %dopar% {
    
    
    makePlots(CUT_RATE[SELECTED_THEME], NO_CUT_RATE[SELECTED_THEME], NO_DELETION[SELECTED_THEME],
              IN_FRAME[SELECTED_THEME], FRAMESHIFT[SELECTED_THEME], MULTIPLES[SELECTED_THEME],
              DELETION_COLOR[SELECTED_THEME], CUT_COLOR[SELECTED_THEME], POSITION_WINDOW[SELECTED_THEME], START[SELECTED_THEME],
              END[SELECTED_THEME], DUPLICATE_INTO_FOLDERS,
              PLOTBACK, PLOTGRID, LINETYPE, FREQONE, ACCUMULATE, FILTERENDS, SIZE,
              processorsSubDataframes[[j]], j, ANALYSIS_PATH, plotPath, dataPath)
  
    
  }
}

# Clsoe the master log file descriptor
close(logFileConn)

# Stop the cluster and go home
stopCluster(cl)

print("Finish!")

print(as.numeric(Sys.time()-currentTime, units = "secs"))

print(paste("Look for your results in",plotPath))