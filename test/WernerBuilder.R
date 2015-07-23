# Werner Builder

FWDSEQUENCES <- "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.1.fastq"
FWDBARCODES  <- "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.barcode_1.fastq"

wenerBuildSequence <- function(readsFilePath){

  unzipReadsFileName <- NULL
  
  # Check that the fileName has .gz in it.
  # If it does we need to unzip it. Otherwise we can start reading it directly.
  
  if(length(grep(".gz", readsFilePath, ignore.case = TRUE))==1){
    
    print(paste("Unzipping",readsFilePath))
    
    # Unzip the file
    gunzip(readsFilePath, skip = TRUE, overwrite = FALSE, remove = FALSE) #If it is already unzipped, use that
    unzipReadsFileName <- sub(".gz", "", readsFilePath, fixed=TRUE)
    
  }
  else{
    
    #print(paste("The file",fileName,"do not need to be unzipped"))
    unzipReadsFileName <- readsFilePath
    
  }
  
  print("A.1!")
  
  # Print feedback of the reading for the user, very large files can take a while to read
  print(paste("Reading file of size",round(file.info(unzipReadsFileName)[["size"]]/1048576,2),"MiBs"))
  
  # Read the file into a string variable
  textFile <- readLines(toString(unzipReadsFileName),encoding="UTF-8")
  
  # The number of lines must be multiple of 4 ALWAYS!
  if(length(textFile) %% 4 != 0){
    print(fileName)
    # TODO: Don't stop this, just skip whole line.
    stop("INCORRECT FILE: This file has a number of lines that is not multiple of 4!")
  }
  
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
  matrixTable <- matrix(ncol=14, nrow=totalRowsTable)
  
  # Fill the matrix with forward reads
  #print("Filling matrix with reads")    
  
  z <- 0
  pBar <- txtProgressBar(style=3)
  for (k in 1:length(textFile)){
    
    # Move forward the progress bar
    setTxtProgressBar(pBar, k/length(textFile))
    
    # Check if we are in the METADATA, SEQUENCE, + , QUALITY line
    z <- k %% 4
    
    # ID line
    if(z==1){
      
      # Separate fist by space, the left side is the Instrument until Y-coordinate (7 values), the right side
      # is from member until index (4 values); everything separated by ":"
      leftRight <- unlist(strsplit(textFile[k], " ", fixed = TRUE))
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
      
      matrixTable[ceiling(k/4),1] <- instrument
      matrixTable[ceiling(k/4),2] <- run
      matrixTable[ceiling(k/4),3] <- flowcellID
      matrixTable[ceiling(k/4),4] <- flowcellLine
      matrixTable[ceiling(k/4),5] <- tile
      matrixTable[ceiling(k/4),6] <- x
      matrixTable[ceiling(k/4),7] <- y
      matrixTable[ceiling(k/4),8] <- member
      matrixTable[ceiling(k/4),9] <- filtered
      matrixTable[ceiling(k/4),10] <- control
      matrixTable[ceiling(k/4),11] <- index
      
    }
    
    # Sequence line
    if(z==2){
      matrixTable[ceiling(k/4),12] <- textFile[k]
    }
    
    # + line
    if(z==3){
      matrixTable[ceiling(k/4),13] <- textFile[k]
    }
    
    # Quality line
    if(z==0){
      matrixTable[ceiling(k/4),14] <- textFile[k]
    }
    
  }
  
  # Transform the matrix into a dataframe
  #TODO: Check by making the dataframe from the begining
  table.df <- data.frame(matrixTable)  
  colnames(table.df) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Member","Filtered","Control","Index","Sequence","Extra","Quality")  
  
  return (table.df)

}

wenerBuildBarcodes <- function(barcodeFilePath){
  
  unzipReadsFileName <- NULL
  
  # Check that the fileName has .gz in it.
  # If it does we need to unzip it. Otherwise we can start reading it directly.
  
  if(length(grep(".gz", barcodeFilePath, ignore.case = TRUE))==1){
    
    print(paste("Unzipping",barcodeFilePath))
    
    # Unzip the file
    gunzip(barcodeFilePath, skip = TRUE, overwrite = FALSE, remove = FALSE) #If it is already unzipped, use that
    unzipReadsFileName <- sub(".gz", "", barcodeFilePath, fixed=TRUE)
    
  }
  else{
    
    #print(paste("The file",fileName,"do not need to be unzipped"))
    unzipReadsFileName <- barcodeFilePath
    
  }
  
  print("B.1!")
  
  # Print feedback of the reading for the user, very large files can take a while to read
  print(paste("Reading file of size",round(file.info(unzipReadsFileName)[["size"]]/1048576,2),"MiBs"))
  
  # Read the file into a string variable
  textFile <- readLines(toString(unzipReadsFileName),encoding="UTF-8")
  
  # The number of lines must be multiple of 4 ALWAYS!
  if(length(textFile) %% 4 != 0){
    print(fileName)
    # TODO: Don't stop this, just skip whole line.
    stop("INCORRECT FILE: This file has a number of lines that is not multiple of 4!")
  }
  
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
  # [12] Sequence line (Barcode)
  # [13] + line
  # [14] Quality
  matrixTable <- matrix(ncol=14, nrow=totalRowsTable)
  
  # Fill the matrix with forward reads
  #print("Filling matrix with reads")    
  
  z <- 0
  pBar <- txtProgressBar(style=3)
  for (k in 1:length(textFile)){
    
    # Move forward the progress bar
    setTxtProgressBar(pBar, k/length(textFile))
    
    # Check if we are in the METADATA, SEQUENCE, + , QUALITY line
    z <- k %% 4
    
    # ID line
    if(z==1){
      
      # Separate fist by space, the left side is the Instrument until Y-coordinate (7 values), the right side
      # is from member until index (4 values); everything separated by ":"
      leftRight <- unlist(strsplit(textFile[k], " ", fixed = TRUE))
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
      
      matrixTable[ceiling(k/4),1] <- instrument
      matrixTable[ceiling(k/4),2] <- run
      matrixTable[ceiling(k/4),3] <- flowcellID
      matrixTable[ceiling(k/4),4] <- flowcellLine
      matrixTable[ceiling(k/4),5] <- tile
      matrixTable[ceiling(k/4),6] <- x
      matrixTable[ceiling(k/4),7] <- y
      matrixTable[ceiling(k/4),8] <- member
      matrixTable[ceiling(k/4),9] <- filtered
      matrixTable[ceiling(k/4),10] <- control
      matrixTable[ceiling(k/4),11] <- index
      
    }
    
    # Sequence line
    if(z==2){
      matrixTable[ceiling(k/4),12] <- textFile[k]
    }
    
    # + line
    if(z==3){
      matrixTable[ceiling(k/4),13] <- textFile[k]
    }
    
    # Quality line
    if(z==0){
      matrixTable[ceiling(k/4),14] <- textFile[k]
    }
    
  }
  
  # Transform the matrix into a dataframe
  #TODO: Check by making the dataframe from the begining
  table.df <- data.frame(matrixTable)  
  colnames(table.df) <- c("Instrument_ID","Run_ID","Flowcell_ID","Flowcell_Line","Tile","X","Y","Member","Filtered","Control","Index","Barcode","Extra","Quality")  
  
  return (table.df)
  
}

wernerBuildReadFile <- function(readsFilePath, barcodeFilePath) {
 
  print("B!")
  
  # Now we read the file with the barcodes
  barcodesDataframe <- wenerBuildBarcodes(barcodeFilePath)
  
  print("A!")
  
  # First we read the file with the forward/reverse reads
  sequencesDataframe <- wenerBuildSequence(readsFilePath)
  

  
  # Lets create a folder where the data is going to be contained
  resultsFolder <- paste(getwd(), "/wernerFiles/", sep = '')
  dir.create(file.path(resultsFolder), showWarnings = TRUE)
  
  # Now, we add the barcode column to the sequences dataframe (maybe we should force a inner join?)
  sequencesDataframe$Barcode <- barcodesDataframe$Barcode
}


wernerBuildReadFile(FWDSEQUENCES,FWDBARCODES)