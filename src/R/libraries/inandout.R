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
#   If the file is zipped, you must have write rights.
#
# The function takes the following parameters:
#   
#   fileName: (String) A String with the absolute path to the file.
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
#   - Line with the sequence
#   - '+ line'
#   - Line with the quality of the nucleotides
# 
#   This means that the number of lines in that file must be divisible by 4.
# 
#   Worth mentioning that this function is suppose to be use for the forward
#   and reverse reads. The number of lines in each of those FASTQ files must
#   be the same.

}
getReadsFileDELETE <- function(fileName) {
  
  unzipReadsFileName <- NULL
  table.df           <- NULL
  
  # Check that the fileName has .gz in it.
  # If it does we need to unzip it. Otherwise we can start reading it directly.
  
  if(length(grep(".gz", fileName, ignore.case = TRUE))==1){
    
    # Unzip the file
    gunzip(fileName, skip = TRUE, overwrite = FALSE, remove = FALSE) #If it is already unzipped, use that
    unzipReadsFileName <- sub(".gz", "", fileName, fixed=TRUE)
    
  }
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
        matrixTable[ceiling(k/4),1] <- textFile[k]
      }
      
      # Quality line
      if(z==0){
        matrixTable[ceiling(k/4),2] <- textFile[k]
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