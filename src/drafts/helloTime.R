#This is some time testing stuff:

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
#   cl <- makeCluster(total_processors, outfile="") #outfile tells the workers to print on the terminal
#   registerDoParallel(cl)
}

#Make some variables where to put the times
matrixTiming <- NULL
matrixToy <- NULL
t1 <- Sys.time()
t2 <- Sys.time()
t3 <- Sys.time()
t4 <- Sys.time()
td <- difftime(t2,t1)
td2 <- difftime(t4,t3)


#Some constant variables

sequence_size = 200          # This is the size of the sequences in the first example

total_mutations = 50         # This is the total of mutations you want to introduce.
                             # A mutation can be a deletion, so your second sequence will have size equal
                             # or less than the first one.

                             # Size SeqA >= Size SeqB >= Size SeqA - total mutations


############################################
# Align three random sequences many times and time all the things, even the times
############################################
pBar <- txtProgressBar(style=3)
for(i in 1:100){

  t3 <- Sys.time()
  {
  
  setTxtProgressBar(pBar, i/100)
    
  t1 <- Sys.time()
  #Make a random genome
  genomeString <- paste(sample(DNA_ALPHABET[1:4], sequence_size, replace = TRUE), collapse = "")
  genomeStringReversed <- paste(rev(substring(genomeString,1:nchar(genomeString),1:nchar(genomeString))),collapse="") 
  genome <- DNAStringSet(genomeString)
  genomeReversed <- DNAStringSet(genomeStringReversed)

  #Make a random reverse Sequence
  positions <- sample(sequence_size, total_mutations) #Take total mutations amount of numbers, from 0 to SIZE, no repetitions.
  ranges <- IRanges(positions, width=1) #Get them into a table form.
  DNASet <- c(DNA_ALPHABET[1:4], "") #Get the DNA characters and add ""(nothing) into an array of characters.
  newNucleotides <- sample(DNASet, total_mutations, replace = TRUE) #Get a totally new sequence of size total mutations.
  candidateReverseSequence <- replaceAt(genomeReversed,at=ranges,newNucleotides) #Replace 

  #Make a random forward Sequence
  positions <- sample(sequence_size, total_mutations) #Take total mutations amount of numbers, from 0 to SIZE, no repetitions.
  ranges <- IRanges(positions, width=1) #Get them into a table form.
  DNASet <- c(DNA_ALPHABET[1:4], "") #Get the DNA characters and add ""(nothing) into an array of characters.
  newNucleotides <- sample(DNASet, total_mutations, replace = TRUE) #Get a totally new sequence of size total mutations.
  candidateForwardSequence <- replaceAt(genome,at=ranges,newNucleotides) #Replace 
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Make_random_sequences",td))

  t1 <- Sys.time()
  for(j in 1:10){
    #Find the reverse complementary of the reverse sequence for the multiple alignment
    candidateComplementarySequence <- reverseComplement(DNAStringSet(candidateReverseSequence))
    
    #Second method
    #Get the forward and reverse candidate and make a DNASet with them
    #forrevSequences <- DNAStringSet(c(toString(candidateComplementarySequence), candidateForwardSequence))
    forrevSequences <- DNAStringSet(c(toString(candidateComplementarySequence)[1], toString(candidateForwardSequence)))
    
    #That DNA set is pairwise aligned with the genome
    myPairwise <- pairwiseAlignment(forrevSequences, genome)
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Pairwise_FandR",td))

  t1 <- Sys.time()
  #Do some random timing and time it
  for(j in 1:100){
    t4 <- Sys.time()  
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Silly_Timing",td))
  
  #Read a file many times
  t1 <- Sys.time()
  for(j in 1:1){
    textForwardFile <- readLines("1_S1_L001_R1_001.fastq",encoding="UTF-8")
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Reading_FASTQ_File",td))
  
  #Create a random dataframe
  t1 <- Sys.time()
  randomDataframe <- data.frame(replicate(100,sample(1:1000,1000,rep=TRUE)))
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Random_Dataframe",td))  
  
  #Get the fake ID and barcode
  t1 <- Sys.time()
  for(j in 1:100){
    currentID = randomDataframe[j,j]
    currentBarcode = randomDataframe[j,j]
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Get_info",td))
  
  #Progress feedback for the user
  t1 <- Sys.time()
  for(j in 1:10){
    print("________________________________________")
    print(paste("Working on",currentID, currentBarcode))
    print(paste("Line",j,"of",j,"% for",j))
    print("----------------------------------------")
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Printing",td))
  
  
  #Do some rounding
  t1 <- Sys.time()
  for(j in 1:100){
    temp <- round(j/100,2)
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Rounding",td))
  
  
  #Flow control test
  t1 <- Sys.time()
  for(k in 1:100){
    for(j in 1:100){
      if(randomDataframe[j,j]<500){
        temp <- 0
      }
      else{
        temp <- 1
      }
    }
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Flow_control",td))
  
  
  #Pasting strings
  t1 <- Sys.time()
  for(j in 1:1000){
    
    paste(randomDataframe[j,j],"Hello",randomDataframe[j,j],"Paste")
    
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Pasting_strings",td))
  
  #Create some folders
  t1 <- Sys.time()
  for(j in 1:100){
    dir.create(file.path("/Home/ii/rafaelc/Desktop/copy/R/garbage/testdir"), showWarnings = FALSE)
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Create_Folders",td))
  
  
  #Unzip the file
  t1 <- Sys.time()
  gunzip("1_S1_L001_R1_001.fastq.gz", overwrite = TRUE, remove = FALSE)
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Unzipping_File",td))
  
  #Modulus operation
  t1 <- Sys.time()
  for(j in 1:100){
    temp <- sample(1:10, 1)%%4
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Modulus",td))
  
  
  #Playing with Strings
  t1 <- Sys.time()
  for(j in 1:100){
    unlist(strsplit("Random String to Unlist","",fixed=TRUE))
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Unlist_Split",td))
  
  
  #Binding into matrix
  t1 <- Sys.time()
  for(j in 1:100){
    rbind(matrixToy, randomDataframe[j,j])
  }
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Binding",td))
  
  
  
  }
  t4 <- Sys.time()
  td2 <- as.numeric(t4-t3, units = "secs")
  matrixTiming <- rbind(matrixTiming, c("Complete_line",td2))
  
  
  
  

}



#Get the timing into a dataframe and plot stuff
timing <- data.frame(matrixTiming)
colnames(timing) <- c("Concept","Time_s")
timing[, 2] <- as.numeric(as.character( timing[, 2] ))

allConcepts <- levels(factor(timing$Concept))
matrixSubsets<- matrix(ncol=1, nrow=length(allConcepts))
for (i in 1:nrow(matrixSubsets)){
  print(allConcepts[i])  
  print(quantile(subset(timing, Concept==allConcepts[i])$Time_s))
  print(nrow(subset(timing, Concept==allConcepts[i])))
}





#Find all the lines that said "Config_Line_Complete"
configLinesReferences <- subset(timing, Concept=="Complete_line")

#We need to se what happends between all of those.
matrixPercentage <- NULL
for (i in 2:nrow(configLinesReferences)){
  
  #Get the total time
  totalTime = configLinesReferences[i,"Time_s"]
  
  #Get all the concept between those two
  boundaryRowA = as.integer(row.names(configLinesReferences[i,])) -1
  boundaryRowB = as.integer(row.names(configLinesReferences[i-1,])) +1
  
  interval <- timing[boundaryRowB:boundaryRowA,]
  
  #Sum all the timing of each concept and Find the % of the total time
  #Fill the data into a new matrix that will have columns [Concept, % time between config lines]
  subSetConcepts <- levels(factor(interval$Concept))
  
  conceptSumVector <- matrix(ncol=length(subSetConcepts), nrow=1)
  
  for(j in 1:length(subSetConcepts)){
    
    conceptSumVector[j] <- sum(subset(interval, Concept==subSetConcepts[j])$Time_s)
    conceptSumVector[j] <- conceptSumVector[j]/totalTime
    matrixPercentage <- rbind(matrixPercentage, c(subSetConcepts[j],conceptSumVector[j]))
    
  }
  
}
#When everything is finish show statistics of those concepts and %
percentages <- data.frame(matrixPercentage)
colnames(percentages) <- c("Concept","Time_s")
percentages[, 2] <- as.numeric(as.character( percentages[, 2] ))
allConcepts <- levels(factor(percentages$Concept))
matrixSubsets<- matrix(ncol=1, nrow=length(allConcepts))
for (i in 1:length(allConcepts)){
  print(allConcepts[i])
  
  print(quantile(subset(percentages, Concept==allConcepts[i])$Time_s))
  print(nrow(subset(percentages, Concept==allConcepts[i])))
}

#Transform into a boxplot graphic
test.m <- melt(percentages)

myBoxplot <- ggplot(test.m, aes(factor(variable), value)) + geom_boxplot(aes(fill = Concept))
myBoxplot

