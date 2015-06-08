#Example of Alignments in R for dummies:

#Some code from http://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.R

decipher_installed <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.


parallel_installed <- TRUE  # Set this to TRUE if you already install everything
                            # otherwise FALSE and the thing should install on its own.

total_processors <- 4       # Set this to the number of processors you want to use.
                            # You can set it to NULL if you want all that are available.

sequence_size = 200          # This is the size of the sequences

total_mutations = 5         # This is the total of mutations you want to introduce.
                            # A mutation can be a deletion, so your second sequence will have size equal
                            # or less than the first one.

                            # Size SeqA >= Size SeqB >= Size SeqA - total mutations

thisIsATest = TRUE          # Make the examples way more smaller to test the code

############################################
#Install the packages if necessary
############################################

if(!decipher_installed){
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("DECIPHER")

}


if(!parallel_installed){
  
  install.packages("foreach")
  install.packages("doMC")
  
}



############################################
#Summons the libraries and get to use functions from there.
############################################
library(DECIPHER)
library(foreach)
library(doMC)
registerDoMC(total_processors)


############################################
# 1.- Align two random sequences
############################################

matrixTiming <- NULL
t1 <- Sys.time()
t2 <- Sys.time()
td <- difftime(t2,t1)

pBar <- txtProgressBar(style=3)
for(i in 1:100){
#Get a random sequences:
sequenceA <- DNAStringSet(paste(sample(DNA_ALPHABET[1:4], sequence_size, replace = TRUE), collapse = ""))

#Get a second sequence, based on the first one with some random mutations
#By chance, the second sequence could mutate and be equal to the first one again ((1 - 4/5)^total_mutations chance)
positions <- sample(20, total_mutations) #Take total mutations amount of numbers, from 0 to SIZE, no repetitions.
# cat(positions)

ranges <- IRanges(positions, width=1) #Get them into a table form.
# print(ranges)

DNASet <- c(DNA_ALPHABET[1:4], "") #Get the DNA characters and add ""(nothing) into an array of characters.
# print(DNASet)

newNucleotides <- sample(DNASet, total_mutations, replace = TRUE) #Get a totally new sequence of size total mutations.
# cat(newNucleotides)

sequenceB <- replaceAt(sequenceA,at=ranges,newNucleotides) #Replace 
# print(sequenceB)

# print("Your first DNA sample is this one:")
# print(sequenceA)
# print("Your second DNA sample is this one:")
# print(sequenceB)
  
setTxtProgressBar(pBar, i/100)

t1 <- Sys.time()
myAlignment <- pairwiseAlignment(sequenceA, sequenceB)
t2 <- Sys.time()
td <- as.numeric(t2-t1, units = "secs")
matrixTiming <- rbind(matrixTiming, c("Read_Config_File",td))

}

print("Your alignment is this:")
print(myAlignment)

#Get the timing into a dataframe and plot stuff
timing <- data.frame(matrixTiming)
colnames(timing) <- c("Concept","Time_s")
timing[, 2] <- as.numeric(as.character( timing[, 2] ))

allConcepts <- levels(factor(timing$Concept))
matrixSubsets<- matrix(ncol=1, nrow=length(allConcepts))
for (i in 1:nrow(matrixSubsets)){
  print(allConcepts[i])
  #matrixSubsets[i,1]  <- subset(timing, Concept==allConcepts[i])
  #quantile(matrixSubsets[i,1]$Time_s)
  
  print(quantile(subset(timing, Concept==allConcepts[i])$Time_s))
#   print(nrow(subset(timing, Concept==allConcepts[i])))
}

############################################
# 2.- Align many sequences
############################################

# #Get some FASTA file from the Internet 
# fastaFile <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
# fastaData <- readDNAStringSet(fastaFile)
# fastaData # This is still unaligned
# 
# #Get only 50 sequences for testing porpuses (so it doesn't take minutes to construct the later tree)
# test <- DNAStringSet(c(fastaData[1:50]))
# test
# 
# #Form a guide tree using pairwiseAlignment
# if(!thisIsATest){
#     totalSequences <- length(fastaData)
# } else{
#     totalSequences <- length(test)
# }
# 
#   
#   d <- matrix(0, nrow=totalSequences, ncol=totalSequences) #Make a totalSequences x totalSequences new matrix
#   pBar <- txtProgressBar(style=3) #Make a progress bar for the impatients and anxious of heart
#   
#   #Make a Pair Wise Aligment between every sequence  
#   
#   counter <- 0
# 
#   system.time({
#       print("Parallel for")
#       foreach(j=2:totalSequences) %dopar% {
#           d[j, 1:(j - 1)] <- pairwiseAlignment(rep(fastaData[j], j - 1),fastaData[1:(j - 1)],scoreOnly=TRUE)
#           counter <- counter + 1
#           counter
#           setTxtProgressBar(pBar, counter/totalSequences)
#       }
#   })
# 
# 
# 
#   system.time({
#       print("Sequential for")
#       for (j in 2:totalSequences) {
#           d[j, 1:(j - 1)] <- pairwiseAlignment(rep(fastaData[j], j - 1),fastaData[1:(j - 1)],scoreOnly=TRUE)
#           setTxtProgressBar(pBar, j/totalSequences)
#       }
#   })
#   
#  
#   close(pBar)
#   
#   # rescale the distance scores from 0 to 1
#   m <- max(d)
#   
#   d[lower.tri(d)] <- d[lower.tri(d)] - m
#   m <- min(d)
#   d[lower.tri(d)] <- d[lower.tri(d)]/m
#   
#    
#   # form a guide tree from the distance matrix
#   gT <- IdClusters(d, cutoff=seq(0, 1, 0.01))
# 
#   system.time({ print("Timing the aligments") #Dummy System time; why are they not working?
#   })
# 
#   if(thisIsATest){
#       
#     # use the guide tree as input for alignment
#     system.time({
#       DNA <- AlignSeqs(test, guideTree=gT) # align directly
#     })
#     
#     system.time({
#       DNA <- AlignSeqs(test, guideTree=gT, processors = total_processors) # align directly
#     })
#     
#     DNA <- AlignTranslation(test, guideTree=gT) # align by translation
#   
#   }else{
#     
#     # use the guide tree as input for alignment
#     system.time({
#       DNA <- AlignSeqs(dna, guideTree=gT) # align directly
#     })
#     
#     system.time({
#       DNA <- AlignSeqs(dna, guideTree=gT, processors = total_processors) # align directly
#     })
#     
#     DNA <- AlignTranslation(dna, guideTree=gT) # align by translation
#     
#   }
# 
# print(DNA)