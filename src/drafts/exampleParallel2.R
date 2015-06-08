#Example of Alignments in R for dummies:

#Some code from http://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.R

decipher_installed <- TRUE  # Set this to TRUE if you already install everything
# otherwise FALSE and the thing should install on its own.


parallel_installed <- TRUE  # Set this to TRUE if you already install everything
# otherwise FALSE and the thing should install on its own.

total_processors <- 2       # Set this to the number of processors you want to use.
# You can set it to NULL if you want all that are available.

sequence_size = 20          # This is the size of the sequences in the first example

total_mutations = 5         # This is the total of mutations you want to introduce.
                            # A mutation can be a deletion, so your second sequence will have size equal
                            # or less than the first one.

                            # Size SeqA >= Size SeqB >= Size SeqA - total mutations

thisIsATest = TRUE          # Make the examples way more smaller to test the code

testSize = 50               # Set the number of sequences you want to try in the test (The total of the FASTA
                            # file is 318 (slow!) it will make a diagonal matrix of testSize^2)

timeSamples = 10            # How many time samples you want to take.

doFirstExample = FALSE      # Run the first example (or not)

############################################
#Install the packages if necessary
############################################

if(!decipher_installed){
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("DECIPHER")
  
}

if(!parallel_installed){
  
  install.packages("doParallel")
  
}

############################################
#Summons the libraries and get to use functions from there.
############################################
library(DECIPHER)
library(doParallel)
cl <- makeCluster(total_processors, outfile="") #outfile tells the workers to print on the terminal
registerDoParallel(cl)


############################################
# 1.- Align two random sequences
############################################

if (doFirstExample){

  #Get a random sequences:
  sequenceA <- DNAStringSet(paste(sample(DNA_ALPHABET[1:4], 20, replace = TRUE), collapse = ""))
  
  #Get a second sequence, based on the first one with some random mutations
  #By chance, the second sequence could mutate and be equal to the first one again ((1 - 4/5)^total_mutations chance)
  positions <- sample(20, total_mutations) #Take total mutations amount of numbers, from 0 to SIZE, no repetitions.
  cat(positions)
  
  ranges <- IRanges(positions, width=1) #Get them into a table form.
  ranges
  
  DNASet <- c(DNA_ALPHABET[1:4], "") #Get the DNA characters and add ""(nothing) into an array of characters.
  DNASet
  
  newNucleotides <- sample(DNASet, total_mutations, replace = TRUE) #Get a totally new sequence of size total mutations.
  cat(newNucleotides)
  
  sequenceB <- replaceAt(sequenceA,at=ranges,newNucleotides) #Replace 
  sequenceB
  
  print("Your first DNA sample is this one:")
  print(sequenceA)
  print("Your second DNA sample is this one:")
  print(sequenceB)
  
  myAlignment <- pairwiseAlignment(string1, string2)
  
  print("Your alignment is this:")
  print(myAlignment)

}

############################################
# 2.- Align many sequences
############################################

#Get some arrays to store the timing data
parallel_for <- rep(0, timeSamples)
sequential_for <- rep(0, timeSamples)
parallel_aligment <- rep(0, timeSamples)
sequential_aligment <- rep(0, timeSamples)

#Get the sequence size for each test
sequence_size <- seq(5, testSize, ceiling(testSize/timeSamples))

#Get some FASTA file from the Internet 
fastaFile <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
fastaData <- readDNAStringSet(fastaFile)
fastaData # This is still unaligned


for (i in 1:length(sequence_size)){
  
  print("You are now in the iteration: ")
  print(i)
  
  #Get only 50 sequences for testing porpuses (so it doesn't take minutes to construct the later tree)
  test <- DNAStringSet(c(fastaData[1:sequence_size[i]]))
  test
  
  #Form a guide tree using pairwiseAlignment
  if(!thisIsATest){
    totalSequences <- length(fastaData)
  } else{
    totalSequences <- sequence_size[i]
  }
  
  
  d <- matrix(0, nrow=totalSequences, ncol=totalSequences) #Make a totalSequences x totalSequences new matrix
  
  pBar <- txtProgressBar(style=3) #Make a couples of progress bar for the impatients and anxious of heart
  pb <- txtProgressBar(min = 1, max = totalSequences, style=3)
  
  #Make a Pair Wise Aligment between every sequence  
  
  counter <- 0
  
  parallel_for[i] <- system.time({
    print("Parallel for")
    foreach(j=2:totalSequences , .packages="DECIPHER") %dopar% {
      d[j, 1:(j - 1)] <- pairwiseAlignment(rep(fastaData[j], j - 1),fastaData[1:(j - 1)],scoreOnly=TRUE)
      #counter <- counter + 1
      #counter
      #setTxtProgressBar(pBar, counter/totalSequences)
      
      setTxtProgressBar(pb, j)
    }
  })
  
  sequential_for[i] <- system.time({
    print("Sequential for")
    foreach (j=2:totalSequences) %do% {
      d[j, 1:(j - 1)] <- pairwiseAlignment(rep(fastaData[j], j - 1),fastaData[1:(j - 1)],scoreOnly=TRUE)
      setTxtProgressBar(pBar, j/totalSequences)
    }
  })
  
  close(pBar)
  
  # rescale the distance scores from 0 to 1
  m <- max(d)
  
  d[lower.tri(d)] <- d[lower.tri(d)] - m
  m <- min(d)
  d[lower.tri(d)] <- d[lower.tri(d)]/m
  
  
  # form a guide tree from the distance matrix
  gT <- IdClusters(d, cutoff=seq(0, 1, 0.01))
  
  
  if(thisIsATest){
    
    # use the guide tree as input for alignment
    
    sequential_aligment[i] <- system.time(AlignSeqs(test, guideTree=gT)) # align directly
    
    parallel_aligment[i] <- system.time(AlignSeqs(test, guideTree=gT, processors = total_processors)) # align directly
    
    DNA <- AlignTranslation(test, guideTree=gT) # align by translation
    
  }else{
    
    # use the guide tree as input for alignment
    DNA <- AlignSeqs(dna, guideTree=gT) # align directly
    
    DNA <- AlignSeqs(dna, guideTree=gT, processors = total_processors) # align directly
    
    DNA <- AlignTranslation(dna, guideTree=gT) # align by translation
    
  }
  
  print(DNA)

}

print("Parallel for timing")
parallel_for
print("Sequential for timing")
sequential_for
print("Parallel aligment timing")
parallel_aligment
print("Sequential aligment timing")
sequential_aligment

print("Sequence sizes")
sequence_size