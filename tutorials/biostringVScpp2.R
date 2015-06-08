
############################################
# This script compare the running time of the sequence aligning algorithm that are in:
# Biostring in R
# Our developed Gotoh in C++
#
# The algorithm generate a random sequence of N size; called sequence A
# Then it modifies the sequence randomly. The new sequence is of size <= N
# and have a total of MUTATION_RATIO % changes from the original.

# For our developed version, we are using the function that allocate memory for two
# sequences only; not the one that allocate for a lot of sequences only one time.
# Since allocation is time expensive, the other version will be faster than the result
# of this test.

# The results is a graphic where the X-axis represent the sequence size and the Y-axis
# represent the time taken to accomplish that alignment.

############################################

library(DECIPHER)
library(Biostrings)
library(Rcpp)

# Lets start the seed
set.seed(as.numeric(as.POSIXct(Sys.time())))

# Set the boundaries
UPPER <- 10000
LOWER <- 9

setwd("~/Desktop/git/ampliCan/src/tutorials") #Workaround for testing in cwd, because RStudio hates me
sourceCpp("nw.cpp")

#This variable represent the sample for each of size
#In total, we will make samples_size * size(sequence_array) alignments
SAMPLES_SIZE = 10000

#This variable tell us the mutation rate of the two sequences. It shouldn't be
#in between [0 and 100). If you set it to 100 the program has an infenitesimal
#chance of create an empty sequence and crash.
MUTATION_RATIO = 10

# Here are the matrices where we are going to store the timing results
matrixResults <- NULL
matrixTiming <- NULL

# File where to write the raw results
rawFD <-  file(paste(getwd(),"/RAWData.txt",collapse = '',sep = '') , open="w")

#Some timing variables
t1 <- Sys.time()
t2 <- Sys.time()
td <- difftime(t2,t1)

# Make some variables for the C++ function.
scoring_matrix = "NUC44"
gapOpening = 50
gapExtension = 0
gapEnding = TRUE

#A progress bar for the impatients
pBar <- txtProgressBar(style=3)

#For each of the sequence sizes:
for(j in 1:SAMPLES_SIZE){

  sequence_size <- ceiling(runif(1, LOWER, UPPER))
  
  # Get the size for this test and the amount of mutations
  current_sequence_size <- sequence_size
  total_mutations <- ceiling(current_sequence_size/MUTATION_RATIO)
  
  # Initialize the result string where the C++ return the result and
  # set the matrix for timing back to 0.
  results <- ""
  
  # Advance the progress bar
  setTxtProgressBar(pBar, j/SAMPLES_SIZE)
    
  # Get a random sequences:
  sequenceA <- DNAStringSet(paste(sample(DNA_ALPHABET[1:4], current_sequence_size, replace = TRUE), collapse = ""))
    
  # Get a second sequence, based on the first one with some random mutations
  # By chance, the second sequence could mutate and be equal to the first one
  # again ((1 - 4/5)^total_mutations chance)
  positions <- sample(current_sequence_size, total_mutations) #Take total_mutations amount of numbers, from 0 to SIZE, no repetitions.    
  ranges <- IRanges(positions, width=1) #Get them into a table form.    
  DNASet <- c(DNA_ALPHABET[1:4], "") #Get the DNA characters and add ""(nothing) into an array of nucleotides characters.
  newNucleotides <- sample(DNASet, total_mutations, replace = TRUE) #Get a totally new sequence of size total mutations.
    
  sequenceB <- replaceAt(sequenceA,at=ranges,newNucleotides) #Replace in A the new nucleotides specify in ranges


  # Make an alignment using the Biostring/Dechiper algorithm
  # And put the result in the timing matrix with the label Aling_Biostring
  t1 <- Sys.time()
  myPairwise <- pairwiseAlignment(sequenceA, sequenceB, type = "global")
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")
        
  writeLines(paste("Biostring",td,sequence_size, sep='\t'), rawFD)
    
  # Make some variables for the C++ function.
  pattern <- toString(sequenceA)
  subject <-  toString(sequenceB)
        
  # Make an alignment using the C++ algorithm
  # And put the result in the timing matrix with the label Aling_CPP
  t1 <- Sys.time()
  result <- nwRCPP(pattern, subject, scoring_matrix, gapOpening, gapExtension, gapEnding)
  t2 <- Sys.time()
  td <- as.numeric(t2-t1, units = "secs")

  writeLines(paste("CPP",td,sequence_size, sep='\t'), rawFD)
  
}

close(rawFD)