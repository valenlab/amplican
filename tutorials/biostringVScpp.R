
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

setwd("~/Desktop/git/ampliCan/src/tutorials") #Workaround for testing in cwd, because RStudio hates me
sourceCpp("nw.cpp")

#This array represent the sequence size for each of the timing comparison.
SEQUENCE_ARRAY <- c(10,20,50,100,200,500,1000,2000,5000)

#This variable represent the sample for each of size
#In total, we will make samples_size * size(sequence_array) alignments
SAMPLES_SIZE = 100

#This variable tell us the mutation rate of the two sequences. It shouldn't be
#in between [0 and 100). If you set it to 100 the program has an infenitesimal
#chance of create an empty sequence and crash.
MUTATION_RATIO = 10

# Here are the matrices where we are going to store the timing results
matrixResults <- NULL
matrixTiming <- NULL

# File where to write the raw results
sourceFD <-  file(paste(htmlPath,"/htmlResults_",i,".html",collapse = '',sep = '') , open="w")

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
for(j in 1:length(SEQUENCE_ARRAY)){

  # Get the size for this test and the amount of mutations
  current_sequence_size = SEQUENCE_ARRAY[j]
  total_mutations <- ceiling(current_sequence_size/MUTATION_RATIO)
  
  # Initialize the result string where the C++ return the result and
  # set the matrix for timing back to 0.
  results <- ""
  matrixTiming <- NULL
  
  #Repeat the following experiment many times:
  for(i in 1:SAMPLES_SIZE){

    # Advance the progress bar
    setTxtProgressBar(pBar, i/SAMPLES_SIZE)
    
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
    matrixTiming <- rbind(matrixTiming, c("Align_Biostring",td))
        
    writeLines(htmlSource, sourceFD)
    
    # Make some variables for the C++ function.
    pattern <- toString(sequenceA)
    subject <-  toString(sequenceB)
        
    # Make an alignment using the C++ algorithm
    # And put the result in the timing matrix with the label Aling_CPP
    t1 <- Sys.time()
    result <- nwRCPP(pattern, subject, scoring_matrix, gapOpening, gapExtension, gapEnding)
    t2 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    matrixTiming <- rbind(matrixTiming, c("Align_CPP",td))
    
  }
  
  # Now that the experiment is over, we have a matrix in the form of
  # Concept , Time
  # Biostring 7.6
  # Biostring 4.5
  # ...
  # CPP       1.3
  # CPP       2.3
  
  # Put that into a dataframe
  timing <- data.frame(matrixTiming)
  colnames(timing) <- c("Concept","Time_s")
  timing[, 2] <- as.numeric(as.character( timing[, 2] ))
  
  # Get the result for each of the groups.
  resultsBio <- summary(subset(timing, Concept=="Align_Biostring")$Time_s)
  resultsCPP <- summary(subset(timing, Concept=="Align_CPP")$Time_s)
  
  # Finally, put the results into the matrix of the results.
  # In here, we annotate the Quartiles and the minimum and maximum for that test.
  
  matrixResults <- rbind(matrixResults, c("Align_CPP","MIN",current_sequence_size,resultsCPP[1]))
  matrixResults <- rbind(matrixResults, c("Align_CPP","Q1",current_sequence_size,resultsCPP[2]))
  matrixResults <- rbind(matrixResults, c("Align_CPP","Q2",current_sequence_size,resultsCPP[3]))
  matrixResults <- rbind(matrixResults, c("Align_CPP","Q3",current_sequence_size,resultsCPP[5]))
  matrixResults <- rbind(matrixResults, c("Align_CPP","MAX",current_sequence_size,resultsCPP[6]))
  
  matrixResults <- rbind(matrixResults, c("Align_Bio","MIN",current_sequence_size,resultsBio[1]))
  matrixResults <- rbind(matrixResults, c("Align_Bio","Q1",current_sequence_size,resultsBio[2]))
  matrixResults <- rbind(matrixResults, c("Align_Bio","Q2",current_sequence_size,resultsBio[3]))
  matrixResults <- rbind(matrixResults, c("Align_Bio","Q3",current_sequence_size,resultsBio[5]))
  matrixResults <- rbind(matrixResults, c("Align_Bio","MAX",current_sequence_size,resultsBio[6]))

  # Give some feedback to the user of where are we at the moment
  print("Finish for:")
  print(SEQUENCE_ARRAY[j])
  print("Biostring:")
  print(resultsBio)
  print("C++")
  print(resultsCPP)
  
}

# Transform the matrix with the result into a nice looking graphic with the comparison

#First, we make a dataframe
results <- data.frame(matrixResults)
colnames(results) <- c("Algorithm","Concept","Sequence_size","Timing")
results[, 3] <- as.numeric(as.character(results[, 3]))
results[, 4] <- as.numeric(as.character(results[, 4]))