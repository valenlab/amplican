############################################
#Summons the libraries and get to use functions from there.
############################################
library(DECIPHER)
library(doParallel)
library(R.utils) #Unzipping the .gz files
library(Biostrings)
cl <- makeCluster(total_processors, outfile="") #outfile tells the workers to print on the terminal
registerDoParallel(cl)

#Check if there is a folder call /result in cwd, if not create a new one
# dir.create(file.path(getwd(), "/results"), showWarnings = FALSE)

#Make a folder call /run_POSIX_TIME where we are going to write the results for this run.
# currentTime <- Sys.time()
# timeStamp <-  strftime(currentTime,"%Y%m%d%H%M%S")
# resultsFolder <- paste(getwd(), "/results/", timeStamp, sep = '')
# dir.create(file.path(resultsFolder), showWarnings = TRUE)

#The folder structure is like this:
#data:
#Data_folder_1 (ie: Cas9_data)
# ...
#Data_folder_X
#results:
#run_??:
#Data_folder_1_results (ie: Cas9_data_results)
#----Config_file_1_results (ie: run11_results)
#-------- Results 1 (ie: 1_S1_L001_R1_001,1_S1_L001_R2_001)
#------------ Alignments
#---------------- Alignment_1.txt
#---------------- Alignment_Z.txt
#------------ StatsForwardsTable.txt [Instrument ID,..., X, Y, # IN, # DEL, # MISS, Alignment txt file]
#------------ StatsReverseTable.txt  [Instrument ID,..., X, Y, # IN, # DEL, # MISS, Alignment txt file]
#-------- Results Y
#-------- Config_file_1_results.txt (ie: run11_results.txt) [ID,Barcode,FR_file,RR_file,Insertions,Deletions,Mismatches]


# exists("variable")
# 
# print("Here 1")
# 
# variable <- 10
# 
# print("Here 2")
# 
# exists("variable")
# 
# print("Here 3")
# 
# print(variable)
# 
# grep("A", c("b","A","c"), fixed=TRUE)
# 
# grep("A", c("bAaa","A","caaa"), fixed=TRUE)
# 
# list1 <- grep("^data$", c("data","bee","gyy"), fixed=FALSE, ignore.case=FALSE)
# length(list1)
# 
# 
# list2 <- grep("^data$", c("sata","bee","gyy","mydata"), fixed=FALSE, ignore.case=FALSE)
# length(list2)


# 
# 
# print((0-ceiling(0/4)*4))
# print((1-ceiling(1/4)*4))
# print((2-ceiling(2/4)*4))
# print((3-ceiling(3/4)*4))
# print((4-ceiling(4/4)*4))
# print((5-ceiling(5/4)*4))
# print((6-ceiling(6/4)*4))
# 
# 0 %% 4
# 1 %% 4
# 2 %% 4
# 3 %% 4
# 4 %% 4
# 5 %% 4
# 6 %% 4
# 7 %% 4
# 8 %% 4
# 
# a <- unlist(strsplit("abc", "", fixed = TRUE))
# 
# print(a)
# print(a[1])
# print(a[2])
# print(a[3])
# 
# 
# regexpr("asdfas","ttttttttttttttttt",fixed=TRUE)

# print("-------------")
# 
# a <- regexpr("GCCTGTGTCGTCTCCCTTATG" , "AGTTCATGTGTTTGCGTTTTGCTATCGCTCCTGTTTCTGTCCGTATCAGCAGATCAGTGTGACAATGACAACACAGCTGACAACTGCAATAAATGGACCGAAACATCTCACAATCAGCCTCGAAATACAGGTAAAGTAGACTTTAACACC" , fixed=TRUE)
# b <- regexpr("GTGTTGGTGTGGAGGAGGA" , "AGTTCATGTGTTTGCGTTTTGCTATCGCTCCTGTTTCTGTCCGTATCAGCAGATCAGTGTGACAATGACAACACAGCTGACAACTGCAATAAATGGACCGAAACATCTCACAATCAGCCTCGAAATACAGGTAAAGTAGACTTTAACACC" , fixed=TRUE)
# c <- regexpr("GGGGCGGCACTTCAGGGCAGTGG" , "AGTTCATGTGTTTGCGTTTTGCTATCGCTCCTGTTTCTGTCCGTATCAGCAGATCAGTGTGACAATGACAACACAGCTGACAACTGCAATAAATGGACCGAAACATCTCACAATCAGCCTCGAAATACAGGTAAAGTAGACTTTAACACC" , fixed=TRUE)
# 
# print(a[1])
# print(b[1])
# print(c[1])
# 
# a <- regexpr("GCCTGTGTCGTCTCCCTTATG" , "GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTATGAATGTGTTTTTTGTGAATGTAATAATAT" , fixed=TRUE)
# b <- regexpr("TCCTCCTCCACACCAACAC" , "TCCTCCTCCACACCAACACAACTACACAAAAAGAAAAAACATATATTATTACATTCACAAAAAACACATTCATAAGTACAAGATTTAGCTTTTAACAGTCATTTAAAAATTTTGCAACCTCACCTGCCATGANNNNNNNNNNNNNNNNNN" , fixed=TRUE)
# 
# print(a[1])
# print(b[1])
# 
# a[1]+2

#dna <- DNAStringSet("GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTATGAATGTGTTTTTTGTGAATGTAATAATAT")
# 
# dna <- DNAStringSet(c("GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTATGAATGTGTTTTTTGTGAATGTAATAATAT","TCCTCCTCCACACCAACACAACTACACAAAAAGAAAAAACATATATTATTACATTCACAAAAAACACATTCATAAGTACAAGATTTAGCTTTTAACAGTCATTTAAAAATTTTGCAACCTCACCTGCCATGANNNNNNNNNNNNNNNNNN","gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga"))
# print(dna)
# 
# dna[[1]]
# 
# dna[[1]]
# 
# as.character(dna)[1]
# 
# 
# 
genome <- "gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga"
candidateForwardSequence <- "GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTATGAATGTGTTTTTTGTGAATGTAATAATAT"
candidateComplementarySequence <- reverseComplement(DNAStringSet( "TCCTCCTCCACACCAACACAACTACACAAAAAGAAAAAACATATATTATTACATTCACAAAAAACACATTCATAAGTACAAGATTTAGCTTTTAACAGTCATTTAAAAATTTTGCAACCTCACCTGCCATGANNNNNNNNNNNNNNNNNN"))

# print(as.character(candidateComplementarySequence)[1])
# print(as.character(genome)[1])
# print(as.character(candidateForwardSequence)[1])
# 
# allCandidates <- DNAStringSet(c(genome,candidateForwardSequence,as.character(candidateComplementarySequence)[1]))
# 
# print(allCandidates)

# myPairwise <- pairwiseAlignment(candidateForwardSequence, genome)
# print(myPairwise)
# 
# show(myPairwise)
# 
# attributes(myPairwise)[-5] # Returns the components of an alignment object.
# 
# print("Matrix")
# as.matrix(myPairwise)
# print("Insertion")
# insertion(myPairwise)
# print("Deletion")
# deletion(myPairwise)
# print("Indel")
# indel(myPairwise)
# print("Nindel")
# nindel(myPairwise)
# print("Score")
# score(myPairwise)
# 
# writePairwiseAlignments(myPairwise, file="test.txt", block.width=500) # Writes sequences to file in FASTA format.

# 
# 
# myAlignment <- AlignSeqs(allCandidates, verbose = FALSE)
# print(myAlignment)
# 
# 

# 
# attributes(myAlignment)[-5] # Returns the components of an alignment object.
# 
# print("Matrix")
# as.matrix(myAlignment)
# # print("Insertion")
# # insertion(myAlignment)
# # print("Deletion")
# # deletion(myAlignment)
# # print("Indel")
# # indel(myAlignment)
# # print("Nindel")
# # nindel(myAlignment)
# # print("Score")
# # score(myAlignment)
# 


# 
# dna1 <- DNAString("ACCGC")
# dna2 <- DNAString("TAC")
# print (dna1)
# print (dna2)
# 
# dna3 <- DNAStringSet(c("ACCGC", "TAC"))
# print(dna3)

print(candidateComplementarySequence)
print(candidateForwardSequence)

dna4 <- DNAStringSet(c(toString(candidateComplementarySequence), candidateForwardSequence))
print(dna4)


myPairwise2 <- pairwiseAlignment(dna4, genome)
print(myPairwise2)

writePairwiseAlignments(myPairwise2, file="test2.txt", block.width=80) # Writes sequences to file in FASTA format.

print("Matrix")
as.matrix(myPairwise2)
# print("Insertion")
# insertion(myPairwise2)
# length(insertion(myPairwise2))
# length(insertion(myPairwise2)[[1]])
# length(insertion(myPairwise2)[[2]])
# print("Deletion")
# deletion(myPairwise2)
# length(deletion(myPairwise2))
# length(deletion(myPairwise2)[[1]])
# length(deletion(myPairwise2)[[2]])
print("Mistmatch")
mismatchTable(myPairwise2)
print("Mistmatch2")
nmatch(myPairwise2)
length(nmatch(myPairwise2))
print("Mistmatch3")
nmismatch(myPairwise2)[1]
nmismatch(myPairwise2)[2]
length(nmatch(myPairwise2))

length(grep("N", strsplit(toString(candidateComplementarySequence),"")[[1]], fixed=TRUE))
# pmatch("N", toString(candidateComplementarySequence), nomatch = 0)


