library(amplican)
context("directory helper functions")
# Invariant: unifyFiles
#
#   The returned string must have the same length as the input string.

#
# Invariant: getReadsFile
#
#   The unzipped file contain a plain text file. This file is formatted with a
#   FASTQ format, meaning that is repeating the following pattern over an over
#   again:
#   - Line with the ID and other metadata
#   - Line with the sequence [01]
#   - '+ line'
#   - Line with the quality of the nucleotides [02]
#
#   This means that the number of lines in that file must be divisible by 4.
#
#   Worth mentioning that this function is suppose to be use for the forward
#   and reverse reads. The number of lines in each of those FASTQ files must
#   be the same.

#
# Invariant: getConfigFolder
#
#   The file path must be Unix style

# Invariant: getConfigName
#
#   The file path must be Unix style

# Invariant:  checkReady
#  The file path must be Unix style
