library(amplican)
context("Warning functions")

# Invariant: checkTarget
#
#   The file descriptor must be set to append.s

# Invariant: checkPrimers
#
#   The file descriptor must be set to append.

# Invariant: checkPositions
#
#   The file descriptor must be set to append.
#   The alignmentPositions array can also be a single value. For R, every single variable is an array of size one.
#   The alignmentPositions must be 0 or less to indicate that is not found. The first character in the amplicon is
#   the character one, not zero.

# Invariant: checkConfigFile
#
#    The dataframe must have the following columns:
#      "ID","Barcode","Forward_Reads_File","Reverse_Reads_File","Experiment_Type", "Target_Primer",
#      "Forward_Primer","Reverse_Primer","Strand","Genome".
#    The config table is not empty.
