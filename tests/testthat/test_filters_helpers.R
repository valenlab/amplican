library(amplican)
context("Filter functions")

# Invariant for (filterByNucleotides)
#
#   The dataframes must have the same length, so the array with the valid rows.
#
#   The dataframes must have a column named Sequence.
#
#   The returned array will have more FALSEs than validRows.
#
#   If validRows have a FALSE in position x, the returned array must have a
#   FALSE in x.

# test_that("str_length is number of characters", {
#   expect_equal(str_length("a"), 1)
#   expect_equal(str_length("ab"), 2)
#   expect_equal(str_length("abc"), 3)
# })

# Invariant for filterByQuality
#
#   The dataframes must have the same length, so the array with the valid rows.
#
#   The dataframes must have a column named Quality.
#
#   The quality is defined by the FASTQ format, being '!' = 33 the lowest quality and '~' = 126 the highest. However, the
#   arguments for functionare given from 0 to 93. If you set them to 94 or bigger, nothing will be valid. And if you set
#   them to 0 or less, everything will be valid.
#
#   The returned array will have more FALSEs than validRows.
#
#   If validRows have a FALSE in position x, the returned array must have a
#   FALSE in x.
