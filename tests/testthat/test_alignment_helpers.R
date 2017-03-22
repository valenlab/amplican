library(amplican)
context("alignment helper functions")

test_that("upperGroups returns correct IRanges", {
  expect_identical(upperGroups("aaaccTTTTGGggg"), IRanges::IRanges(6, 11))
  expect_identical(upperGroups("aaaccTTTTGGgggAAAccc"), IRanges::IRanges(c(6, 15), c(11, 17)))
  expect_identical(upperGroups("aaaccttttggggg"), IRanges::IRanges())
})

test_that("getEventInfo returns correct GRanges", {
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTG"),
                                          Biostrings::DNAString("ACTG"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges())

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAGT"),
                                          Biostrings::DNAString("ACTGAGT"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                        ranges = IRanges::IRanges(4, 4),
                                                                        strand = "+",
                                                                        originally = "",
                                                                        replacement = "",
                                                                        type = "deletion"))

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAAGT"),
                                          Biostrings::DNAString("ACTGAGT"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                        ranges = IRanges::IRanges(4, 4),
                                                                        strand = "+",
                                                                        originally = "G",
                                                                        replacement = "A",
                                                                        type = "mismatch"))

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAAAGT"),
                                          Biostrings::DNAString("ACTGAGT"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                        ranges = IRanges::IRanges(c(6, 4), c(6, 4)),
                                                                        strand = c("+", "+"),
                                                                        originally = c("", "G"),
                                                                        replacement = c("A", "A"),
                                                                        type = c("insertion", "mismatch")))

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("AGGTAAAGTAATTTGTGTGTAA"),
                                          Biostrings::DNAString("ACTGAGTTTTTTAGTGTGTTAA"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                        ranges = IRanges::IRanges(c(21, 3, 2, 5, 7, 8, 10, 11, 14),
                                                                                                  c(21, 3, 2, 5, 7, 8, 10, 11, 14)),
                                                                        strand = rep("+", 9),
                                                                        originally = c("", "", "C", "G", "G", "T", "T", "T", "A"),
                                                                        replacement = c("", "G", "G", "A", "A", "G", "A", "A", "T"),
                                                                        type = c("deletion",
                                                                                 "insertion",
                                                                                 rep("mismatch", 7))))
  expect_identical(getEventInfo(events, "test", 22, strand_info = "-"), GenomicRanges::GRanges(seqnames = "test",
                                                                                           ranges = IRanges::IRanges(c(21, 3, 2, 5, 7, 8, 10, 11, 14),
                                                                                                                     c(21, 3, 2, 5, 7, 8, 10, 11, 14)),
                                                                                           strand = rep("-", 9),
                                                                                           originally = c("", "", "C", "G", "G", "T", "T", "T", "A"),
                                                                                           replacement = c("", "G", "G", "A", "A", "G", "A", "A", "T"),
                                                                                           type = c("deletion",
                                                                                                    "insertion",
                                                                                                    rep("mismatch", 7))))
})
