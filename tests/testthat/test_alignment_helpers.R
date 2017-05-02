library(amplican)
library(testthat)
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

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("AGGGTAAAGTCCATGGCCCCAATTTGTGTGTAG"),
                                          Biostrings::DNAString("AGTGAAGTCAAACATGGAATTAGTGTGTTAA"), type = "global")
  events_ranges <- GenomicRanges::GRanges(seqnames = "test",
                                          ranges = IRanges::IRanges(c(5, 18, 10, 29, 3, 22, 31),
                                                                    c(6, 21, 12, 29, 3, 22, 31)),
                                          strand = rep("+", 7),
                                          originally = c("", "", "", "", "T", "A", "A"),
                                          replacement = c("TA", "CCCC", "", "", "G", "T", "G"),
                                          type = c("insertion", "insertion",
                                                   "deletion", "deletion",
                                                   rep("mismatch", 3)))
  expect_identical(getEventInfo(events, "test", 1), events_ranges )
  strand(events_ranges) <- "-"
  expect_identical(getEventInfo(events, "test", 31, strand_info = "-"), events_ranges)

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("AGGGTAAAAGTCCATGGCCCAATTTGTGTGTAG"),
                                          Biostrings::DNAString("CCCCCCCCCCCAGTGAAGTCAAACATGGAATTAGTGTGTTAA"), type = "global")
  events_ranges <- GenomicRanges::GRanges(seqnames = "test",
                                          ranges = IRanges::IRanges(c(16, 29, 21, 40, 14, 33, 42),
                                                                    c(18, 31, 23, 40, 14, 33, 42)),
                                          strand = rep("+", 7),
                                          originally = c("", "", "", "", "T", "A", "A"),
                                          replacement = c("TAA", "CCC", "", "", "G", "T", "G"),
                                          type = c("insertion", "insertion",
                                                   "deletion", "deletion",
                                                   rep("mismatch", 3)))
  expect_identical(getEventInfo(events, "test", 12), events_ranges)

  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTGGGGGGGGGGACTGGGGGGGGGGACT"),
                                          Biostrings::DNAString("ACTACTACT"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                           ranges = IRanges::IRanges(c(4, 7), c(13, 16)),
                                                                           strand = rep("+", 2),
                                                                           originally = c("", ""),
                                                                           replacement = c("GGGGGGGGGG", "GGGGGGGGGG"),
                                                                           type = c("insertion", "insertion")))
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTACTTCT"),
                                          Biostrings::DNAString("ACTGGGGGGGGGGACTACT"), type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges(seqnames = "test",
                                                                           ranges = IRanges::IRanges(c(4, 17), c(13, 17)),
                                                                           strand = rep("+", 2),
                                                                           originally = c("", "A"),
                                                                           replacement = c("", "T"),
                                                                           type = c("deletion", "mismatch")))
})
