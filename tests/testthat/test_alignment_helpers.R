library(amplican)
library(testthat)
context("alignment helper functions")

test_that("comb_along returns correct number of elements", {
  expect_equal(length(comb_along("AC")), 4^2) # 4 letters of alphabet
  expect_equal(length(comb_along("AAA")), 37)
  expect_equal(length(comb_along("AAA", 3)), 4^3)
  expect_equal(length(comb_along("AAA", 1)), 10)
})

test_that("upperGroups returns correct IRanges", {
  expect_identical(upperGroups("aaaccTTTTGGggg"),
                   IRanges::IRanges(6, 11))
  expect_identical(upperGroups("aaaccTTTTGGgggAAAccc"),
                   IRanges::IRanges(c(6, 15), c(11, 17)))
  expect_identical(upperGroups("aaaccttttggggg"), IRanges::IRanges())
})

test_that("getEventInfo returns correct GRanges", {
  # no events
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTG"),
                                          Biostrings::DNAString("ACTG"),
                                          type = "global")
  expect_identical(getEventInfo(events, "test", 1), GenomicRanges::GRanges())

  # simple deletion
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAGT"),
                                          Biostrings::DNAString("ACTGAGT"),
                                          type = "global")
  gr1 <- GenomicRanges::GRanges(seqnames = "test",
                               ranges = IRanges::IRanges(4, 4),
                               strand = "+",
                               originally = "",
                               replacement = "",
                               type = "deletion",
                               read_id = "1")
  names(gr1) <- "1"
  test_gr1 <- getEventInfo(events, "test", 1)
  test_gr1$score <- NULL
  expect_identical(test_gr1, gr1)

  # simple mismatch
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAAGT"),
                                          Biostrings::DNAString("ACTGAGT"),
                                          type = "global")
  gr2 <- GenomicRanges::GRanges(seqnames = "test",
                               ranges = IRanges::IRanges(4, 4),
                               strand = "+",
                               originally = "G",
                               replacement = "A",
                               type = "mismatch",
                               read_id = "1")
  names(gr2) <- "1"
  test_gr2 <- getEventInfo(events, "test", 1)
  test_gr2$score <- NULL
  expect_identical(test_gr2, gr2)

  #simple mismatch + insertion
  events <- Biostrings::pairwiseAlignment(Biostrings::DNAString("ACTAAAGT"),
                                          Biostrings::DNAString("ACTGAGT"),
                                          type = "global")
  gr3 <- GenomicRanges::GRanges(seqnames = "test",
                               ranges = IRanges::IRanges(c(6, 4), c(6, 4)),
                               strand = c("+", "+"),
                               originally = c("", "G"),
                               replacement = c("A", "A"),
                               type = c("insertion", "mismatch"),
                               read_id = "1")
  names(gr3) <- c("1", "1")
  test_gr3 <- getEventInfo(events, "test", 1)
  test_gr3$score <- NULL
  expect_identical(test_gr3, gr3)

  # ins + del + ins + del
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString("AGGGTAAAGTCCATGGCCCCAATTTGTGTGTAG"),
    Biostrings::DNAString("AGTGAAGTCAAACATGGAATTAGTGTGTTAA"), type = "global")
  gr4 <- GenomicRanges::GRanges(
    seqnames = "test",
    ranges = IRanges::IRanges(c(5, 18, 10, 29, 3, 22, 31),
                              c(6, 21, 12, 29, 3, 22, 31)),
    strand = rep("+", 7),
    originally = c("", "", "", "", "T", "A", "A"),
    replacement = c("TA", "CCCC", "", "", "G", "T", "G"),
    type = c("insertion", "insertion",
             "deletion", "deletion",
             rep("mismatch", 3)),
    read_id = "1")
  names(gr4) <- rep("1", 7)
  test_gr4 <- getEventInfo(events, "test", 1)
  test_gr4$score <- NULL
  expect_identical(test_gr4, gr4)

  # minus strand
  GenomicRanges::strand(gr4) <- "-"
  test_gr4 <- getEventInfo(events, "test", 1, strand_info = "-")
  test_gr4$score <- NULL
  expect_identical(test_gr4, gr4)

  # overhang + ins + del + ins + del
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString("AGGGTAAAAGTCCATGGCCCAATTTGTGTGTAG"),
    Biostrings::DNAString("CCCCCCCCCCCAGTGAAGTCAAACATGGAATTAGTGTGTTAA"),
    type = "global")
  gr5 <- GenomicRanges::GRanges(
    seqnames = "test",
    ranges = IRanges::IRanges(c(16, 29, 21, 40, 14, 33, 42),
                              c(18, 31, 23, 40, 14, 33, 42)),
    strand = rep("+", 7),
    originally = c("", "", "", "", "T", "A", "A"),
    replacement = c("TAA", "CCC", "", "", "G", "T", "G"),
    type = c("insertion", "insertion",
             "deletion", "deletion",
             rep("mismatch", 3)),
    read_id = "1")
  names(gr5) <- rep("1", 7)
  test_gr5 <- getEventInfo(events, "test", 1)
  test_gr5$score <- NULL
  expect_identical(test_gr5, gr5)

  # ins + ins
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString("ACTGGGGGGGGGGACTGGGGGGGGGGACT"),
    Biostrings::DNAString("ACTACTACT"), type = "global")
  gr6 <- GenomicRanges::GRanges(seqnames = "test",
                               ranges = IRanges::IRanges(c(4, 7), c(13, 16)),
                               strand = rep("+", 2),
                               originally = c("", ""),
                               replacement = c("GGGGGGGGGG", "GGGGGGGGGG"),
                               type = c("insertion", "insertion"),
                               read_id = "1")
  names(gr6) <- c("1", "1")
  test_gr6 <- getEventInfo(events, "test", 1)
  test_gr6$score <- NULL
  expect_identical(test_gr6, gr6)

  # insertion + mismatch
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString("ACTACTTCT"),
    Biostrings::DNAString("ACTGGGGGGGGGGACTACT"), type = "global")
  gr7 <- GenomicRanges::GRanges(seqnames = "test",
                               ranges = IRanges::IRanges(c(4, 17), c(13, 17)),
                               strand = rep("+", 2),
                               originally = c("", "A"),
                               replacement = c("", "T"),
                               type = c("deletion", "mismatch"),
                               read_id = "1")
  names(gr7) <- c("1", "1")
  test_gr7 <- getEventInfo(events, "test", 1)
  test_gr7$score <- NULL
  expect_identical(test_gr7, gr7)

  # multiple reads at the same time
  # first pair returns no alignments
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAStringSet(c("AGTG", "ACTAGT", "ACTAAGT", "ACTAAAGT",
                               "AGGGTAAAGTCCATGGCCCCAATTTGTGTGTAG",
                               "AGGGTAAAAGTCCATGGCCCAATTTGTGTGTAG")),
    Biostrings::DNAString("CCCCCCCCCCCAGTGAAGTCAAACATGGAATTAGTGTGTTAA"),
    type = "global")
  read_ids <- c("5", "6", "6", "3", "5", "5", "5", "6", "6", "1", "2", "3",
                "4", "2", "3",
                "3", "4", "4", "5","5", "5", "5", "5", "5", "6", "6", "6")
  gr8 <- GenomicRanges::GRanges(
    seqnames = "test",
    ranges = IRanges::IRanges(c(29, 16, 29, 4, 4, 21, 40, 21, 40, 16, 36, 20,
                                20, 31, 1, 3, 13,
                                15, 1, 2, 3, 15, 33, 42, 14, 33, 42),
                              c(32, 18, 31, 15, 12, 23, 40, 23, 40, 42, 42, 42,
                                42, 31, 1, 3,
                                13, 15, 1, 2, 3, 15, 33, 42, 14, 33, 42)),
    strand = rep("+", 27),
    originally = c(rep("", 13), "T", "C", "C", "G", "G", "C",
                   "C", "C", "G", "A", "A", "T", "A", "A"),
    replacement = c("CCCC", "TAA", "CCC", rep("", 10), "C", "A", "T", "C", "A",
                    "A", "G", "G", "A", "T", "G", "G", "T", "G"),
    type = c(rep("insertion", 3), rep("deletion", 10), rep("mismatch", 14)),
    read_id = read_ids)
  names(gr8) <- read_ids
  test_gr8 <- getEventInfo(events, "test", 1)
  test_gr8$score <- NULL
  expect_identical(sort(test_gr8), sort(gr8))

  # minus strand
  strand(gr8) <- "-"
  ranges(gr8)[10:13] <- IRanges::IRanges(c(1, 1, 1, 1), c(11, 29, 11, 11))
  gr8$read_id[12:13] <- c("4", "6")
  names(gr8)[12:13] <- c("4", "6")
  test_gr8 <- getEventInfo(events, "test", 1, "-")
  test_gr8$score <- NULL
  expect_identical(sort(test_gr8), sort(gr8))

  # test shift     #10                       35     42
  s1 <- "ACTAGT"   #AC--------------------TAGT
  ampl <- "CCCCCCCCCCCAGTGAAGTCAAACATGGAATTAGTGTGTTAA"
  fwdPrPos <- stringr::str_locate(ampl, "CCAGT")[1, 1]
  events <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString(s1),
    Biostrings::subseq(ampl, start = fwdPrPos),
    type = "global")
  gr9 <- GenomicRanges::GRanges(
    seqnames = "test",
    ranges = IRanges::IRanges(c(12, 36, 10),
                              c(31, 42, 10)),
    strand = rep("+", 3),
    originally = c("", "", "C"),
    replacement = c("", "", "A"),
    type = c(rep("deletion", 2), "mismatch"),
    read_id = rep("1", 3))
  names(gr9) <- rep("1", 3)
  test_gr9 <- getEventInfo(events, "test", fwdPrPos)
  test_gr9$score <- NULL
  expect_identical(sort(test_gr9), sort(gr9))
})
