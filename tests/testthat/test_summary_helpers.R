library(amplican)
library(testthat)
context("summary helper functions")

# 3 - types ins, del, mm
# 4 - situations (in agreement, disagreement, empty vs event, partial vs event)
# 2 - for each strand fwd and reverse
# 2 - overlapping case and non-overlapping
# 2 - ID_1 and ID_2

aln <- fread(system.file("test_data", "test_aln.csv", package = "amplican"))
cfgT <- fread(system.file("test_data", "test_cfg.csv", package = "amplican"))

test_that("amplicanOverlap works as intended when not relative", {
  expect_true(all(amplicanOverlap(aln, cfgT, cut_buffer = 0,
                                  relative = FALSE) == aln$overlaps))
})

test_that("amplicanConsensus works as intended when not relative", {
  expect_true(all(amplicanConsensus(aln, cfgT,
                                    promiscuous = FALSE) == aln$consensus))
})

aln <- amplicanMap(aln, cfgT)
aln <- as.data.frame(aln)

test_that("amplicanOverlap works as intended when relative", {
  expect_true(all(amplicanOverlap(aln, cfgT, cut_buffer = 0,
                                  relative = TRUE) == aln$overlaps))
})

test_that("amplicanConsensus works as intended when relative", {
  expect_true(all(amplicanConsensus(aln, cfgT,
                                    promiscuous = FALSE) == aln$consensus))
})

test_that("buffered overlap uses the same origin in absolute and relative coordinates", {
  config <- data.frame(ID = c("A", "B", "whole"), Direction = 0L,
    Amplicon = c("aaaaaGGGGGaaaaa", "aaaCCCaaaaaTTaaaa", "aaaaaaaaaa"),
    ampl_len = c(15L, 17L, 10L), stringsAsFactors = FALSE)
  events <- data.frame(seqnames = c(rep("A", 6), rep("B", 6), rep("whole", 4)),
    start = c(3, 4, 5, 6, 12, 13, 2, 8, 9, 10, 15, 16, 0, 1, 10, 11), width = 1L)
  events$end <- events$start
  expected <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)
  expect_identical(suppressMessages(amplicanOverlap(events, config, cut_buffer = 2)), expected)
  origins <- c(A = 6L, B = 4L, whole = 1L)[events$seqnames]
  relative_events <- transform(events, start = start - origins, end = end - origins)
  expect_identical(suppressMessages(amplicanOverlap(relative_events, config, cut_buffer = 2, relative = TRUE)), expected)
  expect_identical(amplicanOverlap(events[FALSE, ], config), logical())
  expect_identical(suppressMessages(amplicanOverlap(events[16:1, ], config[3:1, ], cut_buffer = 2)), rev(expected))
})

test_that("consensus decisions preserve score ties, modes, and experiment identity", {
  config <- data.table::data.table(ID = c("A", "B"),
    Amplicon = paste0(strrep("a", 40), strrep("G", 20), strrep("a", 140)),
    Direction = 0L, fwdPrPosEnd = 10L, rvePrPos = 180L)
  pair <- function(start = c(50, 50), width = c(1, 1), replacement = c("C", "C"),
                   score = c(10, 20), id = c("A", "A"), overlaps = c(TRUE, TRUE),
                   type = c("mismatch", "mismatch")) {
    data.table::data.table(seqnames = id, read_id = 1L, start = start,
      end = start + width - 1, width = width, replacement = replacement,
      originally = "A", type = type, strand = c("+", "-"), score = score,
      counts = 1L, overlaps = overlaps)
  }
  cases <- list(
    agreement = list(pair(), c(TRUE, FALSE), c(TRUE, FALSE)),
    forward_wins = list(pair(replacement = c("C", "G"), score = c(30, 10)), c(TRUE, FALSE), c(TRUE, FALSE)),
    reverse_wins = list(pair(replacement = c("C", "G")), c(FALSE, TRUE), c(FALSE, TRUE)),
    tie = list(pair(replacement = c("C", "G"), score = c(10, 10)), c(TRUE, FALSE), c(TRUE, FALSE)),
    unconfirmed = list(pair(start = c(50, 55)), c(FALSE, FALSE), c(TRUE, TRUE)),
    experiments = list(pair(id = c("A", "B")), c(FALSE, FALSE), c(TRUE, TRUE)),
    outside_cut = list(pair(start = c(80, 85), overlaps = c(FALSE, FALSE)), c(FALSE, FALSE), c(FALSE, FALSE)),
    truncated = list(pair(start = c(15, 1), width = c(36, 20), replacement = c("", ""),
                          type = c("deletion", "deletion"), overlaps = c(TRUE, FALSE)), c(TRUE, FALSE), c(TRUE, FALSE)))
  for (name in names(cases)) {
    case <- cases[[name]]
    expect_identical(amplicanConsensus(data.table::copy(case[[1]]), config, promiscuous = FALSE), case[[2]], info = name)
    expect_identical(amplicanConsensus(data.table::copy(case[[1]]), config, promiscuous = TRUE), case[[3]], info = name)
    expect_identical(amplicanConsensus(data.table::copy(case[[1]]), config), case[[3]], info = name)
  }
  empty <- pair()[0]
  expect_identical(amplicanConsensus(empty, config), logical())
  single <- pair()[1]
  expect_identical(amplicanConsensus(data.table::copy(single), config, promiscuous = FALSE), FALSE)
  expect_identical(amplicanConsensus(data.table::copy(single), config, promiscuous = TRUE), TRUE)
})

test_that("overlap hits retain original indices and do not mix experiment/read pairs", {
  forward <- data.table::data.table(seqnames = c("A", "B", "A", "A"),
    read_id = c(1L, 1L, 2L, 1L), start = c(10L, 10L, 10L, 30L), end = c(20L, 20L, 20L, 40L))
  reverse <- data.table::data.table(seqnames = c("A", "A", "B", "A"),
    read_id = c(1L, 2L, 1L, 1L), start = c(40L, 20L, 15L, 21L), end = c(45L, 25L, 15L, 29L))
  hits <- amplican:::getHits(forward, reverse)
  pairs <- data.frame(from = S4Vectors::from(hits), to = S4Vectors::to(hits))
  pairs <- pairs[order(pairs$from, pairs$to), ]
  rownames(pairs) <- NULL
  expect_identical(pairs, data.frame(from = c(2L, 3L, 4L), to = c(3L, 2L, 1L)))
  expect_identical(amplican:::getHits(forward[0], reverse), S4Vectors::Hits())
  expect_identical(amplican:::getHits(forward, reverse[0]), S4Vectors::Hits())
})

test_that("two clipped reads resolve their shared disagreement by alignment score", {
  config <- data.table::data.table(ID = "A", Direction = 0L,
    Amplicon = paste0(strrep("a", 40), strrep("G", 20), strrep("a", 140)),
    fwdPrPosEnd = 10L, rvePrPos = 180L)
  events <- data.table::data.table(seqnames = "A", read_id = 1L,
    start = c(55, 50, 1, 50), end = c(200, 50, 45, 50),
    width = c(146, 1, 45, 1), strand = c("+", "+", "-", "-"),
    originally = c("", "A", "", "A"), replacement = c("", "C", "", "G"),
    type = c("deletion", "mismatch", "deletion", "mismatch"),
    score = c(30, 30, 10, 10), counts = 1L, overlaps = TRUE)
  for (mode in c(FALSE, TRUE)) {
    expect_identical(amplicanConsensus(data.table::copy(events), config, promiscuous = mode),
                     c(FALSE, TRUE, FALSE, FALSE))
  }
  reverse_clipped <- data.table::copy(events[c(2, 3, 4)])
  reverse_clipped$score <- c(10, 30, 30)
  forward_clipped <- events[c(1, 2, 4)]
  for (one in list(reverse_clipped, forward_clipped)) {
    expect_identical(amplicanConsensus(data.table::copy(one), config, promiscuous = FALSE), rep(FALSE, 3))
  }
  expect_identical(amplicanConsensus(data.table::copy(reverse_clipped), config), c(TRUE, FALSE, FALSE))
  expect_identical(amplicanConsensus(data.table::copy(forward_clipped), config), c(FALSE, FALSE, TRUE))
})
