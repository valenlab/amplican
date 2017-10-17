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
  expect_true(all(amplicanOverlap(aln, cfgT, cut_buffer = 0) == aln$overlaps))
})

test_that("amplicanConsensus works as intended when not relative", {
  expect_true(all(amplicanConsensus(aln, cfgT) == aln$consensus))
})

aln <- amplicanMap(aln, cfgT)
aln <- as.data.frame(aln)

test_that("amplicanOverlap works as intended when relative", {
  expect_true(all(amplicanOverlap(aln, cfgT, cut_buffer = 0,
                                  relative = TRUE) == aln$overlaps))
})

test_that("amplicanConsensus works as intended when relative", {
  expect_true(all(amplicanConsensus(aln, cfgT) == aln$consensus))
})
