library(amplican)
library(testthat)
context("general helper functions")

config <- read.csv(system.file("extdata", "config.csv", package="amplican"),
                   stringsAsFactors = FALSE)
events <- read.csv(system.file("extdata", "results", "alignments",
                               "raw_events.csv", package="amplican"),
                   stringsAsFactors = FALSE)

test_that("amplicanMap works correctly", {
  events <- amplicanMap(
    GenomicRanges::GRanges(events[events$seqnames == "ID_1",]), config)
  expect_equal(GenomicRanges::start(events)[1], 42)
  expect_equal(GenomicRanges::end(events)[1], 61)
  expect_type(events$type, "character")
})
