library(amplican)
library(testthat)
context("general helper functions")

config <- read.csv(system.file("extdata", "config.csv", package="amplican"),
                   stringsAsFactors = FALSE)
events <- read.csv(system.file("extdata", "results", "alignments_events.csv", package="amplican"),
                   stringsAsFactors = FALSE)

test_that("map_to_relative works correctly", {
  events <- map_to_relative(GRanges(events), config)
  expect_equal(start(events)[1], -32)
  expect_equal(end(events)[1], 51)
  expect_type(events$type, "character")
})
