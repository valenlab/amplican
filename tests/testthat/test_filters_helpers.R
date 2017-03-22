library(amplican)
context("Filter functions")

reads <- ShortRead::readFastq(system.file("extdata", "R1_001.fastq", package="amplican"))

test_that("goodAvgQuality filters below threshold", {
  expect_false(goodAvgQuality(reads, 30)[19]) # 19th read is bad average
  expect_equal(goodAvgQuality(reads, 35), c(rep(TRUE, 8), rep(FALSE, 4), rep(TRUE, 6), FALSE, FALSE))
  expect_true(all(goodAvgQuality(reads, 20)))
  expect_length(goodAvgQuality(reads), length(reads))
  expect_identical(goodAvgQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})

test_that("goodBaseQuality filters below threshold", {
  expect_false(goodBaseQuality(reads, 25)[19])
  expect_false(goodBaseQuality(reads, 25)[20])
  expect_equal(goodBaseQuality(reads, 35), c(rep(TRUE, 8), rep(FALSE, 4), rep(TRUE, 6), FALSE, FALSE))
  expect_true(all(goodBaseQuality(reads, 20)))
  expect_length(goodBaseQuality(reads), length(reads))
  expect_identical(goodBaseQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})

test_that("alphabetQuality filters N", {
  expect_identical(alphabetQuality(reads), c(rep(TRUE, 10), rep(FALSE, 2), rep(TRUE, 5), FALSE, TRUE, TRUE))
  expect_length(goodBaseQuality(reads), length(reads))
  expect_identical(goodBaseQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})
