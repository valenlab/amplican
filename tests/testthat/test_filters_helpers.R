library(amplican)
context("Filter functions")

reads <- ShortRead::readFastq(system.file("extdata", "R1_001.fastq", package="amplican"))

test_that("goodAvgQuality filters below threshold", {
  expect_false(goodAvgQuality(reads, 30)[8]) # 19th read is bad average
  expect_equal(goodAvgQuality(reads, 35), c(rep(TRUE, 7), rep(FALSE, 4), TRUE, rep(FALSE, 2), rep(TRUE, 6)))
  expect_true(all(goodAvgQuality(reads, 20)))
  expect_length(goodAvgQuality(reads), length(reads))
  expect_identical(goodAvgQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})

test_that("goodBaseQuality filters below threshold", {
  expect_false(goodBaseQuality(reads, 25)[13])
  expect_false(goodBaseQuality(reads, 25)[14])
  expect_equal(goodBaseQuality(reads, 35), c(rep(TRUE, 7), rep(FALSE, 4), TRUE, FALSE, FALSE, rep(TRUE, 6)))
  expect_true(all(goodBaseQuality(reads, 20)))
  expect_length(goodBaseQuality(reads), length(reads))
  expect_identical(goodBaseQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})

test_that("alphabetQuality filters N", {
  expect_identical(alphabetQuality(reads), c(rep(TRUE, 9), rep(FALSE, 3), rep(TRUE, 8)))
  expect_length(goodBaseQuality(reads), length(reads))
  expect_identical(goodBaseQuality(c(TRUE, FALSE, TRUE)), c(TRUE, FALSE, TRUE))
})
