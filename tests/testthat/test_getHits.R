library(amplican)
library(testthat)
library(data.table)
context("getHits foverlaps implementation")

# Reference implementation using GRanges (the old approach)
getHits_granges <- function(aln_fwd, aln_rve) {
  if (nrow(aln_fwd) == 0 | nrow(aln_rve) == 0) return(S4Vectors::Hits())
  suppressWarnings(GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(
      seqnames = paste0(aln_fwd$seqnames, "_", aln_fwd$read_id),
      ranges = IRanges::IRanges(start = aln_fwd$start,
                                end = aln_fwd$end),
      strand = "*"),
    GenomicRanges::GRanges(
      seqnames = paste0(aln_rve$seqnames, "_", aln_rve$read_id),
      IRanges::IRanges(start = aln_rve$start,
                       end = aln_rve$end),
      strand = "*"),
    type = "any", select = "all"))
}

# Helper to compare Hits objects regardless of ordering
hits_equal <- function(h1, h2) {
  if (length(h1) != length(h2)) return(FALSE)
  if (length(h1) == 0) return(TRUE)
  p1 <- data.table(from = S4Vectors::from(h1), to = S4Vectors::to(h1))
  p2 <- data.table(from = S4Vectors::from(h2), to = S4Vectors::to(h2))
  setkey(p1, from, to)
  setkey(p2, from, to)
  identical(p1, p2)
}

make_aln <- function(seqnames, read_id, start, end) {
  data.table(seqnames = seqnames, read_id = read_id,
             start = as.integer(start), end = as.integer(end))
}

test_that("getHits: exact matching intervals", {
  fwd <- make_aln("ID_1", "r1", 10, 20)
  rve <- make_aln("ID_1", "r1", 10, 20)
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 1)
})

test_that("getHits: partial overlap", {
  fwd <- make_aln("ID_1", "r1", 10, 20)
  rve <- make_aln("ID_1", "r1", 15, 25)
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 1)
})

test_that("getHits: no overlap same group", {
  fwd <- make_aln("ID_1", "r1", 10, 20)
  rve <- make_aln("ID_1", "r1", 30, 40)
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 0)
})

test_that("getHits: overlapping coords but different group", {
  fwd <- make_aln("ID_1", "r1", 10, 20)
  rve <- make_aln("ID_1", "r2", 10, 20)
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 0)
})

test_that("getHits: overlapping coords but different seqnames", {
  fwd <- make_aln("ID_1", "r1", 10, 20)
  rve <- make_aln("ID_2", "r1", 10, 20)
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 0)
})

test_that("getHits: empty inputs", {
  empty <- make_aln(character(0), character(0), integer(0), integer(0))
  fwd <- make_aln("ID_1", "r1", 10, 20)
  expect_equal(length(amplican:::getHits(empty, fwd)), 0)
  expect_equal(length(amplican:::getHits(fwd, empty)), 0)
  expect_equal(length(amplican:::getHits(empty, empty)), 0)
})

test_that("getHits: multiple events per read", {
  fwd <- make_aln(c("ID_1", "ID_1"), c("r1", "r1"), c(10, 30), c(20, 40))
  rve <- make_aln(c("ID_1", "ID_1"), c("r1", "r1"), c(15, 50), c(25, 60))
  # fwd[1] overlaps rve[1], fwd[2] does not overlap rve[2]
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 1)
})

test_that("getHits: one fwd overlaps multiple rve", {
  fwd <- make_aln("ID_1", "r1", 10, 50)
  rve <- make_aln(c("ID_1", "ID_1"), c("r1", "r1"), c(15, 40), c(25, 55))
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 2)
})

test_that("getHits: mixed groups with overlaps", {
  fwd <- make_aln(
    c("ID_1", "ID_1", "ID_2"),
    c("r1",   "r2",   "r1"),
    c(10,     10,     10),
    c(20,     20,     20))
  rve <- make_aln(
    c("ID_1", "ID_1", "ID_2"),
    c("r1",   "r2",   "r1"),
    c(15,     30,     12),
    c(25,     40,     18))
  # ID_1/r1: overlap; ID_1/r2: no overlap; ID_2/r1: overlap
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 2)
})

test_that("getHits: scale test with many unique groups", {
  n <- 500
  fwd <- make_aln(
    paste0("ID_", seq_len(n)),
    paste0("r_", seq_len(n)),
    rep(10L, n),
    rep(20L, n))
  rve <- make_aln(
    paste0("ID_", seq_len(n)),
    paste0("r_", seq_len(n)),
    rep(15L, n),
    rep(25L, n))
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), n)
})

test_that("getHits: from/to indices are correct", {
  fwd <- make_aln(c("ID_1", "ID_2"), c("r1", "r1"), c(10, 100), c(20, 200))
  rve <- make_aln(c("ID_2", "ID_1"), c("r1", "r1"), c(150, 5), c(250, 25))
  # fwd[1](ID_1/r1 10-20) overlaps rve[2](ID_1/r1 5-25) → from=1, to=2
  # fwd[2](ID_2/r1 100-200) overlaps rve[1](ID_2/r1 150-250) → from=2, to=1
  h_new <- amplican:::getHits(fwd, rve)
  h_ref <- getHits_granges(fwd, rve)
  expect_true(hits_equal(h_new, h_ref))
  expect_equal(length(h_new), 2)
})
