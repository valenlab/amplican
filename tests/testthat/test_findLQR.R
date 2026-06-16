library(amplican)
library(testthat)
library(data.table)
context("findLQR off-target filter")

fix <- function(case) {
  as.data.frame(fread(system.file("extdata", "findLQR", case, "events.csv",
                                  package = "amplican")))
}

n_culled <- function(l, aln) {
  if (!any(l)) return(0L)
  as.integer(sum(unique(as.data.table(aln)[l, ], by = "read_id")$counts))
}

test_that("findLQR disables itself at high editing (>25% candidate cull)", {
  # nlgn4a @ 90% editing: the edited majority is ~51% of reads, so the
  # candidate cluster trips the 25% safety cap and nothing is removed.
  aln <- fix("01_high_edit_falsepos")
  expect_warning(out <- findLQR(aln),
                 regexp = "off-target detection algorithm would remove")
  expect_false(any(out))
  expect_equal(n_culled(out, aln), 0L)
})

test_that("findLQR still removes a genuine small junk minority", {
  # Synthetic WT majority + small low-score junk cluster (~5% of reads):
  # below the 25% cap, so the junk is still removed.
  aln <- fix("06_genuine_junk")
  expect_silent(out <- findLQR(aln))
  expect_true(any(out))
  expect_gt(n_culled(out, aln), 0L)
})
