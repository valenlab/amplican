library(amplican)
library(testthat)
context("plot helper functions")

config <- read.csv(system.file("extdata", "results", "config_summary.csv",
                               package="amplican"), stringsAsFactors = FALSE)
events <- read.csv(system.file("extdata", "results", "alignments",
                               "events_filtered_shifted_normalized.csv",
                               package="amplican"), stringsAsFactors = FALSE)

test_that("flipRanges returns whole integers", {
  flipped_events <- flipRanges(events, config)
  expect_false(any(flipped_events$start %% 1 != 0))
  expect_false(any(flipped_events$end %% 1 != 0))
})

test_that("Amplicon plot is a plot class.",{
  p <- plot_amplicon("aCACTGACTGACTAGACs", 1, 18)
  expect_is(p$layers[[1]], "ggproto")
  expect_equal(p$data$count, rep(1, length(p$data$count)))
  expect_identical(p$layers[[1]]$geom_params$parse, FALSE)
  expect_identical(p$layers[[1]]$geom_params$check_overlap, FALSE)
  expect_identical(p$layers[[1]]$geom_params$na.rm, FALSE)
})

test_that("Mismatch plot is a construct of grobs.",{
  p <- plot_mismatches(events, config, config$ID[1])
  expect_true(gtable::is.gtable(p))
})

test_that("When no mismatch to plot returns message.",{
  expect_match(plot_mismatches(events[events$type != "mismatch",],
                               config, config$ID[2]),
               "No mismatches to plot.")
})

test_that("Deletions plot is a construct of grobs.",{
  p <- plot_deletions(events, config, config$ID[1])
  expect_true(gtable::is.gtable(p))
})

test_that("When no deletions to plot returns text.",{
  expect_match(plot_deletions(events[events$type != "deletion",],
                              config, config$ID[2]),
               "No deletions to plot.")
})

test_that("Insertions plot is a construct of grobs.",{
  p <- plot_insertions(events, config, config$ID[3])
  expect_true(gtable::is.gtable(p))
})

test_that("When no insertions to plot returns text.",{
  expect_match(plot_insertions(events[events$type != "insertion",],
                               config, config$ID[2]),
               "No insertions to plot.")
})

test_that("Cuts plot is returning a plot.",{
  p <- plot_cuts(events, config, c(config$ID[1], config$ID[3]))
  expect_true(ggplot2::is.ggplot(p))
})

test_that("When no cuts to plot returns text.",{
  expect_match(plot_cuts(events[events$type != "deletion",],
                         config, config$ID[2]),
               "No cuts to plot.")
})
