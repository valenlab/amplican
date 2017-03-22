library(amplican)
context("plot helper functions")

config <- read.csv(system.file("extdata", "config.csv", package="amplican"),
                   stringsAsFactors = FALSE)
events <- read.csv(system.file("extdata", "results", "alignments_events.csv", package="amplican"),
                   stringsAsFactors = FALSE)

test_that("flipRanges returns whole integers", {
  flipped_events <- flipRanges(events, config)
  expect_false(any(flipped_events$start %% 1 != 0))
  expect_false(any(flipped_events$end %% 1 != 0))
})

test_that("Amplicon plot is a plot class.",{
  p <- amplican_plot_amplicon("aCACTGACTGACTAGACs")
  expect_is(p$layers[[1]], "ggproto")
  expect_equal(p$data$count, rep(1, length(p$data$count)))
  expect_identical(p$layers[[1]]$geom_params$parse, FALSE)
  expect_identical(p$layers[[1]]$geom_params$check_overlap, FALSE)
  expect_identical(p$layers[[1]]$geom_params$na.rm, FALSE)
})

grobs_of_3_construct <- function(pp) {
  expect_is(pp, "PlotList")
  expect_equal(length(pp), 3)

  expect_identical(pp[[1]]$layers[[1]]$geom_params$width, NULL)
  expect_identical(pp[[1]]$layers[[1]]$geom_params$na.rm, FALSE)
  expect_identical(pp[[3]]$layers[[1]]$geom_params$width, NULL)
  expect_identical(pp[[3]]$layers[[1]]$geom_params$na.rm, FALSE)
}

test_that("Mismatch plot is a construct of grobs.",{
  p <- amplican_plot_mismatches(events, config, config$ID[1])
  grobs_of_3_construct(slot(p, "grobs"))
})

test_that("When no mismatch to plot returns message.",{
  expect_message(amplican_plot_mismatches(events[events$type != "mismatch",], config, config$ID[2]))
})

test_that("Deletions plot is a construct of grobs.",{
  p <- amplican_plot_deletions(events, config, config$ID[1])
  grobs_of_3_construct(slot(p, "grobs"))
})

test_that("When no deletions to plot returns warning.",{
  expect_message(amplican_plot_deletions(events[events$type != "deletion",], config, config$ID[2]))
})

test_that("Insertions plot is a construct of grobs.",{
  p <- amplican_plot_insertions(events, config, config$ID[3])
  grobs_of_3_construct(slot(p, "grobs"))
})

test_that("When no insertions to plot returns warning.",{
  expect_message(amplican_plot_insertions(events[events$type != "insertion",], config, config$ID[2]))
})

test_that("Cuts plot is returning a plot.",{
  p <- amplican_plot_cuts(events, config, c(config$ID[1], config$ID[3]))
  expect_is(p$layers[[1]], "ggproto")
  expect_equal(class(p$data), "waiver")
  expect_identical(as.character(p$layers[[1]]$mapping[["colour"]]), "seqnames")
})

test_that("When no cuts to plot returns warning.",{
  expect_message(amplican_plot_cuts(events[events$type != "deletion",], config, config$ID[2]))
})
