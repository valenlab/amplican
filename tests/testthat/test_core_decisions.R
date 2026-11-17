library(amplican)
library(testthat)

# Small explicit tables keep decision tests independent of alignment extraction.
decision_cfg <- function(ids = c("A", "B")) {
  data.table::data.table(
    ID = ids, Amplicon = paste0(strrep("a", 40), strrep("G", 20), strrep("a", 140)),
    Direction = 0L, fwdPrPosEnd = 10L, rvePrPos = 180L,
    Forward_Primer = strrep("A", 10), Reverse_Primer = strrep("T", 10))
}

decision_events <- function(id = "A", read_id = 1L, start = 30L, width = 1L,
                            type = "mismatch", strand = "+", score = 100,
                            counts = 1L, originally = "A", replacement = "C") {
  data.table::data.table(
    seqnames = id, start = start, end = start + width - 1L, width = width,
    strand = strand, originally = originally, replacement = replacement,
    type = type, read_id = read_id, score = score, counts = counts,
    consensus = TRUE, overlaps = TRUE)
}

decision_rows <- function(events) {
  events <- as.data.frame(events)
  events <- events[order(events$seqnames, events$read_id, events$start, events$type), ]
  rownames(events) <- NULL
  events
}

test_that("primer and dimer filters respect boundaries and experiment identity", {
  cfg <- decision_cfg()
  cfg$rvePrPos[2] <- 160L
  events <- decision_events(
    id = c("A", "A", "B", rep("A", 7)), read_id = c(1, 1, 1, 2, 2, 3, 4, 5, 6, 7),
    start = c(20, 30, 30, 9, 50, 10, 180, 181, 180, 181),
    width = c(151, 1, 1, 1, 1, 150, 4, 1, 1, 1),
    type = c("deletion", rep("mismatch", 4), "deletion", "insertion", "insertion", "mismatch", "mismatch"))
  expect_identical(findEOP(events, cfg), c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))
  expect_identical(findPD(events, cfg, PRIMER_DIMER = 30), c(TRUE, rep(FALSE, 9)))
  expected <- events[c(3, 5, 6, 7, 9)]
  expect_equal(decision_rows(amplicanFilter(data.table::copy(events), cfg, PRIMER_DIMER = 30)),
               decision_rows(expected))
  expect_equal(decision_rows(amplicanFilter(data.table::copy(events[0]), cfg, PRIMER_DIMER = 30)),
               decision_rows(events[0]))
  expect_identical(findEOP(events[0], cfg), logical())
  expect_identical(findPD(events[0], cfg), logical())
  valid <- decision_events(id = c("A", "B"), read_id = c(77, 88))
  expect_equal(decision_rows(amplicanFilter(data.table::copy(valid), cfg, PRIMER_DIMER = 30)), decision_rows(valid))

  # Relative coordinates require a negative event to identify their coordinate system.
  relative <- data.table::copy(events)
  relative[, `:=`(start = start - 41L, end = end - 41L)]
  expect_identical(findEOP(relative, cfg), findEOP(events, cfg))
  expect_identical(findEOP(relative[1], cfg), FALSE)
  at_origin <- decision_events(start = 0L, width = 2L, type = "insertion", originally = "", replacement = "CC")
  expect_identical(findEOP(at_origin, cfg), FALSE)
  cfg_missing <- data.table::copy(cfg)
  cfg_missing[, `:=`(fwdPrPosEnd = NA_integer_, rvePrPos = NA_integer_)]
  expect_identical(findEOP(decision_events(start = c(1, 200, 201)), cfg_missing), c(FALSE, FALSE, TRUE))
})

test_that("normalization pools controls and preserves threshold equality and group isolation", {
  cfg <- data.table::data.table(
    ID = c("C1", "C2", "T", "other_group", "other_guide"),
    Control = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    Reads_Filtered = c(100, 300, 500, 100, 100),
    Group = c("one", "one", "one", "two", "one"),
    guideRNA = c("g", "g", "g", "g", "h"))
  events <- decision_events(
    id = c("C1", "C1", "C1", "C2", "C1", rep("T", 5), "other_group", "other_guide"),
    read_id = seq_len(12), start = c(30, 40, 50, 50, 60, 30, 40, 50, 60, 70, 50, 50),
    counts = c(1, 4, 2, 3, 10, 7, 11, 13, 17, 19, 23, 29),
    replacement = c(rep("C", 4), "G", rep("C", 7)), score = seq_len(12),
    strand = rep(c("+", "-"), 6))
  normalize <- function(x = events, config = cfg, ...) {
    amplicanNormalize(data.table::copy(x), data.table::copy(config), ...)
  }
  # Control total is 400: 1/400 is below, 4/400 equals, and 5/400 exceeds 1%.
  expect_equal(decision_rows(normalize()), decision_rows(events[-8]))
  expect_equal(decision_rows(normalize(events[12:1], cfg[5:1])), decision_rows(events[-8]))
  expect_equal(decision_rows(normalize(add = character())), decision_rows(events[-c(8, 11, 12)]))
  expect_identical(normalize(add = NULL), events)
  expect_equal(decision_rows(normalize(min_freq = 0)), decision_rows(events[-c(6, 7, 8)]))
  expect_equal(decision_rows(normalize(min_freq = 1)), decision_rows(events))
  expect_equal(decision_rows(normalize(events[0])), decision_rows(events[0]))
  expect_equal(decision_rows(normalize(events[1:5])), decision_rows(events[1:5]))
  extra_config <- data.table::rbindlist(list(cfg, data.table::data.table(
    ID = "C3", Control = TRUE, Reads_Filtered = 600, Group = "two", guideRNA = "g")))
  extra_events <- data.table::rbindlist(list(events,
    decision_events(id = "C3", read_id = 13, start = 50, counts = 5)))
  expect_equal(decision_rows(normalize(extra_events, extra_config)), decision_rows(extra_events[-8]))
  expect_equal(decision_rows(normalize(extra_events, extra_config, add = character())), decision_rows(extra_events))
  duplicate_controls <- data.table::copy(events[1:5])
  duplicate_controls$strand <- ifelse(duplicate_controls$strand == "+", "-", "+")
  duplicate_controls$consensus <- FALSE
  paired <- data.table::rbindlist(list(events, duplicate_controls))
  expect_equal(decision_rows(normalize(paired)), decision_rows(paired[-8]))
  without_ids <- events[, setdiff(names(events), "read_id"), with = FALSE]
  # Aggregated event tables without read IDs retain their documented interface.
  result <- normalize(without_ids)
  expect_equal(result[order(seqnames, start, counts)], without_ids[-8][order(seqnames, start, counts)])

  changed <- data.table::copy(events)
  changed$originally[8] <- "T"
  expect_equal(decision_rows(normalize(changed)), decision_rows(changed))
  changed <- data.table::copy(events)
  changed$type[8] <- "insertion"
  changed$originally[8] <- ""
  expect_equal(decision_rows(normalize(changed)), decision_rows(changed))
  no_controls <- data.table::copy(cfg)
  no_controls$Control <- FALSE
  expect_warning(unchanged <- normalize(config = no_controls), "no TRUE/1")
  expect_identical(unchanged, events)
  before <- data.table::copy(cfg)
  invisible(amplicanNormalize(data.table::copy(events), cfg))
  expect_identical(cfg, before)
})

test_that("quality decisions include exact thresholds and survive batching", {
  reads <- ShortRead::ShortReadQ(
    sread = Biostrings::DNAStringSet(c("ACGT", "ACGT", "ACGT", "ACGT", "ANGT")),
    quality = ShortRead::FastqQuality(c("5555", "4444", "6666", "4655", "5555")),
    id = Biostrings::BStringSet(paste0("read", 1:5)))
  expected <- list(goodBaseQuality = c(TRUE, FALSE, TRUE, FALSE, TRUE),
                   goodAvgQuality = c(TRUE, FALSE, TRUE, TRUE, TRUE),
                   alphabetQuality = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  for (name in names(expected)) {
    fn <- get(name, asNamespace("amplican"))
    for (batch in c(NA_integer_, 1L, 2L, 5L, 10L)) {
      args <- list(reads = reads, batch_size = batch)
      if (name == "goodBaseQuality") args$min <- 20
      if (name == "goodAvgQuality") args$avg <- 20
      expect_identical(do.call(fn, args), expected[[name]], info = paste(name, "batch", batch))
    }
    expect_identical(fn(c(TRUE, FALSE)), c(TRUE, FALSE))
    expect_identical(fn(logical()), logical())
  }
})

lqr_decision_fixture <- function(populations = c(96, 18, 6), events_per_read = c(8L, 12L, 30L),
                                 alignment_scores = c(500, 250, 10), last_positions = c(200L, 200L, 200L)) {
  group <- rep(1:3, populations)
  n <- events_per_read[group]
  positions <- unlist(lapply(seq_along(n), function(i) c(seq_len(n[i] - 1L), last_positions[group[i]])))
  decision_events(read_id = rep(seq_along(group), n), start = positions,
                  score = rep(alignment_scores[group], n))
}

test_that("off-target filtering respects weighted caps and fixed-seed decisions", {
  events <- lqr_decision_fixture()
  expected <- events$read_id >= 115L
  expect_identical(findLQR(data.table::copy(events), seed = 0), expected)
  expect_identical(findLQR(data.table::copy(events), seed = 7, max_remove_filterLQR = 0.05), expected)
  expect_identical(findLQR(data.table::copy(events), max_remove_filterLQR = 0.051), expected)
  expect_warning(out <- findLQR(data.table::copy(events), max_remove_filterLQR = 0.049), "more than 5%.*sample A")
  expect_identical(out, rep(FALSE, nrow(events)))
  weighted <- data.table::copy(events)
  weighted$counts[expected] <- 50L
  expect_warning(out <- findLQR(weighted), "more than 25%.*sample A")
  expect_identical(out, rep(FALSE, nrow(events)))
  reversed <- events[nrow(events):1]
  expect_identical(findLQR(reversed), rev(expected))
  expect_identical(findLQR(events[1:999]), rep(FALSE, 999))
  small <- decision_events(read_id = rep(1:50, each = 30), start = rep(1:30, 50))
  expect_identical(findLQR(small), rep(FALSE, 1500))
  flat <- decision_events(read_id = rep(1:100, each = 10), start = rep(1:10, 100))
  expect_identical(findLQR(flat), rep(FALSE, 1000))
  boundary <- lqr_decision_fixture(c(80, 15, 5), c(8L, 12L, 36L))
  expect_identical(findLQR(boundary), boundary$read_id >= 96L)
  paired_junk <- data.table::copy(events[expected])
  paired_junk$strand <- "-"
  mixed_strands <- data.table::rbindlist(list(events, paired_junk))
  expect_identical(findLQR(mixed_strands, max_remove_filterLQR = 0.075),
                   mixed_strands$read_id >= 115L)
  varied <- lqr_decision_fixture(events_per_read = c(8L, 16L, 30L), alignment_scores = c(500, 250, 200))
  expect_identical(findLQR(varied), varied$read_id >= 115L)
  for (ids in list(97:114, 115:120)) {
    extra <- data.table::copy(varied[read_id %in% ids])
    extra$strand <- "-"
    duplicated <- data.table::rbindlist(list(varied, extra))
    expect_identical(findLQR(duplicated), duplicated$read_id >= 115L)
  }
  lengths <- lqr_decision_fixture(events_per_read = c(8L, 16L, 60L),
    alignment_scores = c(500, 250, 200), last_positions = c(200L, 300L, 100L))
  expect_identical(findLQR(lengths), lengths$read_id >= 115L)
  lengths <- lqr_decision_fixture(events_per_read = c(8L, 16L, 30L),
    alignment_scores = c(500, 250, 200), last_positions = c(200L, 300L, 100L))
  expect_identical(findLQR(lengths), lengths$read_id >= 115L)
  cfg <- decision_cfg()
  cfg[, `:=`(fwdPrPosEnd = 1L, rvePrPos = 200L)]
  other <- data.table::copy(events)
  other$seqnames <- "B"
  other$read_id <- ifelse(other$read_id >= 115L, other$read_id - 114L, other$read_id + 6L)
  mixed <- data.table::rbindlist(list(events, other))
  expected_events <- mixed[!((seqnames == "A" & read_id >= 115L) | (seqnames == "B" & read_id <= 6L))]
  expect_equal(decision_rows(amplicanFilter(data.table::copy(mixed), cfg, PRIMER_DIMER = 30)),
               decision_rows(expected_events))
})

test_that("off-target filtering prefers the simpler model on a silhouette tie", {
  events <- lqr_decision_fixture()
  tied_clara <- function(x, k, ...) {
    structure(list(
      silinfo = list(widths = matrix(c(1, 2, 0.5), nrow = 1,
        dimnames = list(NULL, c("cluster", "neighbor", "sil_width")))),
      medoids = matrix(c(1, 0, 0.5, 0.5, 0, 1), ncol = 2, byrow = TRUE),
      clustering = ifelse(x[, 2] > 0.9, 3L, ifelse(x[, 2] > 0.1, 2L, 1L))),
      class = "partition")
  }
  testthat::with_mocked_bindings({
    expect_identical(findLQR(events), rep(FALSE, nrow(events)))
  }, clara = tied_clara, .package = "cluster")
})

test_that("off-target sampling uses the documented R-seed backend", {
  events <- lqr_decision_fixture()
  original_clara <- cluster::clara
  seeded_clara <- function(..., rngR = FALSE) {
    # Without this backend contract the public seed argument does not control CLARA.
    expect_true(rngR)
    original_clara(..., rngR = rngR)
  }
  testthat::with_mocked_bindings({
    expect_identical(findLQR(events, seed = 7), events$read_id >= 115L)
  }, clara = seeded_clara, .package = "cluster")
})

test_that("identical mate evidence leaves per-read clustering inputs unchanged", {
  events <- lqr_decision_fixture(events_per_read = c(8L, 16L, 30L), alignment_scores = c(500, 250, 200))
  mate <- data.table::copy(events[read_id %in% c(1:8, 97:114, 115:120)])
  mate$strand <- "-"
  captured <- list()
  original_clara <- cluster::clara
  capture_clara <- function(x, ...) {
    captured[[length(captured) + 1L]] <<- x
    original_clara(x, ...)
  }
  testthat::with_mocked_bindings({
    invisible(findLQR(data.table::copy(events)))
    invisible(findLQR(data.table::rbindlist(list(events, mate))))
  }, clara = capture_clara, .package = "cluster")
  expect_identical(captured[[1]], captured[[3]])
})

test_that("off-target selection uses both coordinates of the cluster distance", {
  sizes <- c(rep(10L, 10), rep(70L, 10), rep(100L, 80), 1L)
  positions <- unlist(lapply(sizes, function(n) c(seq_len(n - 1L), 1000L)))
  events <- decision_events(read_id = rep(seq_along(sizes), sizes), start = positions,
    score = rep(c(rep(0, 10), rep(500, 10), rep(1000, 81)), sizes))
  partition <- function(x, k, ...) {
    # Three observed medoids are (0, .091), (.5, .697), and (1, 1).
    # Their distances from (1, 0) are about 1.004, .858, and 1.
    widths <- matrix(c(1, 2, if (k == 3) 0.8 else 0.2), nrow = 1,
      dimnames = list(NULL, c("cluster", "neighbor", "sil_width")))
    structure(list(silinfo = list(widths = widths),
      medoids = x[if (k == 3) c(1, 11, 21) else c(1, 21), , drop = FALSE],
      clustering = if (k == 3) c(rep(1L, 10), rep(2L, 10), rep(3L, 81)) else c(rep(1L, 20), rep(2L, 81))),
      class = "partition")
  }
  testthat::with_mocked_bindings({
    expect_identical(findLQR(events), events$read_id <= 10L)
  }, clara = partition, .package = "cluster")
})
