library(amplican)
library(testthat)

test_that("all FASTQ modes and batch sizes produce the explicit pipeline truth", {
  folder <- tempfile("pipeline-modes-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  for (mode in c(0, 0.5, 1, 2)) {
    for (batch in c(1, 7, 1000)) {
      results <- file.path(folder, paste(mode, batch, sep = "-"))
      up_pipeline(fixture, results, fastqfiles = mode, batch_size = batch)
      output <- up_outputs(results)
      expected <- up_summary()
      expect_equal(up_rows(output$summary, names(expected)), up_rows(expected), info = paste(mode, batch))
      events <- output$events[consensus & overlaps]
      expect_equal(up_rows(events, names(up_events())), up_rows(up_events()), info = paste(mode, batch))
    }
  }
})

test_that("input order, partial resume and legacy alignments preserve numerical outputs", {
  folder <- tempfile("pipeline-resume-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  results <- file.path(folder, "baseline")
  up_pipeline(fixture, results)
  truth <- up_outputs(results)
  compare <- function(path) {
    output <- up_outputs(path)
    for (part in names(truth)) expect_equal(up_rows(output[[part]]), up_rows(truth[[part]]), info = paste(path, part))
  }
  up_pipeline(fixture, results)
  compare(results)
  stages <- list(
    alignments = "alignments/temp",
    raw = c("alignments/temp", "alignments/raw_events.csv", "alignments/experiment_data.rds",
            "barcode_reads_filters.csv"),
    filtered = c("alignments/temp", "alignments/raw_events.csv", "alignments/experiment_data.rds",
                 "alignments/events_filtered_shifted.csv", "config_summary.csv", "barcode_reads_filters.csv"))
  for (stage in names(stages)) {
    path <- file.path(folder, stage); dir.create(path)
    for (item in stages[[stage]]) {
      dir.create(dirname(file.path(path, item)), recursive = TRUE, showWarnings = FALSE)
      file.copy(file.path(results, item), dirname(file.path(path, item)), recursive = TRUE)
    }
    writeLines("interrupted", file.path(path, "alignments", "stale.temp"))
    up_pipeline(fixture, path)
    expect_false(file.exists(file.path(path, "alignments", "stale.temp")))
    compare(path)
  }
  legacy <- file.path(folder, "legacy"); dir.create(file.path(legacy, "alignments"), recursive = TRUE)
  object <- up_align(fixture)
  saveRDS(object, file.path(legacy, "alignments", "AlignmentsExperimentSet.rds"))
  up_pipeline(fixture, legacy)
  compare(legacy)
  write.csv(fixture$cfg[4:1, ], fixture$config, row.names = FALSE)
  for (id in names(fixture$reads)) {
    f <- ShortRead::readFastq(file.path(folder, paste0(id, "_1.fastq")))
    r <- ShortRead::readFastq(file.path(folder, paste0(id, "_2.fastq")))
    up_fastq(file.path(folder, paste0(id, "_1.fastq")), rev(as.character(ShortRead::sread(f))))
    up_fastq(file.path(folder, paste0(id, "_2.fastq")), rev(as.character(ShortRead::sread(r))))
  }
  reordered <- file.path(folder, "reordered")
  up_pipeline(fixture, reordered)
  compare(reordered)
})

test_that("both exported pipeline defaults apply their own normalization threshold", {
  folder <- tempfile("pipeline-defaults-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  results <- file.path(folder, "conservative")
  up_pipeline(fixture, results, pipeline = amplicanPipelineConservative)
  output <- up_outputs(results)
  expect_equal(up_rows(output$events[consensus & overlaps], names(up_events(TRUE))), up_rows(up_events(TRUE)))
  expect_equal(up_rows(output$summary, names(up_summary())), up_rows(up_summary()))
})

test_that("pipeline filter wiring removes whole reads and updates weighted denominators", {
  folder <- tempfile("pipeline-filter-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  fixture$cfg$Donor <- ""
  write.csv(fixture$cfg, fixture$config, row.names = FALSE)
  calls <- character()
  filter <- function(aln, ...) {
    calls <<- c(calls, unique(aln$seqnames))
    aln$type == "insertion" & aln$seqnames == "T"
  }
  testthat::local_mocked_bindings(findLQR = filter, .package = "amplican")
  results <- file.path(folder, "enabled")
  up_pipeline(fixture, results, event_filter = TRUE)
  output <- up_outputs(results)
  expect_setequal(calls, c("T", "R", "C"))
  expect_equal(output$summary$Low_Score, c(5, 0, 0, 0))
  expect_equal(output$summary$Reads_Filtered, c(16, 17, 10, 0))
  expect_false(any(output$events$seqnames == "T" & output$events$type == "insertion"))
  expect_true(any(output$events$seqnames == "R" & output$events$type == "insertion"))
  calls <- character()
  up_pipeline(fixture, file.path(folder, "disabled"), event_filter = FALSE)
  expect_length(calls, 0)
})

test_that("pipeline summaries exclude indels rejected by strict consensus", {
  folder <- tempfile("pipeline-consensus-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  cfg <- fixture$cfg[1, ]; cfg$Donor <- ""
  write.csv(cfg, fixture$config, row.names = FALSE)
  reverse <- fixture$reads$T
  reverse[reverse == fixture$insertion] <- fixture$reference
  up_fastq(file.path(folder, "T_2.fastq"), paste0("TTCC",
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reverse)))))
  results <- file.path(folder, "strict")
  up_pipeline(fixture, results, normalize = NULL, promiscuous_consensus = FALSE)
  output <- up_outputs(results)
  expect_equal(output$summary$Reads_Edited, 7)
  expect_equal(output$summary$Reads_In, 0)
  expect_equal(output$summary$Reads_Frameshifted, 0)
  expect_true(any(output$events$type == "insertion" & !output$events$consensus))
})

test_that("two actual workers preserve alignment, extraction and pipeline truth", {
  skip_on_os("windows")
  folder <- tempfile("pipeline-workers-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  fixture <- up_fixture(folder)
  previous <- BiocParallel::registered()
  on.exit(for (param in rev(previous)) BiocParallel::register(param), add = TRUE)
  param <- BiocParallel::MulticoreParam(workers = 2, fallback = FALSE)
  on.exit(BiocParallel::bpstop(param), add = TRUE)
  BiocParallel::register(param)
  workers <- unlist(BiocParallel::bplapply(1:2, function(i) Sys.getpid(), BPPARAM = param))
  expect_length(unique(workers), 2)
  expect_true(all(workers != Sys.getpid()))
  serial <- up_align(fixture)
  parallel <- up_align(fixture, use_parallel = TRUE)
  expect_equal(experimentData(parallel), experimentData(serial))
  expect_equal(readCounts(parallel), readCounts(serial))
  expect_equal(up_rows(extractEvents(parallel, use_parallel = TRUE)), up_rows(extractEvents(serial)))
  results <- file.path(folder, "parallel")
  expect_no_warning(up_pipeline(fixture, results, use_parallel = TRUE))
  output <- up_outputs(results)
  expect_equal(up_rows(output$summary, names(up_summary())), up_rows(up_summary()))
  expect_equal(up_rows(output$events[consensus & overlaps], names(up_events())), up_rows(up_events()))
})
