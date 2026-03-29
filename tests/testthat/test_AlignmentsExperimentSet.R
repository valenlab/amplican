library(amplican)
library(testthat)
context("AlignmentsExperimentSet S4 Class Operations")

config <- system.file("extdata", "config.csv", package = "amplican")
fastq_folder <- system.file("extdata", package = "amplican")
results_folder <- tempdir()

# Generate real AlignmentsExperimentSet objects to test
dir.create(file.path(results_folder, "alignments"), showWarnings = FALSE)
suppressWarnings(
  amplicanAlign(config, fastq_folder, temp_folder = file.path(results_folder, "alignments"), 
                fastqfiles = 0, primer_mismatch = 0)
)
aln1 <- readRDS(file.path(results_folder, "alignments", "barcode_1_aln.rds"))
aln2 <- readRDS(file.path(results_folder, "alignments", "barcode_2_aln.rds"))

test_that("Initialization and Getters work", {
  expect_s4_class(aln1, "AlignmentsExperimentSet")
  # barcode_1 has two experiments in config.csv: ID_1 and ID_2 (Wait, let's just check > 0)
  expect_true(length(fwdReads(aln1)) > 0)
  expect_true(length(readCounts(aln1)) > 0)
})

test_that("Setters trigger validity checking correctly", {
  # Modify using an invalid state
  expect_error({
    readCounts(aln1) <- list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) # Intentionally invalid length
  })
})

test_that("Subsetting drops unassignedData and barcodeData when taking 1 element", {
  aln_sub <- aln1[1] 
  expect_s4_class(aln_sub, "AlignmentsExperimentSet")
  expect_equal(length(fwdReads(aln_sub)), 1)
  expect_null(aln_sub@unassignedData)
  expect_null(aln_sub@barcodeData)
})

test_that("Concatenation c() works correctly on real objects", {
  aln_combined <- c(aln1, aln2)
  expect_s4_class(aln_combined, "AlignmentsExperimentSet")
  expect_equal(length(fwdReads(aln_combined)), length(fwdReads(aln1)) + length(fwdReads(aln2)))
  expect_equal(nrow(aln_combined@experimentData), nrow(aln1@experimentData) + nrow(aln2@experimentData))
  # Should drop barcode and unassigned when combining according to class definition or maybe it merges.
  # Let's check what it does actually. `c` method says: 
  # barcodeData = as.data.frame(data.table::rbindlist(lapply(args, barcodeData), fill=TRUE))
  expect_equal(nrow(barcodeData(aln_combined)), 2)
})
