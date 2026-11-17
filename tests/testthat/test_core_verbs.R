library(amplican)
library(testthat)
context("Core Pipeline Verbs Operations")

# Use real test data rather than generating mock events
aln <- data.table::fread(system.file("test_data", "test_aln.csv", package="amplican"))
cfgT <- data.table::fread(system.file("test_data", "test_cfg.csv", package="amplican"))

# Provide missing Config columns needed for standard pipeline tracking
if (is.null(cfgT$Forward_Primer)) {
  cfgT$Forward_Primer <- ""
}
if (is.null(cfgT$Reverse_Primer)) {
  cfgT$Reverse_Primer <- ""
}

test_that("amplicanFilter correctly removes EOP and PD from real data", {
  # We test amplicanFilter on the original test_aln
  # Since it is a well-formed real data object without PRIMER DIMERS, we first
  # inject a manual EOP event and a manual PD event to ensure they get removed.
  
  aln_test <- data.table::copy(aln)
  
  # Inject EOP (Event Overlapping Primer) - deletion at the very beginning of Amplicon
  eop_event <- aln_test[1, ]
  eop_event$start <- 1 # Starts before/on Forward_Primer
  eop_event$end <- 10 
  eop_event$type <- "deletion"
  eop_event$read_id <- 9999
  
  # Inject PD (Primer Dimer) - large deletion spanning almost the whole amplicon
  pd_event <- aln_test[1, ]
  pd_event$start <- 10
  pd_event$end <- 180
  pd_event$width <- 171
  pd_event$type <- "deletion"
  pd_event$read_id <- 8888
  
  aln_test <- rbind(aln_test, eop_event, pd_event)
  
  filtered_aln <- amplicanFilter(aln_test, cfgT, PRIMER_DIMER = 30)
  
  expect_false(9999 %in% filtered_aln$read_id) # EOP read is gone
  expect_false(8888 %in% filtered_aln$read_id) # PD read is gone
  
  # Ensure genuine valid reads survived
  expect_true(1 %in% filtered_aln$read_id)
})

test_that("amplicanNormalize successfully scrubs events found in controls", {
  aln_test <- data.table::copy(aln)
  cfgT_test <- data.table::copy(cfgT)
  
  # Ensure we have a valid normalization control setup
  # ID_1 is treatment, ID_2 is control
  cfgT_test$Control <- c(FALSE, TRUE)
  cfgT_test$Reads_Filtered <- c(100, 100) # Ensure frequency calculations look sound
  cfgT_test$Group <- c("A", "A")
  cfgT_test$guideRNA <- c("ACTG", "ACTG")
  
  # Add shared noise (mutation present in both control and treatment)
  noise_treatment <- aln_test[1, ]
  noise_treatment$seqnames <- "ID_1"
  noise_treatment$start <- 50
  noise_treatment$end <- 50
  noise_treatment$type <- "mismatch"
  noise_treatment$counts <- 5 # 5/100 = 0.05 freq (>0.01 threshold)
  noise_treatment$read_id <- 7777
  
  noise_control <- aln_test[1, ]
  noise_control$seqnames <- "ID_2"
  noise_control$start <- 50
  noise_control$end <- 50
  noise_control$type <- "mismatch"
  noise_control$counts <- 5 
  noise_control$read_id <- 6666
  
  # Add unique signal (mutation present ONLY in treatment)
  signal_treatment <- aln_test[1, ]
  signal_treatment$seqnames <- "ID_1"
  signal_treatment$start <- 75
  signal_treatment$end <- 75
  signal_treatment$type <- "mismatch"
  signal_treatment$counts <- 20
  signal_treatment$read_id <- 5555
  
  aln_test <- rbind(aln_test, noise_treatment, noise_control, signal_treatment)
  
  norm_aln <- amplicanNormalize(aln_test, cfgT_test, min_freq = 0.01)
  
  # The control noise read stays in the control dataset
  expect_true(6666 %in% norm_aln$read_id)
  
  # The shared noise gets SCRUBBED from the treatment dataset
  expect_false(7777 %in% norm_aln$read_id)
  
  # The unique treatment signal SURVIVES
  expect_true(5555 %in% norm_aln$read_id)
})

test_that("amplicanSummarize counts weighted reads and net frameshifts exactly", {
  events <- data.frame(
    seqnames = c(rep("A", 10), "B", "B"),
    read_id = c(1, 1, 2, 3, 3, 4, 4, 5, 6, 6, 1, 2),
    type = c("deletion", "deletion", "insertion", "deletion", "insertion",
             "mismatch", "mismatch", "mismatch", "deletion", "insertion",
             "insertion", "deletion"),
    width = c(1, 2, 4, 2, 2, 1, 1, 1, 1, 2, 3, 2),
    counts = c(2, 2, 3, 5, 5, 7, 7, 11, 13, 13, 17, 19),
    readType = c(rep(FALSE, 5), TRUE, rep(FALSE, 5), TRUE),
    stringsAsFactors = FALSE)
  # A: reads 1 and 3 are in-frame; read 6 has net +1; read 4 is HDR.
  # B reuses read IDs from A, and "empty" has no events.
  expected <- data.frame(
    ID = c("B", "A", "empty"), HDR = c(19, 7, 0),
    Reads_Del = c(19, 20, 0), Reads_In = c(17, 21, 0),
    Reads_Edited = c(36, 30, 0), Reads_Frameshifted = c(19, 16, 0),
    stringsAsFactors = FALSE)
  summarize <- function(x) {
    amplicanSummarize(data.table::copy(x), data.frame(ID = expected$ID,
                                                   stringsAsFactors = FALSE))
  }

  expect_equal(summarize(events), expected)
  expect_equal(summarize(events[rev(seq_len(nrow(events))), ]), expected)
})
