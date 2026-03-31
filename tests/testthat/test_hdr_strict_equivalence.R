library(testthat)
library(amplican)
library(data.table)
library(Biostrings)

context("is_hdr_strict equivalence tests")

# Shared setup
scoring_matrix <- pwalign::nucleotideSubstitutionMatrix(
  match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA")
gap_opening <- 25
gap_extension <- 0

# Helper: build cfgT + aln from amplicon/donor/reads
# Returns list(cfgT, aln) ready for is_hdr_strict
make_hdr_test_data <- function(amplicon, donor, reads, id = "ID_1",
                               consensus = TRUE) {
  cfgT <- data.table(
    ID = id, Amplicon = amplicon, Donor = donor,
    fwdPrPos = 1L, rvePrPos = nchar(amplicon), Direction = 0L)

  alns <- pwalign::pairwiseAlignment(
    DNAStringSet(reads), DNAStringSet(amplicon),
    substitutionMatrix = scoring_matrix, type = "overlap",
    gapOpening = gap_opening, gapExtension = gap_extension)

  aln_events <- amplican::getEvents(
    pwalign::pattern(alns), pwalign::subject(alns),
    scores = pwalign::score(alns),
    ID = id, strand_info = "+",
    ampl_start = pwalign::start(pwalign::subject(alns)))

  if (length(aln_events) == 0) {
    aln_dt <- data.table(
      seqnames = character(), start = integer(), end = integer(),
      width = integer(), strand = character(), score = numeric(),
      originally = character(), replacement = character(),
      type = character(), read_id = integer(), counts = integer(),
      consensus = logical(), readType = logical())
    setkey(aln_dt, seqnames)
    return(list(cfgT = cfgT, aln = aln_dt))
  }

  aln_events <- amplican:::amplicanMap(aln_events, cfgT)
  names(aln_events) <- NULL
  aln_dt <- as.data.table(as.data.frame(aln_events))
  aln_dt$seqnames <- id
  aln_dt$consensus <- consensus
  if (!"readType" %in% names(aln_dt)) aln_dt$readType <- FALSE
  setkey(aln_dt, seqnames)
  list(cfgT = cfgT, aln = aln_dt)
}

# ---- Test sequences ----
amplicon <- "ATCGATCGATCGATCGATCG"
# Donor: replaces position 9 T->A
donor <- "ATCGATCGAACGATCGATCG"
read_exact <- donor
read_noise <- "ATCGATCGAACGATCGATCC"  # donor + 1 trailing mismatch
read_wt <- amplicon                     # wild-type, no HDR

# Donor with 2 mutations (positions 9 T->A and 15 T->G)
donor_2mut <- "ATCGATCGAACGATCGGTCG"


test_that("test_single_id_exact_donor: exact donor match -> readType TRUE", {
  td <- make_hdr_test_data(amplicon, donor, read_exact)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_true("readType" %in% names(res))
  hdr_reads <- res[readType == TRUE]$read_id
  expect_true(1L %in% hdr_reads)
})


test_that("test_single_id_no_donor: donor='' -> aln unchanged", {
  td <- make_hdr_test_data(amplicon, "", read_exact)
  aln_before <- copy(td$aln)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  # readType should not have been modified
  expect_identical(res$readType, aln_before$readType)
})


test_that("test_single_id_no_match: WT read has no HDR events -> readType FALSE", {
  # WT read is identical to amplicon, produces no events -> empty aln
  # Use a read that differs from amplicon but NOT in the donor mutation position
  read_non_donor <- "GTCGATCGATCGATCGATCG"  # first base changed A->G, not donor mut
  td <- make_hdr_test_data(amplicon, donor, read_non_donor)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_true(nrow(res) > 0)
  expect_false(any(res$readType))
})


test_that("test_single_id_partial_match: only 1 of 2 HDR events -> FALSE", {
  # Use donor with 2 mutations
  # Read only has the first mutation, not the second
  read_partial <- "ATCGATCGAACGATCGATCG"  # has mut at pos 9 but not pos 15
  td <- make_hdr_test_data(amplicon, donor_2mut, read_partial)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  if ("readType" %in% names(res) && nrow(res) > 0) {
    hdr_reads <- res[readType == TRUE]$read_id
    expect_false(1L %in% hdr_reads)
  }
})


test_that("test_multi_id_mixed: 3 IDs, only donor-matching one gets TRUE", {
  # ID_1: has donor, read matches -> TRUE
  td1 <- make_hdr_test_data(amplicon, donor, read_exact, id = "ID_1")
  # ID_2: has donor, read is WT -> FALSE
  td2 <- make_hdr_test_data(amplicon, donor, read_wt, id = "ID_2")
  # ID_3: no donor -> unchanged
  td3 <- make_hdr_test_data(amplicon, "", read_exact, id = "ID_3")

  cfgT <- rbindlist(list(td1$cfgT, td2$cfgT, td3$cfgT))
  aln <- rbindlist(list(td1$aln, td2$aln, td3$aln), fill = TRUE)
  setkey(aln, seqnames)

  res <- amplican:::is_hdr_strict(aln, cfgT, scoring_matrix,
                                   gap_opening, gap_extension)

  # ID_1 read should be HDR
  id1_hdr <- res[seqnames == "ID_1" & readType == TRUE]$read_id
  expect_true(1L %in% id1_hdr)

  # ID_2 read should NOT be HDR
  id2_hdr <- res[seqnames == "ID_2" & readType == TRUE]$read_id
  expect_length(id2_hdr, 0)

  # ID_3 should have readType FALSE (no donor, untouched)
  expect_false(any(res[seqnames == "ID_3"]$readType))
})


test_that("test_consensus_false_ignored: non-consensus events not considered", {
  td <- make_hdr_test_data(amplicon, donor, read_exact, consensus = FALSE)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  # Even though the read matches donor, consensus=FALSE means it shouldn't be HDR
  if (nrow(res) > 0) {
    expect_false(any(res$readType))
  }
})


test_that("test_non_consensus_rows_get_readType: update applies to all rows for read_id", {
  # Build data with consensus=TRUE
  td <- make_hdr_test_data(amplicon, donor, read_exact)
  # Add a duplicate row with consensus=FALSE for the same read_id
  extra_row <- td$aln[1]
  extra_row$consensus <- FALSE
  aln <- rbindlist(list(td$aln, extra_row))
  setkey(aln, seqnames)

  res <- amplican:::is_hdr_strict(aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)

  # ALL rows for this read_id in ID_1 should have readType updated
  read1_types <- res[read_id == 1L]$readType
  expect_true(all(read1_types))
})


test_that("test_multiple_reads_same_id: one HDR read, one not", {
  reads <- c(read_exact, read_wt)
  td <- make_hdr_test_data(amplicon, donor, reads)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)

  hdr_reads <- unique(res[readType == TRUE]$read_id)
  non_hdr_reads <- unique(res[readType == FALSE]$read_id)

  # Read 1 (exact donor) should be HDR
  expect_true(1L %in% hdr_reads)
  # Read 2 (WT) should not
  # WT read may or may not have events; if no events, it won't appear in aln at all
  # But if it does, it should be FALSE
  if (2L %in% res$read_id) {
    expect_true(2L %in% non_hdr_reads)
  }
})


test_that("test_noise_mismatch_still_hdr: donor + noise is still HDR", {
  # read_noise has the donor mutation + an extra mismatch at the end
  td <- make_hdr_test_data(amplicon, donor, read_noise)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  hdr_reads <- res[readType == TRUE]$read_id
  expect_true(1L %in% hdr_reads)
})


test_that("test_empty_aln: empty alignment table returned unchanged", {
  td <- make_hdr_test_data(amplicon, donor, read_exact)
  empty_aln <- td$aln[0]  # 0-row DT with correct schema
  res <- amplican:::is_hdr_strict(empty_aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_equal(nrow(res), 0)
})


test_that("test_no_events_from_alignment: identical donor/amplicon -> no HDR events -> unchanged", {
  # When donor == amplicon, pairwiseAlignment produces no events
  td <- make_hdr_test_data(amplicon, amplicon, read_exact)
  aln_before <- copy(td$aln)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  # readType should remain FALSE (no HDR events to match against)
  if (nrow(res) > 0) {
    expect_false(any(res$readType))
  }
})
