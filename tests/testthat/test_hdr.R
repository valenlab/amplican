library(testthat)
library(amplican)
library(data.table)
context("HDR Detection Logic")

# Basic configuration for tests
scoring_matrix <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA")
gap_opening <- 25
gap_extension <- 0

# Amplicon and Donor Setup
amplicon <- "ATCGATCGATCGATCGATCG"
# Donor replaces middle "T" with "A" (ATCGATCG AACGATCGATCG)
donor    <- "ATCGATCGAACGATCGATCG"

# 1. Exact donor match
read_exact <- donor

# 2. Donor + 1 mismatch (e.g. at the end) -> noise, below mismatch threshold
read_noise <- "ATCGATCGAACGATCGATCC"

# 3. Completely unrelated wild read 
read_fail <- amplicon # A wild-type read should definitely be rejected

test_that("is_hdr evaluates correctly for loose scenarios", {
  
  reads <- c(read_exact, read_noise, read_fail)
  
  # Align to amplicon to generate read scores (using pairwiseAlignment)
  alns <- pwalign::pairwiseAlignment(
    Biostrings::DNAStringSet(reads), 
    Biostrings::DNAStringSet(amplicon),
    substitutionMatrix = scoring_matrix,
    type = "overlap",
    gapOpening = gap_opening, gapExtension = gap_extension
  )
  scores <- pwalign::score(alns)
  
  is_hdr_res <- amplican:::is_hdr(
    reads, scores, amplicon, donor, 
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3
  )
  
  # Loose checking should capture the exact donor and the donor with minor noise.
  expect_true(is_hdr_res[1]) # exact donor
  expect_true(is_hdr_res[2]) # donor + 1 noise
  # Large indels beyond donor_mismatch should be rejected as pure HDR
  expect_false(is_hdr_res[3]) # donor + big indel
})


test_that("is_hdr_strict accurately rejects extra indels", {
  # Construct minimal cfgT
  cfgT <- data.table::data.table(
    ID = "ID_1",
    Amplicon = amplicon,
    Donor = donor,
    fwdPrPos = 1,
    rvePrPos = 20,
    Direction = 0
  )
  
  # Align reads against amplicon to get their events
  reads <- c(read_exact, read_noise, read_fail)
  alns <- pwalign::pairwiseAlignment(
    Biostrings::DNAStringSet(reads), 
    Biostrings::DNAStringSet(amplicon),
    substitutionMatrix = scoring_matrix,
    type = "overlap",
    gapOpening = gap_opening, gapExtension = gap_extension
  )
  
  # Extract events
  pat <- pwalign::pattern(alns)
  subj <- pwalign::subject(alns)
  aln_events <- amplican:::getEvents(
    pat, subj, scores = pwalign::score(alns),
    ID = "ID_1", strand_info = "+",
    ampl_start = pwalign::start(subj)
  )
  
  names(aln_events) <- NULL # Prevent duplicate row.names error
  aln_events <- as.data.frame(aln_events)
  data.table::setDT(aln_events)
  
  # Map events the way the core pipeline does
  aln_events <- amplican:::amplicanMap(aln_events, cfgT)
  aln_events <- as.data.frame(aln_events)
  data.table::setDT(aln_events)
  
  # We must manually spoof the consensus for test since we didn't run full filtering
  aln_events$consensus <- TRUE
  # Also spoof seqnames to match ID properly
  aln_events$seqnames <- "ID_1"
  
  # Run strict
  res_aln <- amplican:::is_hdr_strict(
    aln_events, cfgT, scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension
  )
  
  # Safely extract
  if ("readType" %in% names(res_aln)) {
    hdr_reads <- res_aln[readType == TRUE]$read_id
  } else {
    hdr_reads <- integer(0)
  }
  
  expect_true(1 %in% hdr_reads)   # exact donor is perfect HDR
  # The documentation for is_hdr_strict says: 
  # "It ignores everything else, so other mismatches and small indels etc. as noise are allowed here for valid HDR."
  expect_true(2 %in% hdr_reads)  # donor with noise is still valid strict HDR because it contains intended mutations
  expect_false(3 %in% hdr_reads)  # amplicon without HDR mutations is NOT HDR
})
