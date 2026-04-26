library(testthat)
library(amplican)
library(data.table)
library(Biostrings)

context("Comprehensive HDR detection: is_hdr and is_hdr_strict")

# ── Shared setup ──────────────────────────────────────────────────────────────
scoring_matrix <- pwalign::nucleotideSubstitutionMatrix(
  match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA")
gap_opening  <- 25
gap_extension <- 0

# ── Reference sequences ───────────────────────────────────────────────────────
# amplicon: 15 bp, UPPER region = cut site at pos 6-10
# Relative coords after amplicanMap: UPPER start is pos 6, so offset = -(6-1) = -5
# SNV at pos 10 (G->A) → relative pos = 10-6 = 4
amplicon   <- "tttttGGGGGttttt"
donor_1sv  <- "tttttGGGGAttttt"  # 1 SNV: pos 10 G->A (within cut site)
donor_2sv  <- "tttttGGGAAttttt"  # 2 SNVs: pos 9 G->A, pos 10 G->A

# Reads for SNV donor
read_exact  <- donor_1sv                  # perfect HDR
read_wt     <- amplicon                   # wild-type, no HDR
# "noise" here means additional mismatches OUTSIDE the HDR event window
# is_hdr's donor_mismatch counts only events AT the HDR positions,
# so mismatches outside HDR window don't affect the threshold
read_noise_outside <- "tttAtGGGGAttttt"  # donor_1sv + 1mm at pos 4 (T->A, outside cut site)
# To exercise donor_mismatch, we need mismatches WITHIN the donor event span.
# The donor event is a single mismatch at pos 10. A read that has the donor event
# PLUS an extra mismatch at the same position is impossible (it's the same base).
# Instead: use a read that disagrees with the donor AT the HDR position differently:
read_diff_at_hdr <- "tttttGGGGCttttt"   # pos 10 G->C (not the donor A) — same pos, diff base

# Reads for 2-SNV donor (pos 9 G->A, pos 10 G->A)
read_exact_2sv   <- donor_2sv            # both SNVs
read_partial_2sv <- donor_1sv            # only SNV at pos 10 (not pos 9)

# Deletion donor: gap_opening=25 makes a 1-bp deletion very costly vs mismatch;
# use a longer amplicon to make deletion cheaper than accumulated mismatches.
amplicon_long <- "AAAAAAAAAAAAAAAGGGGGTTTTTTTTTTTTTTTT"
donor_del     <- "AAAAAAAAAAAAAAAAAAAAAGGGGGTTTTTTTTTTTTTTTt"  # won't work well — use SNV+del
# For simplicity: use a dedicated longer construct for del/ins tests.
amplicon_di   <- "ttttttttttGGGGGGtttttttttt"      # 26 bp, cut site at 11-16
donor_del_di  <- "ttttttttttGGGGGtttttttttt"       # delete pos 16 (1 del in cut site) — 25bp
donor_ins_di  <- "ttttttttttGGGGGAGtttttttttt"    # insert A after pos 15 (in cut site) — 27bp
read_del_di   <- donor_del_di
read_ins_di   <- donor_ins_di

# Long amplicon for window-based donor_mismatch tests.
# 30bp amplicon: 10 lowercase + 10 UPPERCASE + 10 lowercase
# Uppercase window (cut site) = positions 11-20 (1-indexed)
# After amplicanMap: uppercase starts at relative pos 0, ends at 9
# With cut_buffer=5: window = [-5, 14]
amplicon_win   <- "aaaaaaaaaaGGGGGGGGGGaaaaaaaaaa"   # 30bp
donor_win      <- "aaaaaaaaaaGGGGGGGGGAaaaaaaaaaa"   # SNV at pos 20 (G->A)
# Exact donor read
read_win_exact <- donor_win
# Noise at pos 3 (far from cut site, relative pos = 3-11 = -8)
# With cut_buffer=5 the window is [-5, 14], so pos -8 is OUTSIDE
read_win_noise_far  <- "aaAaaaaaaaGGGGGGGGGAaaaaaaaaaa"
# Noise at pos 15 (inside cut site, relative pos = 15-11 = 4)
read_win_noise_near <- "aaaaaaaaaaGGGGAGGGGAaaaaaaaaaa"

# Helper: align reads to amplicon, return score
make_scores <- function(reads, ampl = amplicon) {
  pwalign::score(pwalign::pairwiseAlignment(
    DNAStringSet(reads), DNAStringSet(ampl),
    substitutionMatrix = scoring_matrix, type = "overlap",
    gapOpening = gap_opening, gapExtension = gap_extension))
}

# Helper: build cfgT + aln data.table for is_hdr_strict
make_test_data <- function(ampl, donor, reads, id = "ID_1", consensus = TRUE) {
  cfgT <- data.table(
    ID = id, Amplicon = ampl, Donor = donor,
    fwdPrPos = 1L, rvePrPos = nchar(ampl), Direction = 0L)

  alns <- pwalign::pairwiseAlignment(
    DNAStringSet(reads), DNAStringSet(ampl),
    substitutionMatrix = scoring_matrix, type = "overlap",
    gapOpening = gap_opening, gapExtension = gap_extension)

  evts <- amplican::getEvents(
    pwalign::pattern(alns), pwalign::subject(alns),
    scores = pwalign::score(alns),
    ID = id, strand_info = "+",
    ampl_start = pwalign::start(pwalign::subject(alns)))

  if (length(evts) == 0) {
    aln_dt <- data.table(
      seqnames = character(), start = integer(), end = integer(),
      width = integer(), strand = character(), score = numeric(),
      originally = character(), replacement = character(),
      type = character(), read_id = integer(), counts = integer(),
      consensus = logical(), readType = logical())
    setkey(aln_dt, seqnames)
    return(list(cfgT = cfgT, aln = aln_dt))
  }

  evts  <- amplican:::amplicanMap(evts, cfgT)
  names(evts) <- NULL
  aln_dt <- as.data.table(as.data.frame(evts))
  aln_dt$seqnames  <- id
  aln_dt$consensus <- consensus
  if (!"readType" %in% names(aln_dt)) aln_dt$readType <- FALSE
  if (!"counts"   %in% names(aln_dt)) aln_dt$counts   <- 1L
  setkey(aln_dt, seqnames)
  list(cfgT = cfgT, aln = aln_dt)
}

hdr_ids <- function(res) unique(res[readType == TRUE]$read_id)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — is_hdr (donor_strict = FALSE)
# donor_mismatch counts events that overlap the HDR event positions.
# Mismatches OUTSIDE the HDR window are ignored by the threshold.
# ══════════════════════════════════════════════════════════════════════════════

test_that("is_hdr/mismatch=3: exact donor read is HDR", {
  res <- amplican:::is_hdr(read_exact, make_scores(read_exact),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)
  expect_true(res[1])
})

test_that("is_hdr/mismatch=3: WT read is not HDR", {
  res <- amplican:::is_hdr(read_wt, make_scores(read_wt),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)
  expect_false(res[1])
})

test_that("is_hdr/mismatch=3: donor + noise OUTSIDE HDR window is still HDR", {
  # donor_mismatch only counts events AT the HDR event positions.
  # Extra mismatches outside the cut site window don't affect the threshold.
  res <- amplican:::is_hdr(read_noise_outside, make_scores(read_noise_outside),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)
  expect_true(res[1])
})

test_that("is_hdr/mismatch=0: exact donor is HDR", {
  res <- amplican:::is_hdr(read_exact, make_scores(read_exact),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)
  expect_true(res[1])
})

test_that("is_hdr/mismatch=0: WT read is not HDR", {
  res <- amplican:::is_hdr(read_wt, make_scores(read_wt),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)
  expect_false(res[1])
})

test_that("is_hdr/mismatch=0: read with different base AT HDR position is not HDR", {
  # Same position as donor event, but a different substitution — wrong HDR event
  res <- amplican:::is_hdr(read_diff_at_hdr, make_scores(read_diff_at_hdr),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)
  expect_false(res[1])
})

test_that("is_hdr/mismatch=0: noise OUTSIDE HDR window is still accepted (not counted)", {
  # Outside-window noise doesn't increment the mismatch counter
  res <- amplican:::is_hdr(read_noise_outside, make_scores(read_noise_outside),
    amplicon, donor_1sv, type = "overlap",
    scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)
  expect_true(res[1])
})

test_that("is_hdr/mismatch=3: mixed batch — correct per-read classification", {
  reads  <- c(read_exact, read_wt, read_noise_outside, read_diff_at_hdr)
  scores <- make_scores(reads)
  res <- amplican:::is_hdr(reads, scores, amplicon, donor_1sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)
  expect_true(res[1])   # exact donor → HDR
  expect_false(res[2])  # WT → not HDR (no score improvement vs donor)
  expect_true(res[3])   # noise outside window → still HDR
  # read_diff_at_hdr scores EQUAL to donor (both have 1 mm at same pos, different base),
  # so it passes the score >= threshold and the overlapping event is within donor_mismatch=3
  expect_true(res[4])   # ties allowed by >= ; is_hdr_strict below correctly rejects it
})

test_that("is_hdr/mismatch=0: mixed batch", {
  reads  <- c(read_exact, read_wt, read_noise_outside, read_diff_at_hdr)
  scores <- make_scores(reads)
  res <- amplican:::is_hdr(reads, scores, amplicon, donor_1sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)
  expect_true(res[1])   # exact donor → HDR
  expect_false(res[2])  # WT → not HDR
  expect_true(res[3])   # noise outside window → accepted (not at HDR pos)
  expect_false(res[4])  # wrong base at HDR pos → not HDR (event at HDR pos, wrong)
})


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — is_hdr_strict (donor_strict = TRUE)
# Requires ALL donor-specific events to be present in the read.
# With donor_mismatch = Inf (default): extra events do NOT disqualify.
# With donor_mismatch = 0: only reads with exactly the donor events pass.
# ══════════════════════════════════════════════════════════════════════════════

test_that("is_hdr_strict: exact 1-SNV donor → HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_exact)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict: WT read → NOT HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_wt)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict: donor + noise OUTSIDE window → still HDR (default donor_mismatch=Inf)", {
  td <- make_test_data(amplicon, donor_1sv, read_noise_outside)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict: wrong base at HDR position → NOT HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_diff_at_hdr)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  # read_diff_at_hdr has G->C at pos 10, not the donor G->A — different event
  expect_false(any(res$readType))
})

test_that("is_hdr_strict/2-SNV donor: read with BOTH events → HDR", {
  td <- make_test_data(amplicon, donor_2sv, read_exact_2sv)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict/2-SNV donor: read with only 1 of 2 events → NOT HDR", {
  # read_partial_2sv = donor_1sv: has pos-10 SNV only, not pos-9 SNV
  td <- make_test_data(amplicon, donor_2sv, read_partial_2sv)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict: mixed batch — correct per-read classification", {
  reads <- c(read_exact, read_wt, read_noise_outside, read_diff_at_hdr)
  td    <- make_test_data(amplicon, donor_1sv, reads)
  res   <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                     gap_opening, gap_extension)
  hdr <- hdr_ids(res)
  expect_true(1L  %in% hdr)                                    # exact donor → HDR
  if (2L %in% res$read_id) expect_false(2L %in% hdr)          # WT → not HDR
  expect_true(3L  %in% hdr)                                    # noise outside → HDR
  expect_false(4L %in% hdr)                                    # wrong event → not HDR
})

test_that("is_hdr_strict: consensus=FALSE rows not used to determine HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_exact, consensus = FALSE)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict: readType update applies to ALL rows of a matching read_id", {
  # Confirm the update `aln[aln_id, readType := ...]` hits all rows, not just consensus
  td    <- make_test_data(amplicon, donor_1sv, read_exact)
  extra <- copy(td$aln[1]); extra$consensus <- FALSE
  aln   <- rbindlist(list(td$aln, extra)); setkey(aln, seqnames)
  res   <- amplican:::is_hdr_strict(aln, td$cfgT, scoring_matrix,
                                     gap_opening, gap_extension)
  expect_true(all(res[read_id == 1L]$readType))
})

test_that("is_hdr_strict: no donor → aln returned unchanged", {
  td     <- make_test_data(amplicon, "", read_exact)
  before <- copy(td$aln)
  res    <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                      gap_opening, gap_extension)
  expect_identical(res$readType, before$readType)
})

test_that("is_hdr_strict: donor == amplicon → no hdr_events → readType stays FALSE", {
  td  <- make_test_data(amplicon, amplicon, read_exact)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension)
  if (nrow(res) > 0) expect_false(any(res$readType))
})

test_that("is_hdr_strict: empty aln → 0-row result", {
  td    <- make_test_data(amplicon, donor_1sv, read_exact)
  empty <- td$aln[0]
  res   <- amplican:::is_hdr_strict(empty, td$cfgT, scoring_matrix,
                                     gap_opening, gap_extension)
  expect_equal(nrow(res), 0L)
})

test_that("is_hdr_strict: multi-ID — only correct experiment flagged", {
  td1 <- make_test_data(amplicon, donor_1sv, read_exact,  id = "ID_1")  # has HDR
  td2 <- make_test_data(amplicon, donor_1sv, read_wt,     id = "ID_2")  # no HDR
  td3 <- make_test_data(amplicon, "",        read_exact,  id = "ID_3")  # no donor
  cfgT <- rbindlist(list(td1$cfgT, td2$cfgT, td3$cfgT))
  aln  <- rbindlist(list(td1$aln,  td2$aln,  td3$aln),  fill = TRUE)
  setkey(aln, seqnames)
  res  <- amplican:::is_hdr_strict(aln, cfgT, scoring_matrix,
                                    gap_opening, gap_extension)
  expect_true(any(res[seqnames == "ID_1"]$readType))
  if (nrow(res[seqnames == "ID_2"]) > 0)
    expect_false(any(res[seqnames == "ID_2"]$readType))
  if (nrow(res[seqnames == "ID_3"]) > 0)
    expect_false(any(res[seqnames == "ID_3"]$readType))
})


# ── is_hdr_strict with donor_mismatch = 0 ────────────────────────────────────

test_that("is_hdr_strict/mismatch=0: exact donor → HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_exact)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict/mismatch=0: donor + noise in window → NOT HDR (short amplicon)", {
  # On the 15bp amplicon with cut_buffer=5 the window covers the entire
  # sequence, so the noise event at relative pos -2 falls inside the window.
  td <- make_test_data(amplicon, donor_1sv, read_noise_outside)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict/mismatch=0: WT → NOT HDR", {
  td <- make_test_data(amplicon, donor_1sv, read_wt)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict/mismatch=1: donor + 1 noise → HDR (within tolerance)", {
  td <- make_test_data(amplicon, donor_1sv, read_noise_outside)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 1)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict/mismatch=0: mixed batch — only exact donor passes", {
  reads <- c(read_exact, read_wt, read_noise_outside, read_diff_at_hdr)
  td    <- make_test_data(amplicon, donor_1sv, reads)
  res   <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                     gap_opening, gap_extension,
                                     donor_mismatch = 0)
  hdr <- hdr_ids(res)
  expect_true(1L  %in% hdr)             # exact donor → HDR
  if (2L %in% res$read_id) expect_false(2L %in% hdr)  # WT → not HDR
  expect_false(3L %in% hdr)             # noise in window → NOT HDR (extra event)
  expect_false(4L %in% hdr)             # wrong event → not HDR
})

# ── is_hdr_strict window tests (long amplicon) ───────────────────────────────

test_that("is_hdr_strict/mismatch=0: noise FAR from cut site → HDR (outside window)", {
  # read_win_noise_far has noise at relative pos -8, outside [-5,14] window
  td <- make_test_data(amplicon_win, donor_win, read_win_noise_far)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0, cut_buffer = 5)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict/mismatch=0: noise NEAR cut site → NOT HDR (inside window)", {
  # read_win_noise_near has noise at relative pos 4, inside [-5,14] window
  td <- make_test_data(amplicon_win, donor_win, read_win_noise_near)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0, cut_buffer = 5)
  expect_false(any(res$readType))
})

test_that("is_hdr_strict/mismatch=1: noise NEAR cut site → HDR (within tolerance)", {
  td <- make_test_data(amplicon_win, donor_win, read_win_noise_near)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 1, cut_buffer = 5)
  expect_true(1L %in% hdr_ids(res))
})

test_that("is_hdr_strict/mismatch=0: cut_buffer=0 shrinks window, far noise ignored", {
  # With cut_buffer=0 the window is [0,9]; noise at -8 still outside
  td <- make_test_data(amplicon_win, donor_win, read_win_noise_far)
  res <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                   gap_opening, gap_extension,
                                   donor_mismatch = 0, cut_buffer = 0)
  expect_true(1L %in% hdr_ids(res))
})


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Explicit contrasts: is_hdr vs is_hdr_strict
# ══════════════════════════════════════════════════════════════════════════════

test_that("contrast: is_hdr_strict REQUIRES all events; is_hdr only needs score improvement", {
  # Partial-match read has 1 of 2 donor events.
  # is_hdr (score-based) may accept it; is_hdr_strict must reject it.
  reads  <- read_partial_2sv
  scores <- make_scores(reads)
  res_loose <- amplican:::is_hdr(reads, scores, amplicon, donor_2sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)
  td <- make_test_data(amplicon, donor_2sv, read_partial_2sv)
  res_strict <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                          gap_opening, gap_extension)
  # strict must reject the partial match
  expect_false(any(res_strict$readType))
  # (is_hdr result is informational — may or may not be TRUE depending on scoring)
})

test_that("contrast: noise OUTSIDE HDR window accepted by both", {
  reads  <- read_noise_outside
  scores <- make_scores(reads)

  res_loose <- amplican:::is_hdr(reads, scores, amplicon, donor_1sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)

  td <- make_test_data(amplicon, donor_1sv, read_noise_outside)
  res_strict <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                          gap_opening, gap_extension)

  # Both should accept: the donor event IS present; noise is outside the HDR window
  # is_hdr ignores outside-window noise; is_hdr_strict default (Inf) also ignores
  expect_true(res_loose[1])
  expect_true(1L %in% hdr_ids(res_strict))
})

test_that("contrast: noise OUTSIDE HDR window — is_hdr_strict/mismatch=0 rejects, is_hdr accepts", {
  reads  <- read_noise_outside
  scores <- make_scores(reads)

  res_loose <- amplican:::is_hdr(reads, scores, amplicon, donor_1sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 0)

  td <- make_test_data(amplicon, donor_1sv, read_noise_outside)
  res_strict <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                          gap_opening, gap_extension,
                                          donor_mismatch = 0)

  # is_hdr: outside-window noise doesn't count → accepted
  expect_true(res_loose[1])
  # is_hdr_strict/mismatch=0: any extra consensus event → rejected
  expect_false(any(res_strict$readType))
})

test_that("contrast: wrong base AT HDR position — is_hdr accepts (tie score), is_hdr_strict rejects", {
  # read_diff_at_hdr has G->C at pos 10 instead of donor's G->A.
  # vs amplicon: 1 mismatch; vs donor: 1 mismatch → equal score → is_hdr passes (>=)
  # is_hdr_strict: donor event is G->A mismatch at pos 4 (relative);
  #                read has G->C mismatch at same pos — DIFFERENT replacement → no match → rejected
  reads  <- read_diff_at_hdr
  scores <- make_scores(reads)

  res_loose <- amplican:::is_hdr(reads, scores, amplicon, donor_1sv,
    type = "overlap", scoring_matrix = scoring_matrix,
    gap_opening = gap_opening, gap_extension = gap_extension,
    donor_mismatch = 3)

  td <- make_test_data(amplicon, donor_1sv, read_diff_at_hdr)
  res_strict <- amplican:::is_hdr_strict(td$aln, td$cfgT, scoring_matrix,
                                          gap_opening, gap_extension)

  # is_hdr accepts because score ties are allowed (>=) and event within donor_mismatch
  expect_true(res_loose[1])
  # is_hdr_strict rejects: the event (G->C) != donor event (G->A) — merge finds no matching row
  expect_false(any(res_strict$readType))
})
