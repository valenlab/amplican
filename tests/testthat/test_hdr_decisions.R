library(amplican)
library(testthat)

hdr_decision_amplicon <- "tttttACGTCAGTACgtcta"
hdr_decision_donor <- hdr_decision_amplicon
substr(hdr_decision_donor, 8, 8) <- "A"
substr(hdr_decision_donor, 12, 12) <- "C"
hdr_decision_matrix <- pwalign::nucleotideSubstitutionMatrix(
  match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA")

hdr_decision_events <- function() {
  read <- c(1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 1, 1)
  start <- c(2, 6, 2, 2, 6, 2, 6, 4, 2, 6, -4, 2, 6, -4, 2, 6)
  data.table::data.table(
    seqnames = c(rep("A", 14), "B", "B"), read_id = read,
    start = start, end = start, width = 1L, strand = "+", score = 100,
    originally = c(rep("G", 7), "C", "G", "G", "T", "G", "G", "T", "G", "G"),
    replacement = c("A", "C", "A", "A", "T", "A", "C", "G", "A", "C", "A", "A", "C", "A", "A", "C"),
    type = "mismatch", counts = 1L, consensus = c(rep(TRUE, 12), FALSE, rep(TRUE, 3)),
    readType = FALSE)
}

test_that("strict HDR requires every donor event and respects noise-window boundaries", {
  cfg <- data.table::data.table(ID = c("A", "B"), Amplicon = hdr_decision_amplicon,
    Donor = c(hdr_decision_donor, ""), Direction = 0L, fwdPrPos = 1L, rvePrPos = 20L)
  events <- hdr_decision_events()
  classify <- function(tolerance = Inf, buffer = 0, x = events, config = cfg) {
    amplican:::is_hdr_strict(data.table::copy(x), config, hdr_decision_matrix,
      donor_mismatch = tolerance, cut_buffer = buffer)
  }
  expected <- function(ids) events$seqnames == "A" & events$read_id %in% ids
  expect_identical(classify()$readType, expected(c(1, 4, 5)))
  expect_identical(classify(0)$readType, expected(c(1, 5)))
  expect_identical(classify(1)$readType, expected(c(1, 4, 5)))
  expect_identical(classify(0, 5)$readType, expected(1))
  expect_identical(classify(0, 4)$readType, expected(1))
  expect_identical(classify(0, 3)$readType, expected(c(1, 5)))
  expect_identical(classify(x = events[16:1])$readType, rev(expected(c(1, 4, 5))))
  expect_identical(classify(x = events[0]), events[0])
  no_donor <- data.table::copy(cfg)
  no_donor$Donor <- ""
  expect_identical(classify(config = no_donor), events)
  identical_donor <- data.table::copy(cfg)
  identical_donor$Donor <- identical_donor$Amplicon
  expect_identical(classify(config = identical_donor), events)
  no_consensus <- data.table::copy(events)
  no_consensus$consensus <- FALSE
  expect_identical(classify(x = no_consensus), no_consensus)
  different <- data.table::copy(cfg)
  reference_b <- tolower(hdr_decision_amplicon)
  substr(reference_b, 8, 8) <- "t"
  substr(reference_b, 12, 12) <- "t"
  substr(reference_b, 10, 14) <- toupper(substr(reference_b, 10, 14))
  different$Amplicon[2] <- reference_b
  events_b <- data.table::copy(events)
  events_b[seqnames == "B", originally := "T"]
  expect_identical(classify(x = events_b, config = different)$readType, expected(c(1, 4, 5)))

  # More experiments than configuration columns catches iteration over the wrong dimension.
  n <- ncol(cfg) + 1L
  many <- cfg[rep(1L, n)]
  many$ID <- paste0("sample", seq_len(n))
  exact <- data.table::rbindlist(lapply(seq_len(n), function(i) {
    x <- data.table::copy(events[1:2])
    x$seqnames <- many$ID[i]
    x
  }))
  expect_identical(classify(x = exact, config = many)$readType, rep(TRUE, nrow(exact)))
})

test_that("permissive HDR tolerates only the declared donor-site differences", {
  partial <- hdr_decision_amplicon
  substr(partial, 8, 8) <- "A"
  wrong <- hdr_decision_donor
  substr(wrong, 12, 12) <- "T"
  outside <- hdr_decision_donor
  substr(outside, 2, 2) <- "A"
  inside <- hdr_decision_donor
  substr(inside, 10, 10) <- "G"
  wrong_both <- hdr_decision_amplicon
  substr(wrong_both, 8, 8) <- "T"
  substr(wrong_both, 12, 12) <- "T"
  wrong_run <- hdr_decision_donor
  substr(wrong_run, 8, 10) <- "GAG"
  reads <- c(hdr_decision_donor, hdr_decision_amplicon, partial, wrong, outside, inside, wrong_both, wrong_run)
  scores <- pwalign::score(pwalign::pairwiseAlignment(
    Biostrings::DNAStringSet(reads), Biostrings::DNAStringSet(hdr_decision_amplicon),
    substitutionMatrix = hdr_decision_matrix, type = "overlap", gapOpening = 25, gapExtension = 0))
  classify <- function(ids = seq_along(reads), tolerance = 0) {
    amplican:::is_hdr(reads[ids], scores[ids], hdr_decision_amplicon, hdr_decision_donor,
      scoring_matrix = hdr_decision_matrix, donor_mismatch = tolerance)
  }
  expect_identical(classify(), c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(classify(tolerance = 1), c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_identical(classify(tolerance = 2), c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_identical(classify(6:1), rev(c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)))
  expect_identical(classify(2), FALSE)
  expect_identical(classify(1), TRUE)
})

test_that("HDR decisions distinguish long indels from tolerated single-base noise", {
  reference <- "TACGATCGTACGATTCGAGTACCTGACATGCTAGCTAGGTCAGTACGATCAGCTAGTCAG"
  donor <- reference
  substr(donor, 30, 30) <- "A"
  substr(donor, 40, 40) <- if (substr(reference, 40, 40) == "A") "C" else "A"
  insert_at <- function(sequence, pos) paste0(substr(sequence, 1, pos-1), "CC", substr(sequence, pos, nchar(sequence)))
  remove_at <- function(sequence, pos, n) paste0(substr(sequence, 1, pos-1), substr(sequence, pos+n, nchar(sequence)))
  reads <- c(donor, reference, insert_at(donor, 10), insert_at(donor, 30), remove_at(donor, 30, 2), remove_at(donor, 30, 1))
  matrix <- pwalign::nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA")
  scores <- pwalign::score(pwalign::pairwiseAlignment(Biostrings::DNAStringSet(reads),
    Biostrings::DNAStringSet(reference), substitutionMatrix = matrix, type = "overlap", gapOpening = 3, gapExtension = 1))
  classify <- function(tolerance) amplican:::is_hdr(reads, scores, reference, donor,
    scoring_matrix = matrix, gap_opening = 3, gap_extension = 1, donor_mismatch = tolerance)
  expect_identical(classify(0), c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(classify(1), c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))

  amplicon <- paste0(tolower(substr(reference, 1, 19)), substr(reference, 20, 45), tolower(substr(reference, 46, nchar(reference))))
  for (type in c("insertion", "deletion")) {
    template <- if (type == "insertion") insert_at(reference, 25) else remove_at(reference, 25, 2)
    cfg <- data.table::data.table(ID = "A", Amplicon = amplicon, Donor = template,
      Direction = 0L, fwdPrPos = 1L, rvePrPos = nchar(reference))
    events <- data.table::data.table(seqnames = "A", read_id = c(1L, 2L), start = 5L,
      end = c(6L, 5L), width = c(2L, 1L), type = c(type, "mismatch"), originally = c("", "A"),
      replacement = c(if (type == "insertion") "CC" else "", "C"),
      consensus = TRUE, readType = FALSE, strand = "+", counts = 1L, score = 100)
    actual <- amplican:::is_hdr_strict(events, cfg, matrix, gap_opening = 3, gap_extension = 1)
    expect_identical(actual$readType, c(TRUE, FALSE), info = type)
  }
})
