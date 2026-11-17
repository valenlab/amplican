library(amplican)
library(testthat)

test_that("mismatch labels stay with their events across insertion columns", {
  actual <- getEvents(
    Biostrings::DNAStringSet(c("ATTAT", "AAGT")),
    Biostrings::DNAStringSet(c("AC-GT", "ACGT")),
    scores = c(11, 22), ID = "fixture")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(3, 2, 2, 3), width = 1),
    strand = "+", originally = c("", "C", "C", "G"),
    replacement = c("T", "A", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = c("1", "2", "1", "1"), score = c(11, 22, 11, 11))
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})

test_that("an insertion cannot reorder mismatch labels within one read", {
  actual <- getEvents(
    Biostrings::DNAStringSet("ATTATAT"),
    Biostrings::DNAStringSet("AC-GTAA"), scores = 17, ID = "fixture")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(3, 6, 2, 3), width = 1),
    strand = "+", originally = c("", "A", "C", "G"),
    replacement = c("T", "T", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = "1", score = 17)
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})

test_that("event content is independent of read order and chunk boundaries", {
  pattern <- Biostrings::DNAStringSet(c("ATTAT", "AAGT", "ACGT"))
  subject <- Biostrings::DNAStringSet(c("AC-GT", "ACGT", "ACGT"))
  scores <- c(11, 22, 33)
  extract <- function(ids) {
    events <- getEvents(pattern[ids], subject[ids], scores[ids], ID = "fixture")
    if (length(events)) {
      events$read_id <- as.character(ids[as.integer(events$read_id)])
      names(events) <- events$read_id
    }
    events
  }
  canonical <- function(events) {
    events[order(as.integer(events$read_id), GenomicRanges::start(events),
                 GenomicRanges::end(events), events$type)]
  }
  whole <- canonical(extract(1:3))

  expect_identical(canonical(extract(3:1)), whole)
  expect_identical(canonical(c(extract(1L), extract(2L), extract(3L))), whole)
  expect_identical(canonical(c(extract(1:2), extract(3L))), whole)
  expect_identical(extract(3L), GenomicRanges::GRanges())
})

test_that("mismatch extraction preserves per-read offsets and strand metadata", {
  actual <- getEvents(
    Biostrings::DNAStringSet(c("ATTAT", "AAGT")),
    Biostrings::DNAStringSet(c("AC-GT", "ACGT")),
    scores = c(11, 22), ID = "fixture", ampl_shift = c(4L, 8L),
    ampl_start = c(2L, 5L), strand_info = "-")
  expected <- GenomicRanges::GRanges(
    seqnames = "fixture", ranges = IRanges::IRanges(c(7, 13, 6, 7), width = 1),
    strand = "-", originally = c("", "C", "C", "G"),
    replacement = c("T", "A", "T", "A"),
    type = c("insertion", "mismatch", "mismatch", "mismatch"),
    read_id = c("1", "2", "1", "1"), score = c(11, 22, 11, 11))
  names(expected) <- expected$read_id

  expect_identical(actual, expected)
})

# Scan columns independently of compareStrings, RleList, and extraction helpers.
reference_events <- function(pattern, subject, read_id, score, offset, strand) {
  query <- strsplit(pattern, "")[[1]]
  ref <- strsplit(subject, "")[[1]]
  events <- list()
  position <- 1L + offset
  column <- 1L
  last_deletion <- 0L
  while (column <= length(ref)) {
    start <- position
    originally <- replacement <- ""
    if (ref[column] == "-") {
      type <- "insertion"
      while (column <= length(ref) && ref[column] == "-") {
        replacement <- paste0(replacement, query[column])
        column <- column + 1L
      }
      end <- start + nchar(replacement) - 1L
    } else {
      position <- position + 1L
      end <- start
      if (query[column] == "-") {
        type <- "deletion"
        # Insertions do not interrupt consecutive deleted reference bases.
        if (last_deletion > 0L && events[[last_deletion]]$end == start - 1L) {
          events[[last_deletion]]$end <- end
          column <- column + 1L
          next
        }
      } else if (query[column] != ref[column]) {
        type <- "mismatch"
        originally <- ref[column]
        replacement <- query[column]
      } else {
        column <- column + 1L
        next
      }
      column <- column + 1L
    }
    events[[length(events) + 1L]] <- data.frame(
      start = start, end = end, originally = originally,
      replacement = replacement, type = type, stringsAsFactors = FALSE)
    if (type == "deletion") last_deletion <- length(events)
  }
  events <- if (length(events)) do.call(rbind, events) else data.frame(
    start = integer(), end = integer(), originally = character(),
    replacement = character(), type = character())
  n <- nrow(events)
  data.frame(
    seqnames = rep("generated", n), start = events$start, end = events$end,
    width = events$end - events$start + 1L, strand = rep(strand, n),
    originally = events$originally, replacement = events$replacement,
    type = events$type, read_id = rep(as.character(read_id), n),
    score = rep(score, n), stringsAsFactors = FALSE)
}

event_cases <- local({
  operations <- c("", unlist(lapply(1:4, function(n) {
    grid <- expand.grid(rep(list(c("M", "X", "I", "D")), n),
                        stringsAsFactors = FALSE)
    apply(grid, 1, paste0, collapse = "")
  })))
  operations <- paste0("M", operations, "M")
  pairs <- lapply(strsplit(operations, ""), function(ops) {
    # Every subject consumes a prefix of the same reference sequence.
    bases <- c("A", "C", "G", "T")
    ref <- bases[(cumsum(ops != "I") - 1L) %% 4L + 1L]
    query <- ref
    query[ops == "X"] <- c(A = "C", C = "G", G = "T", T = "A")[ref[ops == "X"]]
    query[ops == "I"] <- bases[seq_along(ops)[ops == "I"] %% 4L + 1L]
    query[ops == "D"] <- "-"
    ref[ops == "I"] <- "-"
    c(pattern = paste0(query, collapse = ""),
      subject = paste0(ref, collapse = ""))
  })
  data.frame(operations = operations, do.call(rbind, pairs),
             stringsAsFactors = FALSE)
})

canonical_events <- function(events) {
  if (!length(events)) return(GenomicRanges::GRanges())
  events[order(as.integer(events$read_id), GenomicRanges::start(events),
               GenomicRanges::end(events), events$type)]
}

extract_event_cases <- function(ids, offset = 0L, strand = "+") {
  events <- getEvents(
    event_cases$pattern[ids], event_cases$subject[ids], scores = ids * 11,
    ID = "generated", ampl_shift = ids %% 5L + 1L + offset,
    ampl_start = ids %% 7L + 1L, strand_info = strand)
  if (length(events)) {
    events$read_id <- as.character(ids[as.integer(events$read_id)])
    names(events) <- events$read_id
  }
  events
}

test_that("generated alignments match an independent event reference", {
  ids <- seq_len(nrow(event_cases))
  for (strand in c("+", "-")) {
    actual <- canonical_events(extract_event_cases(ids, strand = strand))
    actual <- as.data.frame(actual, row.names = seq_along(actual))
    actual$seqnames <- as.character(actual$seqnames)
    actual$strand <- as.character(actual$strand)
    for (i in ids) {
      expected <- reference_events(
        event_cases$pattern[i], event_cases$subject[i], i, i * 11,
        offset = i %% 5L + i %% 7L, strand = strand)
      read_events <- actual[actual$read_id == as.character(i), ]
      expected <- expected[order(expected$start, expected$end, expected$type), ]
      rownames(read_events) <- rownames(expected) <- NULL
      expect_identical(read_events, expected,
        info = paste("operations", event_cases$operations[i],
                     "pattern", event_cases$pattern[i],
                     "subject", event_cases$subject[i], "strand", strand))
    }
  }
})

test_that("generated events preserve read order, chunk, and shift invariants", {
  ids <- seq_len(nrow(event_cases))
  whole <- canonical_events(extract_event_cases(ids))
  expect_identical(canonical_events(extract_event_cases(rev(ids))), whole)
  chunks <- split(ids, cut(ids, c(0L, 1L, 170L, length(ids)), labels = FALSE))
  expect_identical(
    canonical_events(Reduce(c, lapply(chunks, extract_event_cases))), whole)
  expect_identical(canonical_events(extract_event_cases(ids, offset = 7L)),
                   GenomicRanges::shift(whole, 7L))
})
