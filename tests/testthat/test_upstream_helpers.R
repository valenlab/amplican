library(amplican)
library(testthat)

test_that("primer matches honor overlap and edit limits and prefer the best alignment", {
  primer <- "AGGTCTACGTGGACCCTGCA"
  reads <- c(primer, paste0("TTAA", primer), substr(primer, 11, 20),
             substr(primer, 12, 20), up_change(primer, 6, "A"),
             up_change(up_change(primer, 6, "A"), 12, "A"),
             paste0(substr(primer, 13, 20), "TTTTTTTT", primer), strrep("A", 30))
  expect_equal(amplican:::locate_pr_start(reads, primer, m = 1), c(1, 5, 1, NA, 1, NA, 17, NA))
  expect_equal(amplican:::locate_pr_start(reads[3:4], primer, m = 0, min_overlap = 9), c(1, 1))
  expect_equal(amplican:::locate_pr_start(reads[5:6], primer, m = 0), c(NA_real_, NA_real_))
  expect_equal(amplican:::locate_pr_start(reads[5:6], primer, m = 2), c(1, 1))
  expect_equal(amplican:::locate_pr_start(primer, primer, min_overlap = 20), 1)
  expect_equal(amplican:::locate_pr_start(up_change(primer, 6, "N"), primer, m = 1), 1)
  expect_equal(amplican:::locate_pr_start(character(), primer), numeric())
  expect_error(amplican:::locate_pr_start(primer, primer, min_overlap = 21), "min_overlap")
})

test_that("oriented sequences preserve uppercase groups and swap the correct primers", {
  cfg <- data.frame(ID = c("forward", "reverse", "absent"), Direction = c(0, 1, 0),
    Amplicon = c("aaCGttAAA", "acTGaaCCC", "acgt"), Donor = c("acGT", "acGT", NA),
    Forward_Primer = c("ACG", "TCA", "AC"), Reverse_Primer = c("GAT", "CCA", "GT"))
  expect_identical(get_seq(cfg, "forward"), "aaCGttAAA")
  expect_identical(get_seq(cfg, "reverse"), "GGGttCAgt")
  expect_identical(get_seq(cfg, c("reverse", "forward")), "GGGttCAgt")
  expect_identical(get_seq(cfg, "forward", row = 2), "GGGttCAgt")
  expect_identical(get_seq(cfg, "reverse", "Donor"), "ACgt")
  expect_identical(get_seq(cfg, "absent", "Donor"), "")
  expect_identical(get_seq(cfg, "reverse", "missing"), "")
  expect_identical(amplican:::get_left_primer(cfg, "forward"), "ACG")
  expect_identical(amplican:::get_left_primer(cfg, "reverse"), "CCA")
  expect_identical(amplican:::get_left_primer(cfg, c("reverse", "forward")), "CCA")
  expect_identical(amplican:::get_left_primer(cfg, "forward", row = 2), "CCA")
  expect_identical(amplican:::get_right_primer(cfg, "forward"), "ATC")
  expect_identical(amplican:::get_right_primer(cfg, "reverse"), "TGA")
  expect_identical(amplican:::get_right_primer(cfg, c("reverse", "forward")), "TGA")
  expect_identical(amplican:::get_right_primer(cfg, "forward", row = 2), "TGA")
  expect_identical(amplican:::upperGroups("ACgtTaaGG"), IRanges::IRanges(c(1L, 5L, 8L), c(2L, 5L, 9L)))
  expect_identical(amplican:::upperGroups("acgt"), IRanges::IRanges())
})

test_that("range flipping preserves anchors, bases and metadata for mixed experiments", {
  cfg <- data.frame(ID = c("F", "R", "short"), Direction = c(0, 1, 1),
                    Amplicon = c(strrep("a", 20), strrep("a", 20), strrep("a", 12)))
  events <- data.frame(seqnames = c("R", "F", "R", "short", "R"),
    start = c(3L, 3L, 5L, 2L, 20L), end = c(5L, 5L, 6L, 2L, 20L),
    width = c(3L, 3L, 2L, 1L, 1L), strand = c("+", "-", "-", "+", "*"),
    originally = c("", "", "", "A", "G"), replacement = c("", "", "AC", "C", "T"),
    type = c("deletion", "deletion", "insertion", "mismatch", "mismatch"),
    read_id = c(1, 1, 2, 3, 4), counts = c(2, 3, 5, 7, 11), score = c(10, 20, 30, 40, 50))
  expected <- events
  expected$start <- c(16L, 3L, 16L, 11L, 1L)
  expected$end <- c(18L, 5L, 17L, 11L, 1L)
  expected$strand <- c("-", "-", "+", "-", "*")
  expected$originally <- c("", "", "", "T", "C")
  expected$replacement <- c("", "", "GT", "G", "A")
  flipped <- amplican:::flipRanges(data.table::copy(events), cfg)
  expect_identical(flipped, expected)
  expect_identical(amplican:::flipRanges(data.table::copy(flipped), cfg), events)
  expect_identical(amplican:::flipRanges(data.table::copy(events[FALSE, ]), cfg), events[FALSE, ])
  expect_equal(up_rows(amplican:::flipRanges(events[5:1, ], cfg[3:1, ])), up_rows(expected))
})

test_that("mapping uses each experiment's oriented first uppercase position", {
  cfg <- data.frame(ID = c("A", "B", "none"), Direction = c(0, 1, 0),
    Amplicon = c("aaaCCaaGGaa", "aaCCCaaaaaa", "acgt"))
  events <- data.frame(seqnames = c("B", "A", "B", "A"),
    start = c(7L, 3L, 11L, 8L), end = c(7L, 4L, 11L, 8L),
    width = c(1L, 2L, 1L, 1L), strand = "+", type = c("mismatch", "insertion", "deletion", "mismatch"),
    read_id = c(2, 1, 2, 3), counts = c(2, 3, 2, 5), replacement = c("A", "CC", "", "T"))
  mapped <- as.data.frame(amplicanMap(events, cfg))
  expect_identical(mapped$start, c(0L, -1L, 4L, 4L))
  expect_identical(mapped$end, c(0L, 0L, 4L, 4L))
  expect_equal(mapped[, names(events)[!names(events) %in% c("seqnames", "start", "end", "strand")]],
               events[, names(events)[!names(events) %in% c("seqnames", "start", "end", "strand")]])
  expect_equal(up_rows(as.data.frame(amplicanMap(events[4:1, ], cfg[3:1, ]))), up_rows(mapped))
  no_upper <- events[1, ]; no_upper$seqnames <- "none"
  expect_warning(absent <- amplicanMap(no_upper, cfg), "without UPPER case")
  expect_length(absent, 0)
  mixed <- rbind(events[1:2, ], no_upper, events[3:4, ])
  expect_warning(mixed_mapped <- amplicanMap(mixed, cfg), "without UPPER case")
  expect_equal(up_rows(as.data.frame(mixed_mapped)), up_rows(mapped))
  expect_length(amplicanMap(events[FALSE, ], cfg), 0)
})

test_that("alignment events retain read identities, scores and clipped ends", {
  reference <- Biostrings::DNAString("ACGTCAGTACGATCGTACGA")
  reads <- Biostrings::DNAStringSet(c("ACGTCAGTACGATCGTACGA", "ACGTCAGTCCGATCGTACGA",
                                    "ACGTCAGTACGA", "ACGATCGTACGA"))
  alignment <- pwalign::pairwiseAlignment(reads, reference, type = "overlap")
  plus <- as.data.frame(getEventInfo(alignment, "A", 4L, "+"))
  minus <- as.data.frame(getEventInfo(alignment, "A", 4L, "-"))
  cols <- c("read_id", "start", "end", "type", "originally", "replacement")
  expect_equal(up_rows(plus, cols), up_rows(data.frame(read_id = c("2", "3"),
    start = c(12L, 16L), end = c(12L, 23L), type = c("mismatch", "deletion"),
    originally = c("A", ""), replacement = c("C", ""))))
  expect_equal(up_rows(minus, cols), up_rows(data.frame(read_id = c("2", "4"),
    start = c(12L, 1L), end = c(12L, 11L), type = c("mismatch", "deletion"),
    originally = c("A", ""), replacement = c("C", ""))))
  expect_equal(plus$score, as.numeric(pwalign::score(alignment))[as.integer(plus$read_id)])
  expect_equal(minus$score, as.numeric(pwalign::score(alignment))[as.integer(minus$read_id)])
  expect_true(all(plus$strand == "+"))
  expect_true(all(minus$strand == "-"))
  expect_length(getEventInfo(NULL, "A", 1), 0)
  expect_error(getEventInfo(alignment, "A", 0), "less than 1")
})

test_that("extraction attaches multiplicities and HDR flags to each mate and experiment", {
  reference <- Biostrings::DNAString("ACGTCAGTACGATCGTACGA")
  alignment <- pwalign::pairwiseAlignment(Biostrings::DNAStringSet(c(
    "ACGTCAGTCCGATCGTACGA", "ACGTCAGTACGATCGTACGA", "ACGTCAGTACGCTCGTACGA")),
    reference, type = "global")
  obj <- methods::new("AlignmentsExperimentSet",
    fwdReads = list(A = alignment, B = alignment, E = NULL),
    rveReads = list(A = alignment, B = alignment, E = NULL),
    fwdReadsType = list(A = c(TRUE, FALSE, FALSE), B = c(FALSE, FALSE, TRUE), E = NULL),
    rveReadsType = list(A = c(FALSE, FALSE, TRUE), B = c(TRUE, FALSE, FALSE), E = NULL),
    readCounts = list(A = c(2, 3, 5), B = c(7, 11, 13), E = NULL),
    experimentData = data.frame(ID = c("A", "B", "E"), Barcode = "one", Direction = c(0, 1, 0),
      Amplicon = rep(as.character(reference), 3), fwdPrPos = c(4L, 1L, 1L)))
  expect_no_warning(events <- extractEvents(obj))
  cols <- c("seqnames", "read_id", "start", "strand", "originally", "replacement", "counts", "readType")
  expected <- data.frame(seqnames = rep(c("A", "B"), each = 4), read_id = rep(c(1L, 3L), 4),
    start = c(12L, 15L, 12L, 15L, 12L, 9L, 12L, 9L),
    strand = c("+", "+", "-", "-", "-", "-", "+", "+"),
    originally = c(rep("A", 4), rep("T", 4)),
    replacement = c(rep("C", 4), rep("G", 4)),
    counts = c(2, 5, 2, 5, 7, 13, 7, 13),
    readType = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))
  events$read_id <- as.integer(events$read_id)
  expect_equal(up_rows(events, cols), up_rows(expected))
  expect_equal(events$score, as.numeric(pwalign::score(alignment))[events$read_id])
  cfg <- experimentData(obj)
  cfg$fwdPrPos <- NULL
  experimentData(obj) <- cfg
  expect_no_warning(no_offset <- extractEvents(obj))
  no_offset$read_id <- as.integer(no_offset$read_id)
  expected$start[expected$seqnames == "A"] <- rep(c(9L, 12L), 2)
  expect_equal(up_rows(no_offset, cols), up_rows(expected))
  cfg$fwdPrPos <- NA_integer_
  experimentData(obj) <- cfg
  expect_no_warning(missing_offset <- extractEvents(obj))
  missing_offset$read_id <- as.integer(missing_offset$read_id)
  expect_equal(up_rows(missing_offset, cols), up_rows(expected))
  lookups <- 0L
  testthat::with_mocked_bindings({
    default <- extractEvents(obj)
    expect_equal(lookups, 0L)
    explicit <- extractEvents(obj, use_parallel = TRUE)
    expect_equal(lookups, 1L)
    expect_equal(up_rows(default), up_rows(explicit))
  }, bpparam = function(...) {
    lookups <<- lookups + 1L
    BiocParallel::SerialParam()
  }, .package = "BiocParallel")
  empty <- methods::new("AlignmentsExperimentSet")
  expect_equal(nrow(extractEvents(empty)), 0L)
})
