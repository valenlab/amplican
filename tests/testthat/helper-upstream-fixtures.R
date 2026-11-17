# Truth is specified in reference coordinates, independently of event extraction.
up_reference <- function() {
  paste0("acgtcagtacgatcgtacga", "ttcgagcatgacctgatcga",
         "GCTAGTCAGTACCGATGCAT", "cgatgacctagttcgacgta",
         "gatccgtagctagcatcgac")
}

up_change <- function(sequence, position, base) {
  substr(sequence, position, position) <- base
  sequence
}

up_fastq <- function(path, reads, qualities = vapply(reads, function(x) strrep("I", nchar(x)), "")) {
  lines <- unlist(lapply(seq_along(reads), function(i) {
    c(paste0("@read", i), reads[i], "+", qualities[i])
  }), use.names = FALSE)
  writeLines(as.character(lines), path)
}

up_rows <- function(x, columns = names(x)) {
  x <- as.data.frame(x)[, columns, drop = FALSE]
  x[] <- lapply(x, function(column) if (is.factor(column)) as.character(column) else column)
  if (nrow(x)) x <- x[do.call(order, x), , drop = FALSE]
  rownames(x) <- NULL
  x
}

up_fixture <- function(folder) {
  reference <- toupper(up_reference())
  snp <- up_change(reference, 48L, "C")
  insertion <- paste0(substr(reference, 1, 48), "CC", substr(reference, 49, 100))
  deletion <- paste0(substr(reference, 1, 45), substr(reference, 49, 100))
  donor <- up_change(up_change(reference, 55L, "T"), 58L, "A")
  dimer <- paste0(substr(reference, 1, 20), substr(reference, 81, 100))
  cfg <- data.frame(ID = c("T", "R", "C", "E"),
    Barcode = c("treatment", "reverse", "control", "treatment"),
    Forward_Reads = c("T_1.fastq", "R_1.fastq", "C_1.fastq", "T_1.fastq"),
    Reverse_Reads = c("T_2.fastq", "R_2.fastq", "C_2.fastq", "T_2.fastq"),
    Group = c("one", "reverse", "one", "empty"),
    Control = c(FALSE, FALSE, TRUE, FALSE), guideRNA = "GCTAGTCAGTACCGATGCAT",
    Forward_Primer = c(rep("ACGTCAGTACGATCGTACGA", 3), "TTGGCCAATTCGGTACCTGA"),
    Reverse_Primer = c(rep("GTCGATGCTAGCTACGGATC", 3), "TGCACGGTTAGCCTAACGTG"),
    Direction = c(FALSE, TRUE, FALSE, FALSE),
    Amplicon = c(rep(up_reference(), 3), paste0("ttggccaattcggtacctga",
      substr(up_reference(), 21, 80), "cacgttaggctaaccgtgca")),
    Donor = c(donor, "", "", ""), stringsAsFactors = FALSE)
  reads <- list(T = rep(c(reference, snp, insertion, deletion, donor, dimer), c(2, 3, 5, 7, 4, 2)),
                R = rep(c(reference, snp, insertion, deletion), c(2, 3, 5, 7)),
                C = rep(c(reference, snp), c(9, 1)))
  for (id in names(reads)) {
    up_fastq(file.path(folder, paste0(id, "_1.fastq")), paste0("GGAA", reads[[id]]))
    up_fastq(file.path(folder, paste0(id, "_2.fastq")),
      paste0("TTCC", as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reads[[id]])))))
  }
  path <- file.path(folder, "config.csv")
  write.csv(cfg, path, row.names = FALSE)
  list(config = path, cfg = cfg, reads = reads, reference = reference,
       snp = snp, insertion = insertion, deletion = deletion, donor = donor)
}

up_align <- function(fixture, ..., fastqfiles = 0, batch_size = 7, primer_mismatch = 0) {
  suppressMessages(amplicanAlign(fixture$config, dirname(fixture$config),
    fastqfiles = fastqfiles, batch_size = batch_size, primer_mismatch = primer_mismatch, ...))
}

up_pipeline <- function(fixture, results, ..., fastqfiles = 0, batch_size = 7,
                        event_filter = FALSE, pipeline = amplicanPipeline) {
  dir.create(results, showWarnings = FALSE)
  suppressMessages(pipeline(fixture$config, dirname(fixture$config), results,
    knit_reports = FALSE, fastqfiles = fastqfiles, primer_mismatch = 0,
    batch_size = batch_size, event_filter = event_filter, ...))
}

up_summary <- function() {
  data.frame(ID = c("T", "R", "C", "E"), Reads = c(23, 17, 10, 0),
    PRIMER_DIMER = c(2, 0, 0, 0), Low_Score = 0,
    Reads_Filtered = c(21, 17, 10, 0), Reads_Del = c(7, 7, 0, 0),
    Reads_In = c(5, 5, 0, 0), Reads_Edited = c(16, 12, 0, 0),
    Reads_Frameshifted = c(5, 5, 0, 0), HDR = c(4, 0, 0, 0))
}

up_events <- function(conservative = FALSE) {
  expected <- data.frame(seqnames = c("C", rep("T", 4), rep("R", 3)),
    start = c(7, 5, 8, 14, 17, 12, 11, 12),
    width = c(1, 3, 2, 1, 1, 3, 2, 1),
    type = c("mismatch", "deletion", "insertion", "mismatch", "mismatch",
             "deletion", "insertion", "mismatch"),
    originally = c("A", "", "", "A", "C", "", "", "T"),
    replacement = c("C", "", "CC", "T", "A", "", "GG", "G"),
    counts = c(1, 7, 5, 4, 4, 7, 5, 3),
    readType = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
  if (conservative) expected <- rbind(expected, transform(expected[1, ], seqnames = "T", counts = 3))
  expected
}

up_outputs <- function(results) {
  list(summary = data.table::fread(file.path(results, "config_summary.csv")),
       raw = data.table::fread(file.path(results, "alignments", "raw_events.csv")),
       events = data.table::fread(file.path(results, "alignments", "events_filtered_shifted_normalized.csv")),
       barcode = data.table::fread(file.path(results, "barcode_reads_filters.csv")))
}
