library(testthat)
context("amplican main making analysis of example files")

# un-comment line below this comment and below test_that to make vignettes
# after any changes to rmarkdown templates
# changes paths so that system.file points to the actual development folder
# devtools::load_all()

config <- system.file("extdata", "config.csv", package = "amplican")
fastq_folder <- system.file("extdata", package = "amplican")
results_folder <- tempdir()
# results_folder <- system.file("extdata", "results", package = "amplican")
# vignettes_path <- system.file('vignettes', package = 'amplican')
#
# devtools::unload()
library(amplican)

test_that("amplican theme templates are loaded properly", {
  expect_gte(
    nchar(system.file("rmarkdown", "templates",
                      "id_report", package = "amplican")), 1)
})

test_that("amplican runs through example files without any issues", {
  expect_warning(
    amplicanPipeline(config, fastq_folder, results_folder,
                     knit_reports = FALSE,
                     primer_mismatch = 0,
                     fastqfiles = 0,
                     continue = FALSE))
})

test_that("fresh and resumed pipelines preserve known event and read counts", {
  exercise_pipeline <- function() {
    folder <- tempfile("pipeline-truth-")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE))
    cfg <- read.csv(config, stringsAsFactors = FALSE)[1, ]
    cfg$Forward_Reads <- cfg$Reverse_Reads <- "known.fastq"
    cfg$Donor <- ""
    amplicon <- toupper(cfg$Amplicon)
    position <- regexpr("[A-Z]", cfg$Amplicon)[1] + 8L
    deletion <- paste0(substr(amplicon, 1L, position - 1L),
                        substr(amplicon, position + 3L, nchar(amplicon)))
    insertion <- paste0(substr(amplicon, 1L, position), "A",
                         substr(amplicon, position + 1L, nchar(amplicon)))
    mismatch <- amplicon
    substr(mismatch, position, position) <- if (
      substr(amplicon, position, position) == "A") "C" else "A"
    reads <- rep(c(amplicon, deletion, insertion, mismatch), c(2, 3, 5, 7))
    fastq <- unlist(lapply(seq_along(reads), function(i) {
      c(paste0("@read", i), reads[i], "+", strrep("I", nchar(reads[i])))
    }))
    writeLines(fastq, file.path(folder, "known.fastq"))
    write.csv(cfg, file.path(folder, "config.csv"), row.names = FALSE)
    results <- file.path(folder, "results")
    dir.create(results)

    expected_summary <- data.frame(
      Reads = 17, Reads_Filtered = 17, Reads_Del = 3, Reads_In = 5,
      Reads_Edited = 8, Reads_Frameshifted = 5, HDR = 0)
    expected_events <- data.frame(
      type = c("deletion", "insertion", "mismatch"), width = c(3, 1, 1),
      counts = c(3, 5, 7), originally = c("", "", substr(amplicon, position, position)),
      replacement = c("", "A", substr(mismatch, position, position)),
      stringsAsFactors = FALSE)
    for (resume in c(FALSE, TRUE)) {
      suppressMessages(amplicanPipeline(
        file.path(folder, "config.csv"), folder, results, knit_reports = FALSE,
        fastqfiles = 1, primer_mismatch = 0, event_filter = FALSE,
        normalize = NULL, batch_size = 1000, continue = resume))
      summary <- data.table::fread(file.path(results, "config_summary.csv"))
      expect_equal(as.data.frame(summary[, names(expected_summary), with = FALSE]),
                   expected_summary, info = paste("continue", resume))
      for (name in c("raw_events.csv", "events_filtered_shifted_normalized.csv")) {
        events <- data.table::fread(file.path(results, "alignments", name))
        expect_equal(
          as.data.frame(events[order(type), names(expected_events), with = FALSE]),
          expected_events, info = paste("continue", resume, "file", name))
      }
    }
  }
  exercise_pipeline()
})

# # after successful amplican pipeline run, replace auto generated
# # system links in .Rmd reports with system.file type of links
# # also make vignettes out of .Rmd reports
# unlink(file.path(results_folder, "reports"), recursive = TRUE, force = FALSE)
#
# # run this when developing and willing to generate new vignettes
# # for visual quality test and for the users to see example reports
# # fix run param file
# runParam <- readLines(file.path(results_folder, "RunParameters.txt"))
# runParam <- gsub(config,
#                  "full/path/to/config/file/that/has/been/used.csv",
#                  runParam)
# cat(runParam,
#     file = file.path(results_folder, "RunParameters.txt"), sep = "\n")
# # fix paths to the summary config
# conf_file <- system.file('extdata',
#                          'results',
#                          'config_summary.csv',
#                          package = 'amplican')
# conf_df <- read.csv(conf_file,
#                     stringsAsFactors = FALSE)
# conf_df$Reverse_Reads_File <- basename(conf_df$Reverse_Reads_File)
# conf_df$Forward_Reads_File <- basename(conf_df$Forward_Reads_File)
# write.csv(conf_df, file = conf_file, row.names = FALSE)
#
# template_names <- c("id_report", "barcode_report", "group_report",
#                     "guide_report", "amplicon_report", "index_report")
# rmd_names <- c("id_report.Rmd", "barcode_report.Rmd", "group_report.Rmd",
#                "guide_report.Rmd", "amplicon_report.Rmd", "index.Rmd")
#
# for (i in seq_along(rmd_names)) {
#   rmd_file_path <- file.path(vignettes_path, paste0("example_", rmd_names[i]))
#   if (file.exists(rmd_file_path)) {
#     file.remove(rmd_file_path)
#   }
#   rmarkdown::draft(file = rmd_file_path,
#                    template = template_names[i],
#                    package = "amplican",
#                    edit = FALSE)
#
#   rmd_content <- readLines(rmd_file_path)
#   settings_end <- which(rmd_content == "---")[2]-1
#   fixed_rmd_content <- c(rmd_content[1:settings_end],
#                          "vignette: >",
#                          paste0("  %\\VignetteIndexEntry{example ",
#                                 sub("report_", "",
#                                     tools::file_path_sans_ext(rmd_names[i])),
#                                 " report}"),
#                          "  %\\VignetteEngine{knitr::rmarkdown}",
#                          "  %\\VignetteEncoding{UTF-8}",
#                          rmd_content[(settings_end+1):length(rmd_content)])
#
#   if (rmd_names[i] == "index.Rmd") {
#     fixed_rmd_content <-
#       c(fixed_rmd_content[1:13],
#         paste0('  links: \"',
#                '1. [Report by id](./example_id_report.html)\\n',
#                '2. [Report by barcode](./example_barcode_report.html)\\n',
#                '3. [Report by group](./example_group_report.html)\\n',
#                '4. [Report by guide](./example_guide_report.html)\\n',
#                '5. [Report by amplicon](./example_amplicon_report.html)\\n\"'),
#         fixed_rmd_content[15:length(fixed_rmd_content)])
#   }
#   cat(fixed_rmd_content, file = rmd_file_path, sep = "\n")
# }
#
# devtools::build_vignettes()

test_that("paired pipelines preserve control normalization and HDR across resume", {
  exercise <- function() {
    folder <- tempfile("paired-truth-")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE))
    cfg <- read.csv(config, stringsAsFactors = FALSE)[c(1, 1), ]
    cfg$ID <- c("T", "C")
    cfg$Barcode <- c("treatment", "control")
    cfg$Control <- c(FALSE, TRUE)
    cfg$Group <- "paired"
    cfg$Forward_Reads <- c("T_R1.fastq", "C_R1.fastq")
    cfg$Reverse_Reads <- c("T_R2.fastq", "C_R2.fastq")
    reference <- toupper(cfg$Amplicon[1])
    origin <- regexpr("[A-Z]", cfg$Amplicon[1])[1]
    position <- origin + 8L
    change <- function(sequence, pos) {
      substr(sequence, pos, pos) <- if (substr(sequence, pos, pos) == "A") "C" else "A"
      sequence
    }
    shared <- change(reference, position)
    donor <- change(change(reference, origin + 12L), origin + 15L)
    partial <- change(reference, origin + 12L)
    cfg$Donor <- c(donor, "")
    deletion <- paste0(substr(reference, 1L, position - 1L), substr(reference, position + 3L, nchar(reference)))
    insertion <- paste0(substr(reference, 1L, position), "A", substr(reference, position + 1L, nchar(reference)))
    primer_dimer <- paste0(substr(reference, 1L, nchar(cfg$Forward_Primer[1])),
      substr(reference, nchar(reference) - nchar(cfg$Reverse_Primer[1]) + 1L, nchar(reference)))
    samples <- list(rep(c(reference, shared, insertion, deletion, donor, partial, primer_dimer), c(2, 3, 5, 7, 4, 6, 2)),
                    rep(c(reference, shared), c(8, 2)))
    for (i in 1:2) {
      for (reverse in c(FALSE, TRUE)) {
        reads <- samples[[i]]
        if (reverse) reads <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reads)))
        lines <- unlist(lapply(seq_along(reads), function(j) c(paste0("@read", j), reads[j], "+", strrep("I", nchar(reads[j])))))
        writeLines(lines, file.path(folder, if (reverse) cfg$Reverse_Reads[i] else cfg$Forward_Reads[i]))
      }
    }
    write.csv(cfg, file.path(folder, "config.csv"), row.names = FALSE)
    expected_summary <- data.frame(ID = c("C", "T"), Reads = c(10, 29), PRIMER_DIMER = c(0, 2), Reads_Filtered = c(10, 27),
      Reads_Del = c(0, 7), Reads_In = c(0, 5), Reads_Edited = c(0, 16), Reads_Frameshifted = c(0, 5), HDR = c(0, 4))
    expected <- data.frame(seqnames = c("C", rep("T", 5)),
      start = c(8, 8, 9, 12, 15, 12), width = c(1, 3, 1, 1, 1, 1),
      type = c("mismatch", "deletion", "insertion", rep("mismatch", 3)), counts = c(2, 7, 5, 4, 4, 6),
      originally = c(substr(reference, position, position), "", "", substr(reference, origin+12L, origin+12L),
                     substr(reference, origin+15L, origin+15L), substr(reference, origin+12L, origin+12L)),
      replacement = c(substr(shared, position, position), "", "A", substr(donor, origin+12L, origin+12L),
                      substr(donor, origin+15L, origin+15L), substr(partial, origin+12L, origin+12L)),
      readType = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE), stringsAsFactors = FALSE)
    canonical <- function(x) {
      x <- as.data.frame(x)
      x <- x[order(x$seqnames, x$start, x$type, x$counts), ]
      rownames(x) <- NULL
      x
    }
    for (phase in c("loose", "strict", "boundary")) {
      results <- file.path(folder, phase)
      dir.create(results)
      wanted <- expected
      if (phase == "boundary") {
        background <- expected[1, ]
        background$seqnames <- "T"
        background$counts <- 3
        wanted <- rbind(wanted, background)
      }
      for (resume in c(FALSE, TRUE)) {
        suppressMessages(amplicanPipeline(file.path(folder, "config.csv"), folder, results,
          knit_reports = FALSE, fastqfiles = 0, primer_mismatch = 0, event_filter = FALSE,
          donor_strict = phase != "loose", donor_mismatch = 0,
          min_freq = if (phase == "boundary") 0.2 else 0.01, batch_size = 1000, continue = resume))
        summary <- data.table::fread(file.path(results, "config_summary.csv"))
        expect_equal(as.data.frame(summary[order(ID), names(expected_summary), with = FALSE]), expected_summary,
                     info = paste(phase, "continue", resume))
        events <- data.table::fread(file.path(results, "alignments", "events_filtered_shifted_normalized.csv"))
        expect_equal(canonical(events[consensus & overlaps, names(wanted), with = FALSE]), canonical(wanted),
                     info = paste(phase, "continue", resume))
      }
    }
  }
  exercise()
})
