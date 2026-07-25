#' Automated analysis of CRISPR experiments.
#'
#' Main goals:
#' \enumerate{
#' \item Flexible pipeline for analysis of the CRISPR Mi-Seq or Hi-Seq data.
#' \item Compatible with GRanges and data.table style.
#' \item Precise quantification of mutation rates.
#' \item Prepare automatic reports as .Rmd files that are flexible
#' and open for manipulation.
#' \item Provide specialized plots for deletions, insertions, mismatches,
#' variants, heterogeneity of the reads.
#' }
#'
#' To learn more about amplican, start with the vignettes:
#' \code{browseVignettes(package = "amplican")}
#'
#' @docType package
#' @name amplican
#' @useDynLib amplican
#'
#' @import Rcpp ggthemes knitr methods data.table
#' @rawNamespace import(BiocGenerics, except = Position)
#' @importFrom Rcpp sourceCpp
#' @importFrom IRanges coverage
#' @importFrom Biostrings DNAString DNAStringSet extractAt quality
#' @importFrom pwalign pairwiseAlignment writePairwiseAlignments pattern subject unaligned compareStrings
#'
"_PACKAGE"
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0(
      "version: ", utils::packageVersion("amplican"), "\n",
      "Please consider supporting this software by citing:\n\n",
      "Labun et al. 2019\n",
      "Accurate analysis of genuine CRISPR editing events with ampliCan.\n",
      "Genome Res. 2019 Mar 8\n",
      "doi: 10.1101/gr.244293.118\n"
    )
  )
}

amplicanPipe <- function(min_freq_default) {
  function(
    config, fastq_folder, results_folder, knit_reports = TRUE,
    write_alignments_format = "None", average_quality = 30,
    min_quality = 0, filter_n = FALSE, batch_size = 1e7, use_parallel = FALSE,
    scoring_matrix = pwalign::nucleotideSubstitutionMatrix(
      match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA"
    ),
    gap_opening = 25, gap_extension = 0, fastqfiles = 0.5,
    primer_mismatch = 2,
    donor_mismatch = 3, donor_strict = FALSE,
    PRIMER_DIMER = 30,
    event_filter = TRUE, max_remove_filterLQR = 0.25, cut_buffer = 5,
    promiscuous_consensus = TRUE, normalize = c("guideRNA", "Group"),
    min_freq = min_freq_default,
    continue = TRUE, sample = 0, seed = 0
  ) {
    config <- normalizePath(config)
    fastq_folder <- normalizePath(fastq_folder)
    results_folder <- normalizePath(results_folder)

    message("Checking write access...")
    checkFileWriteAccess(results_folder)

    if (!continue) {
      message("continue is FALSE, removeing contents of results folder.")
      unlink(results_folder, recursive = TRUE)
      dir.create(results_folder, showWarnings = FALSE)
    } else {
      temp_files <- list.files(results_folder,
        pattern = "\\.temp$",
        full.names = TRUE, recursive = TRUE
      )
      if (length(temp_files) > 0) {
        message("Cleaning up ", length(temp_files), " stale .temp files...")
        file.remove(temp_files)
      }
    }
    resultsFolder <- file.path(results_folder, "alignments")
    if (!dir.exists(resultsFolder)) {
      dir.create(resultsFolder)
    }
    rds_file <- file.path(resultsFolder, "AlignmentsExperimentSet.rds")
    re_file <- file.path(resultsFolder, "raw_events.csv")
    un_file <- file.path(resultsFolder, "unassigned_reads.csv")
    bd_file <- file.path(results_folder, "barcode_reads_filters.csv")
    cfgT_temp_file <- file.path(resultsFolder, "experiment_data.rds")

    if (file.exists(rds_file)) {
      message("Loading alignments...")
      aln <- readRDS(rds_file)
      if (!"None" %in% write_alignments_format) {
        for (frmt in write_alignments_format) {
          aln_file_frmt <- file.path(
            resultsFolder,
            paste0("alignments.", frmt)
          )
          if (!file.exists(aln_file_frmt)) {
            aln_file_frmt_temp <- paste0(aln_file_frmt, ".temp")
            writeAlignments(aln, aln_file_frmt_temp, frmt)
            file.rename(aln_file_frmt_temp, aln_file_frmt)
          }
        }
      }
      if (!file.exists(un_file)) {
        message("Saving unassigned sequences...")
        unData <- unassignedData(aln)
        if (!is.null(unData)) {
          un_file_temp <- paste0(un_file, ".temp")
          data.table::fwrite(unData, un_file_temp)
          file.rename(un_file_temp, un_file)
        }
      }
      if (!file.exists(bd_file)) {
        message("Saving barcode statistics...")
        bd_file_temp <- paste0(bd_file, ".temp")
        data.table::fwrite(barcodeData(aln), bd_file_temp)
        file.rename(bd_file_temp, bd_file)
      }
      cfgT <- experimentData(aln)
      if (!file.exists(re_file)) {
        message("Translating alignments into events...")
        aln <- extractEvents(aln, use_parallel = use_parallel)
        message("Saving complete events - unfiltered...")
        re_file_temp <- paste0(re_file, ".temp")
        data.table::fwrite(aln, re_file_temp)
        file.rename(re_file_temp, re_file)
        message("Saved complete events - unfiltered.")
        aln <- data.table::as.data.table(aln)
      } else {
        message("Reading complete events - unfiltered.")
        aln <- data.table::fread(re_file)
      }
    } else {
      tempFolder <- file.path(resultsFolder, "temp")
      if (!dir.exists(tempFolder)) dir.create(tempFolder)

      if (!file.exists(re_file) || !file.exists(cfgT_temp_file)) {
        message("Making alignments in chunked mode...")
        aln_paths <- amplicanAlign(
          config = config,
          fastq_folder = fastq_folder,
          use_parallel = use_parallel,
          average_quality = average_quality,
          batch_size = batch_size,
          scoring_matrix = scoring_matrix,
          gap_opening = gap_opening,
          gap_extension = gap_extension,
          min_quality = min_quality,
          filter_n = filter_n,
          fastqfiles = fastqfiles,
          primer_mismatch = primer_mismatch,
          donor_mismatch = donor_mismatch,
          donor_strict = donor_strict,
          temp_folder = tempFolder,
          sample = sample,
          seed = seed
        )

        message("Extracting events and compiling statistics...")
        p <- if (!use_parallel) BiocParallel::SerialParam() else BiocParallel::bpparam()

        chunk_results <- BiocParallel::bplapply(aln_paths, function(path) {
          chunk_aln <- readRDS(path)
          csv_file <- gsub("_aln.rds", "_events.csv", path)

          if (!file.exists(csv_file)) {
            if (!"None" %in% write_alignments_format) {
              for (frmt in write_alignments_format) {
                frmt_file <- paste0(path, ".", frmt)
                frmt_temp <- paste0(frmt_file, ".temp")
                writeAlignments(chunk_aln, frmt_temp, frmt)
                file.rename(frmt_temp, frmt_file)
              }
            }
            chunk_events <- extractEvents(chunk_aln, use_parallel = FALSE)
            csv_file_temp <- paste0(csv_file, ".temp")
            data.table::fwrite(chunk_events, csv_file_temp)
            file.rename(csv_file_temp, csv_file)
          }

          return(list(
            events_file = csv_file,
            unData = unassignedData(chunk_aln),
            bdData = barcodeData(chunk_aln),
            cfgT = experimentData(chunk_aln)
          ))
        }, BPPARAM = p)

        if (!"None" %in% write_alignments_format) {
          for (frmt in write_alignments_format) {
            aln_file_frmt <- file.path(resultsFolder, paste0("alignments.", frmt))
            aln_file_frmt_temp <- paste0(aln_file_frmt, ".temp")
            if (file.exists(aln_file_frmt_temp)) unlink(aln_file_frmt_temp)
            for (path in aln_paths) {
              chunk_frmt <- paste0(path, ".", frmt)
              if (file.exists(chunk_frmt)) {
                file.append(aln_file_frmt_temp, chunk_frmt)
                file.remove(chunk_frmt)
              }
            }
            file.rename(aln_file_frmt_temp, aln_file_frmt)
          }
        }

        unData <- data.table::rbindlist(lapply(chunk_results, function(x) x$unData), fill = TRUE)
        if (!is.null(unData) && nrow(unData) > 0) {
          message("Saving unassigned sequences...")
          un_file_temp <- paste0(un_file, ".temp")
          data.table::fwrite(unData, un_file_temp)
          file.rename(un_file_temp, un_file)
        }

        message("Saving barcode statistics...")
        bdData <- data.table::rbindlist(lapply(chunk_results, function(x) x$bdData), fill = TRUE)
        bd_file_temp <- paste0(bd_file, ".temp")
        data.table::fwrite(bdData, bd_file_temp)
        file.rename(bd_file_temp, bd_file)

        cfgT_chunks <- lapply(chunk_results, function(x) x$cfgT)
        cfgT <- as.data.frame(data.table::rbindlist(cfgT_chunks, fill = TRUE))
        original_config <- data.frame(data.table::fread(config))
        cfgT <- cfgT[match(original_config$ID, cfgT$ID), ]
        cfgT_temp_file_writing <- paste0(cfgT_temp_file, ".temp")
        saveRDS(cfgT, cfgT_temp_file_writing)
        file.rename(cfgT_temp_file_writing, cfgT_temp_file)

        message("Saving complete events - unfiltered...")
        aln <- data.table::rbindlist(lapply(chunk_results, function(x) data.table::fread(x$events_file, na.strings = "")), fill = TRUE)
        re_file_temp <- paste0(re_file, ".temp")
        data.table::fwrite(aln, re_file_temp)
        file.rename(re_file_temp, re_file)
        message("Saved complete events - unfiltered.")
      } else {
        message("Reading complete events - unfiltered.")
        aln <- data.table::fread(re_file, na.strings = "")
        cfgT <- readRDS(cfgT_temp_file)
      }
    }

    logFileName <- file.path(results_folder, "RunParameters.txt")
    if (!file.exists(logFileName)) {
      message("Saving parameters...")
      logFileNameTemp <- paste0(logFileName, ".temp")
      logFileConn <- file(logFileNameTemp, open = "at")
      writeLines(c(
        paste("amplican Version:   ", utils::packageVersion("amplican")),
        paste("Config file:        ", config),
        paste("Average Quality:    ", average_quality),
        paste("Minimum Quality:    ", min_quality),
        paste("Filter N-reads:     ", filter_n),
        paste("Batch size:         ", batch_size),
        paste("Write Alignments:   ", toString(write_alignments_format)),
        paste("Fastq files Mode:   ", fastqfiles),
        paste("Gap Opening:        ", gap_opening),
        paste("Gap Extension:      ", gap_extension),
        paste("Consensus:          ", promiscuous_consensus),
        paste("Normalize:          ", toString(normalize)),
        paste("PRIMER DIMER buffer:", PRIMER_DIMER),
        paste("Event filter:       ", event_filter),
        paste("Max remove filterLQR:", max_remove_filterLQR),
        paste("Cut buffer:", cut_buffer),
        "Scoring Matrix:"
      ), logFileConn)
      utils::write.csv(scoring_matrix, logFileConn, quote = FALSE, row.names = TRUE)
      close(logFileConn)
      file.rename(logFileNameTemp, logFileName)
    }

    data.table::setDT(cfgT)
    data.table::setDT(aln)
    if (nrow(aln) == 0) {
      stop(
        "There are no events.",
        "Check whether you have correct primers in the config file."
      )
    }

    efs_file <- file.path(resultsFolder, "events_filtered_shifted.csv")
    cs_file <- file.path(results_folder, "config_summary.csv")
    if (!file.exists(efs_file) | !file.exists(cs_file)) {
      aln$overlaps <- amplicanOverlap(aln, cfgT, cut_buffer = cut_buffer)
      aln$consensus <- if (fastqfiles <= 0.5) {
        amplicanConsensus(aln, cfgT, promiscuous = promiscuous_consensus)
      } else {
        TRUE
      }

      # filter events overlapping primers
      eOP <- findEOP(aln, cfgT)
      aln <- aln[!eOP, ]

      # find PRIMER DIMERS
      PD <- findPD(aln, cfgT, PRIMER_DIMER = PRIMER_DIMER)

      # summarize how many PRIMER DIMER reads per ID
      onlyPD <- aln[PD, ]
      onlyPD <- unique(onlyPD, by = c("seqnames", "read_id"))
      onlyPD <- data.table::as.data.table(onlyPD)
      summaryPD <- onlyPD[, list(counts = sum(counts)), by = c("seqnames")]
      cfgT$PRIMER_DIMER <- 0
      cfgT$PRIMER_DIMER[match(summaryPD$seqnames, cfgT$ID)] <- summaryPD$counts

      # apply filter - remove all events that come from PD infected reads
      aln <- aln[!onlyPD, on = list(seqnames, read_id)]

      # alignment event filter
      cfgT$Low_Score <- 0
      if (event_filter) {
        data.table::setkey(aln, seqnames)
        bad_reads_list <- lapply(seq_len(nrow(cfgT)), function(i) {
          aln_id <- aln[.(cfgT$ID[i]), nomatch = NULL]
          if (nrow(aln_id) == 0 || cfgT$Donor[i] != "") {
            return(NULL)
          }
          onlyBR <- aln_id[findLQR(aln_id, seed = seed,
                                   max_remove_filterLQR = max_remove_filterLQR), ]
          onlyBR <- unique(onlyBR, by = "read_id")

          if (nrow(onlyBR) > 0) {
            data.table::set(cfgT, i, "Low_Score", sum(onlyBR$counts))
            return(onlyBR[, c("seqnames", "read_id"), with = FALSE])
          }
          return(NULL)
        })

        bad_reads <- data.table::rbindlist(bad_reads_list)
        if (nrow(bad_reads) > 0) {
          aln <- aln[!bad_reads, on = c("seqnames", "read_id")]
        }
      }
      cfgT$Reads_Filtered <- cfgT$Reads - cfgT$PRIMER_DIMER - cfgT$Low_Score

      # shift to relative (most left UPPER case is position 0)
      message("Shifting events as relative...")
      data.table::setDF(aln)
      aln <- data.frame(amplicanMap(aln, cfgT), stringsAsFactors = FALSE)
      message("Saving shifted events - filtered...")
      efs_file_temp <- paste0(efs_file, ".temp")
      data.table::fwrite(aln, efs_file_temp)
      file.rename(efs_file_temp, efs_file)
      message("Saved shifted events - filtered.")
      # revert guides to 5'-3'
      cfgT$guideRNA[cfgT$Direction] <- revComp(cfgT$guideRNA[cfgT$Direction])
      cs_file_temp <- paste0(cs_file, ".temp")
      data.table::fwrite(cfgT, cs_file_temp, nThread = 1)
      file.rename(cs_file_temp, cs_file)
    } else {
      message("Reading shifted events - filtered.")
      aln <- fread(efs_file, na.strings = "")
      cfgT <- fread(cs_file)
    }

    efsn_file <- file.path(
      resultsFolder,
      "events_filtered_shifted_normalized.csv"
    )
    if (!file.exists(efsn_file)) {
      message("Normalizing events...")
      # we remove all N as they are just noise from poor sequencing
      aln <- aln[!is.na(aln$replacement) & aln$replacement != "N", ]
      aln <- amplicanNormalize(aln, cfgT, min_freq = min_freq, add = normalize)
      message("Saving normalized events...")
      efsn_file_temp <- paste0(efsn_file, ".temp")
      data.table::fwrite(aln, efsn_file_temp)
      file.rename(efsn_file_temp, efsn_file)
      message("Saved normalized events.")
    } else {
      message("Reading normalized events.")
      aln <- fread(efsn_file, na.strings = "")
    }

    if (donor_strict) {
      message("HDR detection with strict search...")
      aln <- is_hdr_strict(
        aln, cfgT,
        scoring_matrix, gap_opening,
        gap_extension,
        donor_mismatch = donor_mismatch,
        cut_buffer = cut_buffer
      )
      message("Saving normalized events with HDR...")
      efsn_file_temp <- paste0(efsn_file, ".temp")
      data.table::fwrite(aln, efsn_file_temp)
      file.rename(efsn_file_temp, efsn_file)
      message("Saved normalized events with HDR.")
    }

    # summarize
    cfgT <- amplicanSummarize(aln[aln$consensus & aln$overlaps, ], cfgT)
    cs_file_temp <- paste0(cs_file, ".temp")
    data.table::fwrite(
      cfgT[, c(
        "ID", "Barcode", "Forward_Reads_File", "Reverse_Reads_File",
        "Group", "guideRNA", "Found_Guide", "Control", "Forward_Primer",
        "Reverse_Primer", "Direction", "Amplicon", "Donor", "fwdPrPosEnd",
        "rvePrPos", "Reads", "PRIMER_DIMER", "Low_Score",
        "Reads_Filtered", "Reads_Del", "Reads_In",
        "Reads_Edited", "Reads_Frameshifted", "HDR"
      )], cs_file_temp
    )
    file.rename(cs_file_temp, cs_file)

    # reports
    reportsFolder <- file.path(results_folder, "reports")
    if (dir.exists(reportsFolder)) {
      unlink(reportsFolder, recursive = TRUE)
      dir.create(reportsFolder)
    } else {
      dir.create(reportsFolder)
    }

    message(paste0(
      "Making reports... \nDue to high quality ",
      "figures, it is time consuming. Use .Rmd templates for ",
      "more control."
    ))
    amplicanReport(results_folder,
      knit_reports = knit_reports,
      cut_buffer = cut_buffer,
      report_files = file.path(
        reportsFolder,
        c(
          "id_report",
          "barcode_report",
          "group_report",
          "guide_report",
          "amplicon_report",
          "index"
        )
      )
    )
    message("Finished.")
    invisible(results_folder)
  }
}


#' Wraps main package functionality into one function.
#'
#' amplicanPipeline is convenient wrapper around all functionality of the
#' package with the most robust settings. It will generate all results in the
#' \code{result_folder} and also knit prepared reports into 'reports' folder.
#' @param results_folder (string) Where do you want to store results?
#' The package will create files in that folder so make sure you have writing
#' permissions.
#' @param config (string) The path to your configuration file. For example:
#' \code{system.file("extdata", "config.txt", package = "amplican")}.
#' Configuration file can contain additional columns, but first 11 columns
#' have to follow the example config specification.
#' @param fastq_folder (string) Path to FASTQ files. If not specified,
#' FASTQ files should be in the same directory as config file.
#' @param knit_reports (boolean) whether function should "knit" all
#' reports automatically for you (it is time consuming, be patient), when false
#' reports will be prepared, but not knitted
#' @param use_parallel (boolean) Set to TRUE, if you have registered
#' multicore back-end.
#' @param average_quality (numeric) The FASTQ file have a quality for each
#' nucleotide, depending on sequencing technology there exist many formats.
#' This package uses \code{\link[ShortRead]{readFastq}} to parse the reads.
#' If the average quality of the reads fall below value of
#' \code{average_quality} then sequence is filtered. Default is 0.
#' @param min_quality (numeric)  Similar as in average_quality, but depicts
#' the minimum quality for ALL nucleotides in given read. If one of nucleotides
#' has quality BELLOW \code{min_quality}, then the sequence is filtered.
#' Default is 20.
#' @param filter_n (boolean)  Whether to filter out reads containing N base.
#' @param batch_size (numeric) How many reads to analyze at a time? Needed for
#' filtering of large fastq files.
#' @param write_alignments_format (character vector) Whether
#' \code{amplicanPipeline} should write alignments results to separate files.
#' Alignments are also saved as chunked .rds objects inside the `temp` folder
#' to conserve memory.
#' Possible options are:
#' \describe{
#'  \item{"fasta"}{ outputs alignments in fasta format where header indicates
#' experiment ID, read id and number of reads}
#'  \item{"txt"}{ simple format, read information followed by forward read and
#'  amplicon sequence followed by reverse read with its amplicon sequence
#'  eg.: \cr
#' \preformatted{
#' ID: ID_1 Count: 7
#' ACTGAAAAA--------
#' ACTG-----ACTGACTG
#'
#' ------G-ACTG
#' ACTGACTGACTG
#' }}
#' \item{"None"}{ Don't write any alignments to files.}
#' \item{c("fasta", "txt")}{ There are also possible combinations of
#' above formats, pass a vector to get alignments in multiple formats.}
#' }
#' @param scoring_matrix (matrix) Default is 'NUC44'. Pass desired matrix using
#' \code{\link[pwalign]{nucleotideSubstitutionMatrix}}.
#' @param gap_opening (numeric) The opening gap score.
#' @param gap_extension (numeric) The gap extension score.
#' @param fastqfiles (numeric) Normally you want to use both FASTQ files. But in
#' some special cases, you may want to use only the forward file, or only
#' the reverse file. Possible options:
#' \describe{
#'  \item{0}{ Use both FASTQ files.}
#'  \item{0.5}{ Use both FASTQ files, but only for one of the reads (forward or
#'  reverse) is required to have primer perfectly matched to sequence - eg. use
#'  when reverse reads are trimmed of primers, but forward reads have forward
#'  primer in the sequence.}
#'  \item{1}{ Use only the forward FASTQ file.}
#'  \item{2}{ Use only the reverse FASTQ file.}
#' }
#' @param primer_mismatch (numeric) Decide how many mismatches are allowed
#' during primer matching of the reads, that groups reads by experiments.
#' When \code{primer_mismatch = 0} no mismatches are allowed, which can increase
#' number of unasssigned read.
#' @param donor_mismatch (numeric) Maximum number of width-1 events (single-base
#' mismatches or single-base indels) that are allowed to overlap the
#' donor-vs-amplicon event positions. Only events whose coordinates fall within
#' the donor-event window are counted; mismatches or indels elsewhere in the read
#' are invisible to this threshold. The higher the value the more permissive the
#' HDR calling. Set to 0 to require the donor region to match perfectly (not
#' recommended in practice due to sequencing error rate). Only used when a donor
#' template is provided and \code{donor_strict = FALSE}.
#' @param donor_strict (logical) Applies the strict event-presence algorithm for
#' HDR detection via \code{\link{is_hdr_strict}}. When \code{TRUE}, only reads
#' that contain \emph{every} event distinguishing the donor from the amplicon
#' (matched exactly by coordinate, type, and sequence) are counted as HDR.
#' When \code{donor_mismatch = Inf} (the default), additional events elsewhere
#' in the read (noise) do \strong{not} disqualify a read.  When
#' \code{donor_mismatch} is finite (e.g. 0), reads carrying more than
#' \code{donor_mismatch} extra consensus events beyond the donor events are
#' rejected.
#' Use when your reads should span the full donor-event window and you want
#' exact event-coordinate matching. More time-consuming than the default.
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @param event_filter (logical) Whether detection of offtarget reads,
#' should be enabled. Defaults to \code{FALSE} since at high editing rates the
#' unsupervised detector can mistake the edited majority for off-targets.
#' @param max_remove_filterLQR (numeric) Fraction of reads (0-1) above which the
#' off-target filter (\code{\link{findLQR}}) disables itself with a warning.
#' Only relevant when \code{event_filter = TRUE}.
#' @param cut_buffer The number of bases by which extend expected cut sites
#' (specified as UPPER case letters in the amplicon) in 5' and 3' directions.
#' @param promiscuous_consensus (boolean) Whether rules of
#' \code{\link{amplicanConsensus}} should be \code{promiscuous}. When
#' promiscuous, we allow indels that have no confirmation on the other strand.
#' @param normalize (character vector or NULL) If column 'Control' in config table
#' has all FALSE/0 values then normalization is skipped. Otherwise,
#' normalization is strict, which means events that are
#' found in 'Control' TRUE group will be removed in 'Control' FALSE group.
#' This parameter by default uses columns 'guideRNA' and 'Group' to impose
#' additional restrictions on normalized events eg. only events created by the
#' same 'guideRNA' in the same 'Group' will be normalized. Pass \code{NULL}
#' to skip normalization entirely even when Control rows are present. Pass
#' \code{c()} to normalize globally with no group stratification.
#' @param min_freq (numeric) All events below this frequency are treated as
#' sequencing errors and rejected. This parameter is used during normalization
#' through \code{\link{amplicanNormalize}}.
#' @param continue (boolean) Default TRUE, decides whether to continue failed
#' ampliCan runs. In case of FALSE, all contents in `results` folder will
#' be removed.
#' @param sample (numeric) if user specifies `sample` > 0, we will sample only `sample` reads instead of reading full file, then we will process as normal.
#' @param seed (numeric) random seed used for sampling, only used when `sample` > 0.
#' @include amplicanAlign.R amplicanReport.R
#' @return (invisible) results_folder path
#' @export
#' @family analysis steps
#' @examples
#' # path to example config file
#' config <- system.file("extdata", "config.csv", package = "amplican")
#' # path to example fastq files
#' fastq_folder <- system.file("extdata", package = "amplican")
#' # output folder
#' results_folder <- tempdir()
#'
#' # full analysis, not knitting files automatically
#' amplicanPipeline(config, fastq_folder, results_folder, knit_reports = FALSE)
#'
# config <- system.file("extdata", "config.csv", package = "amplican")
# fastq_folder <- system.file("extdata", package = "amplican")
# results_folder <- tempdir()
# knit_reports = TRUE
# write_alignments_format = "txt"
# average_quality = 30
# min_quality = 0
# use_parallel = FALSE
# scoring_matrix = pwalign::nucleotideSubstitutionMatrix(
#   match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA")
# gap_opening = 25
# gap_extension = 0
# fastqfiles = 0.5
# PRIMER_DIMER = 30
# event_filter = FALSE
# max_remove_filterLQR = 0.25
# cut_buffer = 5
# primer_mismatch = 1
# promiscuous_consensus = TRUE
# normalize = c("guideRNA", "Group")
# donor_mismatch = 3
# min_freq = 0.01
amplicanPipeline <- amplicanPipe(0.01)


#' Wraps main package functionality into one function.
#'
#' amplicanPipelineIndexHopping is identical as amplicanPipeline except that
#' default \code{min_freq} threshold is set to 0.15. Setting this threshold
#' higher will decrease risks of inadequate normalization in cases of potential
#' Index Hopping, potentially decreasing precision of true editing rate calling.
#' Index Hopping can be mitigated with use  of unique dual indexing pooling
#' combinations. However, in cases when you might expect Index Hopping to occur
#' you should use this function instead of amplicanPipeline.
#'
#' \code{result_folder} and also knit prepared reports into 'reports' folder.
#' @inheritParams amplicanPipeline
#' @include amplicanAlign.R amplicanReport.R
#' @return (invisible) results_folder path
#' @export
#' @family analysis steps
#'
amplicanPipelineConservative <- amplicanPipe(0.15)
