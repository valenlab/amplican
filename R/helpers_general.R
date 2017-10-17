seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))


#' Reverse and complement given string or list of strings
#'
#' @keywords internal
#' @param x (string or vector of strings)
#' @return (string or vector of strings) reverse complemented input
#'
revComp <- function(x) {
  return(
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  )
}


#' Get codons for given string - translate
#'
#' @keywords internal
#' @param x (string)
#' @return (string) codons
#'
decode <- function(x) {
  return(
    strsplit(
      as.character(
        Biostrings::translate(
          Biostrings::DNAString(paste0(x, collapse = "")))), "")[[1]])
}


#' amplicon sequence, reverse complemented when needed
#'
#' @keywords internal
#' @param config (data.frame) config table
#' @param id (vector) a vector of id's
#' @return (character) amplicon sequence, reverse complemented if Direction 1
#'
get_amplicon <- function(config, id) {
  amplicon <- as.character(config[which(config$ID == id[1]), "Amplicon"])
  if (config[which(config$ID == id[1]), "Direction"] == 1) {
    # revComp makes upper cases
    groups <- as.data.frame(upperGroups(amplicon))
    old_starts <- groups$start
    groups$start <- nchar(amplicon) - groups$end + 1
    groups$end <- nchar(amplicon) - old_starts + 1
    amplicon <- tolower(revComp(amplicon))
    for (i in seq_len(dim(groups)[1])) { # there may be many UPPER groups
      substr(amplicon, groups$start[i], groups$end[i]) <-
        toupper(substr(amplicon, groups$start[i], groups$end[i]))
    }
  }
  return(amplicon)
}


#' left primer sequence
#'
#' @keywords internal
#' @param config (data.frame) config table
#' @param id (vector) a vector of id's
#' @return (character) left primer sequence
#'
get_left_primer <- function(config, id) {
  if (config[which(config$ID == id[1]), "Direction"] == 1) {
    as.character(config[which(config$ID == id[1]), "Reverse_Primer"])
  } else {
    as.character(config[which(config$ID == id[1]), "Forward_Primer"])
  }
}


#' right primer sequence
#'
#' @keywords internal
#' @param config (data.frame) config table
#' @param id (vector) a vector of id's
#' @return (character) right primer sequence
#'
get_right_primer <- function(config, id) {
  if (config[which(config$ID == id[1]), "Direction"] == 1) {
    revComp(as.character(config[which(config$ID == id[1]), "Forward_Primer"]))
  } else {
    revComp(as.character(config[which(config$ID == id[1]), "Reverse_Primer"]))
  }
}


#' Helper to construct GRanges with additional metadata columns.
#'
#' @keywords internal
#' @param x (IRanges) names(x) indicating read_id
#' @param ID (string)
#' @param type (string)
#' @param score (numeric) scores from the alignments
#' @param strand_info (string) Either '+', '-'
#' @param originally (string) Base pairs on the amplicon.
#' @param replacement (string) Base pairs on the read.
#' @return (GRanges) Object with meta-data
#'
defGR <- function(x,
                  ID,
                  score,
                  strand_info = "+",
                  type = "deletion",
                  originally = "",
                  replacement = "") {

  if (length(x) == 0) return(GenomicRanges::GRanges())
  GenomicRanges::GRanges(
    ranges = x,
    strand = strand_info,
    seqnames = ID,
    originally = as.character(originally),
    replacement = as.character(replacement),
    type = type,
    read_id = names(x),
    score = score
  )
}


#' Cumulative sum to calculate shift
#'
#' @keywords internal
#' @param x (IRanges)
#' @return (numeric vector)
#'
cumsumw <- function(x) {
  c(0, cumsum(IRanges::width(x))[-length(x)])
}


#' This function takes alignments and gives back the events coordinates.
#'
#' @keywords internal
#' @param align (PairwiseAlignmentsSingleSubject)
#' @param ID (string)
#' @param ampl_shift (numeric vector) Shift events additionally by this value.
#' PairwiseAlignmentsSingleSubject returns truncated alignments.
#' @param ampl_len (numeric) Length of the amplicon (subject)
#' @param strand_info (string) Either '+', '-' or default '*'
#' @return (GRanges) Object with meta-data for insertion, deletion, mismatch
#'
getEventInfo <- function(align, ID, ampl_shift, strand_info = "+") {
  if (length(align) == 0) return(GenomicRanges::GRanges())
  if (any(ampl_shift < 1)) stop("Amplicon shift can't be less than 1.")
  scores <- Biostrings::score(align)

  ampl_len <- width(unaligned(subject(align)))
  ampl_end <- end(subject(align))
  ampl_start <- start(subject(align))

  if (strand_info == "+") {
    s_err <- ampl_end >= ampl_len
    start <- ampl_end[!s_err] + ampl_shift
    end <- if (all(s_err)) integer() else ampl_len + ampl_shift - 1
  } else {
    s_err <- ampl_start == 1
    start <- if (all(s_err)) integer() else 1L
    end <- ampl_start[!s_err] + ampl_shift - 2
  }
  sizes <- IRanges::IRanges(start = start, end = end, names = which(!s_err))

  del <- Biostrings::deletion(align)
  ins <- Biostrings::insertion(align)
  mm <- lapply(align, function(x) Biostrings::mismatchSummary(x)$subject)

  ins_sft <- IRanges::shift(ins, IRanges::IntegerList(lapply(ins, cumsumw)))
  ins_seq <- substr(rep(Biostrings::pattern(align),
                        times = sapply(ins_sft, length)),
                    start = unlist(BiocGenerics::start(ins_sft)),
                    stop = unlist(BiocGenerics::end(ins_sft)))
  names(ins) <- seq_along(ins)

  del <- IRanges::shift(del, IRanges::IntegerList(lapply(del, cumsumw)))
  shift_del <- S4Vectors::mendoapply(function(x, y, w) {
    vapply(x, function(x_i, y, w) sum(w[x_i > y]), integer(1), y, w)
  },
  BiocGenerics::start(del), BiocGenerics::end(ins), BiocGenerics::width(ins))
  del <- IRanges::shift(del, -1 * shift_del)
  names(del) <- seq_along(del)

  subj <- as.character(subject(align))
  ins <- IRanges::shift(ins, ampl_shift + ampl_start - 2L)
  del <- IRanges::shift(del, ampl_shift + ampl_start - 2L)

  mism <- lapply(mm, function(x) IRanges::IRanges(x$SubjectPosition, width = 1))
  names(mism) <- seq_along(mism)
  mism <- IRanges::shift(IRanges::IRangesList(mism), ampl_shift - 1L)

  c(defGR(unlist(IRanges::IRangesList(ins)), ID,
          rep(scores, times = sapply(ins, length)),
          strand_info, "insertion", "", ins_seq),
    defGR(unlist(IRanges::IRangesList(del)), ID,
          rep(scores, times = sapply(del, length)),
          strand_info),
    # artificial deletions indicating end of reads
    defGR(sizes, ID, scores[!s_err], strand_info),
    defGR(unlist(mism), ID,
          rep(scores, times = sapply(mism, length)),
          strand_info,
          "mismatch",
          unlist(lapply(mm, function(x) as.character(x$Subject))),
          unlist(lapply(mm, function(x) as.character(x$Pattern)))))
}


#' Detect uppercases as ranges object.
#'
#' For a given string, detect how many groups of uppercases is inside, where
#' are they, and how long they are.
#'
#' For example:
#'    asdkfaAGASDGAsjaeuradAFDSfasfjaeiorAuaoeurasjfasdhfashTTSfajeiasjsf
#'
#' Has 4 groups of uppercases of length 7, 4, 1 and 3.
#' @keywords internal
#' @param candidate (string) A string with the nucleotide sequence.
#' @return (Ranges) A Ranges object with uppercases groups for given candidate
#' string
#'
upperGroups <- function(candidate) {
  return(IRanges::reduce(IRanges::IRanges(
    start = which(stringr::str_detect(strsplit(candidate, "")[[1]],
                                      "[[:upper:]]")),
    width = 1)))
}


#' Reverse complement events that have amplicons with direction 1.
#'
#' @keywords internal
#' @param idR (data.frame) Loaded events.
#' @param cfgT (data.frame) Loaded configuration file.
#' @return (data.frame) Returns input idR, but events for amplicons with
#' direction 1 reverse complemented, "+" and "-" swapped.
#'
flipRanges <- function(idR, cfgT) {

  is_dir <- as.logical(cfgT$Direction)
  to_flip <- cfgT[is_dir, "ID"]
  to_flip <- idR$seqnames %in% to_flip

  if (any(to_flip)) {
    ampl_lengths <- nchar(as.character(cfgT[is_dir, "Amplicon"]))
    ampl_ids <- as.character(cfgT[is_dir, "ID"])
    ids_mapping <- match(idR[to_flip, "seqnames"], ampl_ids)
    ampl_lengths <- ampl_lengths[ids_mapping]

    idR[to_flip, "originally"] <- revComp(idR[to_flip, "originally"])
    idR[to_flip, "replacement"] <- revComp(idR[to_flip, "replacement"])

    strand <- idR[to_flip, "strand"]
    strand_minus <- strand == "-"
    strand[strand == "+"] <- "-"
    strand[strand_minus] <- "+"
    idR[to_flip, "strand"] <- strand

    # mm + del end -> start & start -> end
    # ins start -> start & end -> end
    ins <- idR$type == "insertion"

    old_starts <- idR[to_flip & !ins, "start"]
    idR[to_flip & !ins, "start"] <-
      ampl_lengths[!ins[to_flip]] - idR[to_flip & !ins, "end"] + 1
    idR[to_flip & !ins, "end"] <- ampl_lengths[!ins[to_flip]] - old_starts + 1

    idR[to_flip & ins, "start"] <-
      ampl_lengths[ins[to_flip]] - idR[to_flip & ins, "start"] + 1
    idR[to_flip & ins, "end"] <-
      idR[to_flip & ins, "width"] + idR[to_flip & ins, "start"] - 1
  }
  return(idR)
}


#' Map events to their respective relative coordinates specified with
#' UPPER case.
#'
#' Translate coordinates of GRanges events so that they can be relative to the
#' amplicon. As point zero we assume first left sided UPPER case letter in the
#' amplicon. Be weary that events for amplicons without expected cut sites are
#' filtered. Don't use this function, if you don't have expected cut sites
#' specified and don't use any of the metaplots.
#'
#' @param aln (data.frame) List of events to map to the relative coordinates.
#' @param cfgT (data.frame) config table
#' @return (GRanges) Same as events, but the coordinates are relative to the
#' @export
#' @family analysis steps
#' @examples
#' # example config
#' config <- read.csv(system.file("extdata", "config.csv",
#'                    package = "amplican"))
#' # example events
#' events <- read.csv(system.file("extdata", "results", "alignments",
#'                    "raw_events.csv", package = "amplican"))
#' # make events relative to the UPPER case
#' amplicanMap(events, config)
#'
amplicanMap <- function(aln, cfgT) {
  aln <- GenomicRanges::GRanges(aln)
  no_upper <- FALSE

  for (id in unique(GenomeInfoDb::seqnames(aln))) {
    amplicon <- get_amplicon(cfgT, id)
    zero_point <- upperGroups(amplicon)
    if (length(zero_point) == 0) {
      no_upper <- TRUE
      aln <- aln[GenomeInfoDb::seqnames(aln) != id, ]
      next()
    }
    aln[GenomicRanges::seqnames(aln) == id] <-
      GenomicRanges::shift(aln[GenomeInfoDb::seqnames(aln) == id],
                           shift = -1 * GenomicRanges::start(zero_point)[1])
  }

  if (no_upper) {
    warning("Events for amplicons without UPPER case are filtered.")
  }

  return(aln)
}
