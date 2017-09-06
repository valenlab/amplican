seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))


#' Reverse and complement given string or list of strings
#'
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
  finalGR <- GenomicRanges::GRanges(
    ranges = x,
    strand = S4Vectors::Rle(rep(strand_info, length(x))),
    seqnames = S4Vectors::Rle(rep(ID, length(x))))

  finalGR$originally = as.character(originally)
  finalGR$replacement = as.character(replacement)
  finalGR$type = type
  finalGR$read_id = names(x)
  finalGR$score = score

  return(finalGR)
}


#' Make insertions and deletions relative to the full subject sequence
#'
#' @param x (IRangesList)
#' @param ampl_shift (numeric vector)
#' @param subject (vector of strings) as.character(subject(align))
#' @param strand_info (string) + or -
#' @return (IRangesList) shifted
#'
shift_to_subj <- function(x, ampl_shift, subject, strand_info) {
  if (strand_info == "+") {
    IRanges::shift(x, ampl_shift - 1)
  } else {
    subj_bas <- stringr::str_count(subject, "[ATCG]")
    IRanges::shift(x, ampl_shift - subj_bas)
  }
}


#' Cumulative sum to calculate shift
#'
#' @param x (IRanges)
#' @return (numeric vector)
#'
cumsumw <- function(x) {
  c(0, cumsum(IRanges::width(x))[-length(x)])
}


#' This function takes alignments and gives back the events coordinates.
#'
#' @param align (PairwiseAlignmentsSingleSubject)
#' @param ID (string)
#' @param ampl_shift (numeric vector) Shift events additionally by this value.
#' PairwiseAlignmentsSingleSubject returns truncated alignments.
#' @param ampl_len (numeric) Length of the amplicon (subject)
#' @param strand_info (string) Either '+', '-' or default '*'
#' @return (GRanges) Object with meta-data for insertion, deletion, mismatch
#'
getEventInfo <- function(align, ID, ampl_shift, ampl_len, strand_info = "+") {
  if (length(align) == 0) return(GenomicRanges::GRanges())
  if (any(ampl_shift < 1)) stop("Amplicon shift can't be less than 1.")
  scores <- Biostrings::score(align)
  sizes <- nchar(align)
  if (strand_info == "+") {
    s_err <- sizes + ampl_shift - 1 >= ampl_len
    sizes <- if (all(s_err)) IRanges::IRanges() else
      IRanges::IRanges(start = sizes[!s_err] + ampl_shift,
                       end = ampl_len)
    names(sizes) <- which(!s_err)
  } else {
    s_err <- sizes + abs(ampl_shift - ampl_len) >= ampl_len
    sizes <- if (all(s_err)) IRanges::IRanges() else
      IRanges::IRanges(start = 1, end = ampl_shift - sizes[!s_err])
    names(sizes) <- which(!s_err)
  }

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
  # make deletions to be relative to the subject
  del <- mapply(function(x, y){
    if (length(y) > 0) { # shift when insertions are before deletions
      for (i in seq_along(x)) {
        ins_before <- BiocGenerics::start(x)[i] > BiocGenerics::start(y)
        if (any(ins_before)) {
          x[i] <- IRanges::shift(x[i], -1 *
                                   sum(BiocGenerics::width(y[ins_before])))
        }
      }
    }
    x
  }, del, ins)
  names(del) <- seq_along(del)

  subj <- as.character(subject(align))
  ins <- shift_to_subj(ins, ampl_shift, subj, strand_info)
  del <- shift_to_subj(IRanges::IRangesList(del), ampl_shift, subj, strand_info)

  mism <- lapply(mm, function(x) IRanges::IRanges(x$SubjectPosition, width = 1))
  names(mism) <- seq_along(mism)

  c(defGR(unlist(IRanges::IRangesList(ins)), ID,
          rep(scores, times = sapply(ins, length)),
          strand_info, "insertion", "", ins_seq),
    defGR(unlist(IRanges::IRangesList(del)), ID,
          rep(scores, times = sapply(del, length)),
          strand_info),
    # artificial deletions indicating end of reads
    defGR(sizes, ID, scores[!s_err], strand_info),
    defGR(unlist(IRanges::IRangesList(mism)), ID,
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

    old_starts <- idR[to_flip, "start"]
    idR[to_flip, "start"] <- ampl_lengths - idR[to_flip, "end"] + 1
    idR[to_flip, "end"] <- ampl_lengths - old_starts + 1

    strand <- idR[to_flip, "strand"]
    strand_minus <-strand == "-"
    strand[strand == "+"] <- "-"
    strand[strand_minus] <- "+"
    idR[to_flip, "strand"] <- strand
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
    warning(paste0("Events for amplicons without UPPER case are filtered."))
  }

  return(aln)
}
