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


#' Transform aligned strings into GRanges representation of events.
#'
#' Transforms aligned strings into GRanges representation with
#' events of deletions, insertions and mismatches. Subject should come from
#' one amplicon sequence, after alignment to many sequences (patterns).
#'
#' @param pattern (character) Aligned pattern.
#' @param subject (character) Aligned subject.
#' @param scores (integer) Alignment scores of the pattern and subject.
#' @param ID (character) Will be used as seqnames of output GRanges.
#' @param ampl_shift (numeric) Possible shift of the amplicons.
#' @param ampl_start (numeric) Real amplicon starts.
#' Biostrings::pairwiseAlignment clips alignments, therefore to output
#' GRanges relative to the amplicon sequence (subject) ranges have to be
#' shifted.
#' @param strand_info (character) Strands to assign.
#' @return (GRanges) Same as events.
#' @export
#'
getEvents <- function(pattern, subject, scores, ID = "NA", ampl_shift = 1L,
                      ampl_start = 1L, strand_info = "+") {
  pattern <- DNAStringSet(pattern)
  subject <- DNAStringSet(subject)
  read_id <- seq_along(pattern)
  comparison <- compareStrings(pattern, subject)
  comparison <- IRanges::RleList(strsplit(comparison, split = ""))
  mm <- IRanges::IRangesList(comparison == "?")
  del <- IRanges::IRangesList(comparison[comparison != "+"] == "-")
  ins <- IRanges::IRangesList(comparison == "+")
  ins_rep <- unlist(extractAt(pattern, ins), use.names = FALSE)
  mm_rep <- unlist(extractAt(pattern, mm), use.names = FALSE)
  mm_lw <- width(mm_rep) > 1
  mm_rep <- c(as.character(mm_rep[!mm_lw]),
              unlist(strsplit(as.character(mm_rep[mm_lw]), split = "")))
  mm_org <- unlist(extractAt(subject, mm), use.names = FALSE)
  mm_org <- c(as.character(mm_org[!mm_lw]),
              unlist(strsplit(as.character(mm_org[mm_lw]), split = "")))
  ins <- IRanges::shift(ins, IRanges::IntegerList(
    lapply(ins, function(x) -1L * cumsumw(x))))
  mm <- IRanges::IRangesList(comparison[comparison != "+"] == "?")
  mm <- IRanges::shift(mm, ampl_shift + ampl_start - 2L)
  names(mm) <- read_id
  mm <- unlist(mm, use.names = TRUE)
  mm_l <- mm[IRanges::width(mm) > 1]
  mm <- mm[IRanges::width(mm) == 1]
  mm_ln <- names(mm_l)
  mm_l <- IRanges::tile(mm_l, width = 1L)
  names(mm_l) <- mm_ln
  mm_l <- unlist(mm_l, use.names = TRUE)
  mm <- c(mm, mm_l)
  ins <- IRanges::shift(ins, ampl_shift + ampl_start - 2L)
  names(ins) <- read_id
  del <- IRanges::shift(del, ampl_shift + ampl_start - 2L)
  names(del) <- read_id

  c(defGR(unlist(ins, use.names = TRUE), ID,
          rep(scores, times = lengths(ins)),
          strand_info, "insertion", "",
          as.character(ins_rep)),
    defGR(unlist(del, use.names = TRUE), ID,
          rep(scores, times = sapply(del, length)),
          strand_info),
    defGR(mm, ID, scores[as.integer(names(mm))],
          strand_info, "mismatch", as.character(mm_org),
          as.character(mm_rep)))
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
  scores <- score(align)
  subj <- subject(align)
  pat <- pattern(align)

  ampl_len <- width(unaligned(subj))
  ampl_end <- end(subj)
  ampl_start <- start(subj)

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
  ranges <- getEvents(
    pat, subj, scores = scores, ID = ID, ampl_shift = ampl_shift,
    ampl_start = ampl_start, strand_info = strand_info)

  # artificial deletions indicating end of reads
  scores <- scores[!s_err]
  if (!S4Vectors::isEmpty(sizes)) {
    ranges <- c(
      ranges, defGR(sizes[ampl_len > 0], ID, scores[ampl_len > 0], strand_info))
  }
  return(ranges)
}


#' Transform extended CIGAR strings into GRanges.
#'
#' Transform extended CIGAR strings into GRanges representation with
#' events of deletions, insertions and mismatches. Use with caution as function
#' is being tested.
#'
#' @param cigars (character) Extended CIGARS., , query_seq, ref, read_id, mapq, seqnames, strands
#' @param aln_pos_start (integer) Pos of CIGARS.
#' @param query_seq (character) Aligned query sequences.
#' @param ref (character) Reference sequences used for alignment.
#' @param read_id (numeric) Read id for assignment for each of the CIGARS.
#' @param mapq (numeric) Maping scores.
#' @param seqnames (character) Names of the sequences, potentially ids of
#' the reference sequences.
#' @param counts (integer) Vector of cigar counts, if data collapsed.
#' @param strands (character) Strands to assign.
#' @return (GRanges) Same as events.
#' @export
#'
cigarsToEvents <- function(cigars, aln_pos_start, query_seq, ref, read_id, mapq,
                           seqnames, strands, counts) {

  if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
    stop("Install GenomicAlignments before calling this function.")
  }

  ids <- seq_along(cigars)
  # INS
  ins <- GenomicAlignments::cigarRangesAlongQuerySpace(cigars, ops = "I")
  repl <- Biostrings::extractAt(DNAStringSet(query_seq), ins)
  ins <- GenomicAlignments::cigarRangesAlongPairwiseSpace(cigars, ops = "I")
  csum <- lapply(ins, function(x) -1L * cumsumw(x))
  csum <- IRanges::IntegerList(csum)

  names(ins) <- ids
  names(csum) <- ids
  csum <- unlist(csum, use.names = TRUE)
  ins <- unlist(ins, use.names = TRUE) # empty ranges are droped out
  csum <- csum[names(csum) %in% names(ins)]
  ins <- IRanges::shift(ins, csum) # shift by cumsum of ins

  iids <- as.integer(names(ins))
  ins <- if (length(ins) > 0) {
    GenomicRanges::GRanges(
      seqnames = seqnames[iids],
      ranges = ins,
      strand = strands[iids],
      originally = "",
      replacement = as.character(unlist(repl, use.names = FALSE)),
      type = "insertion",
      read_id = read_id[iids],
      score = mapq[iids],
      counts = counts[iids])
  } else {
    GenomicRanges::GRanges()
  }
  # DEL
  del <- GenomicAlignments::cigarRangesAlongReferenceSpace(
    cigars, ops = c("D", "N"), pos = aln_pos_start)
  origin <- Biostrings::extractAt(DNAStringSet(ref), del)
  names(del) <- ids
  del <- unlist(del, use.names = TRUE)
  iids <- as.integer(names(del))

  del <- if (length(del) > 0) {
    GenomicRanges::GRanges(
      seqnames = seqnames[iids],
      ranges = del,
      strand = strands[iids],
      originally = as.character(unlist(origin, use.names = FALSE)),
      replacement = "",
      type = "deletion",
      read_id = read_id[iids],
      score = mapq[iids],
      counts = counts[iids])
  } else {
    GenomicRanges::GRanges()
  }

  # MISMATCH - X #TODO MD tags contain mm too if X not present in CIGARS
  mm <- GenomicAlignments::cigarRangesAlongQuerySpace(cigars, ops = "X")
  repl <- Biostrings::extractAt(DNAStringSet(query_seq), mm)
  repl <- unlist(repl, use.names = FALSE)
  repl_l <- unlist(strsplit(as.character(repl[width(repl) > 1]), ""))
  repl <- as.character(repl[width(repl) == 1])
  mm <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigars, ops = "X")
  names(mm) <- ids
  mm <- unlist(mm, use.names = TRUE)
  mm_l <- mm[width(mm) > 1]
  mm <- mm[width(mm) == 1]
  mm_ln <- names(mm_l)
  mm_l <- IRanges::tile(mm_l, width = 1L)
  names(mm_l) <- mm_ln
  mm_l <- unlist(mm_l, use.names = TRUE)
  iids <- as.integer(c(names(mm), names(mm_l)))

  mm <- if (length(mm) > 0) {
    GenomicRanges::GRanges(
      seqnames = seqnames[iids],
      ranges = c(mm, mm_l),
      strand = strands[iids],
      originally = "",
      replacement = c(repl, repl_l),
      type = "mismatch",
      read_id = read_id[iids],
      score = mapq[iids],
      counts = counts[iids])
  } else {
    GenomicRanges::GRanges()
  }
  seqnames <- unique(as.character(seqnames))
  GenomeInfoDb::seqlevels(mm) <- seqnames
  GenomeInfoDb::seqlevels(del) <- seqnames
  GenomeInfoDb::seqlevels(ins) <- seqnames
  c(mm, del, ins)
}


#' Read "pair" format of EMBOSS needle into GRanges as events.
#'
#' Parse EMBOSS needle (or needleall) "pair" format into GRanges representation
#' with events of deletions, insertions and mismatches. Make sure that each file
#' corresponds to single subject (single amplicon). Assumes that bottom sequence
#' "-bsequence" corresponds to the "subject" and full sequence alignment is
#' returned.
#'
#' @param file (character) File path.
#' @param ID (character) ID of the experiment, will be used as seqnames of the
#' reutner ranges.
#' @param strand_info (character) Strand to assign.
#' @return (GRanges) Same as events.
#' @export
#'
pairToEvents <- function(file, ID = "NA", strand_info = "+") {
  file <- readLines(file)

  # header
  header <- grep("########################################", file)[2]
  file_header <- file[seq_len(header)]
  file <- file[-seq_len(header)]
  if (file_header[2] != "# Program: needleall" |
      file_header[grep("-aformat3", file_header)] != "#    -aformat3 pair") {
    stop("This function can only parse needleall 'pair' format.")
  }
  gapopen <- as.numeric(gsub("#    -gapopen ", "",
                             file_header[grep("-gapopen ", file_header)]))
  gapext <- as.numeric(gsub("#    -gapextend ", "",
                            file_header[grep("-gapextend ", file_header)]))
  seq1 <- gsub("# 1: ", "", file[grep("# 1: ", file)])
  seq2 <- gsub("# 2: ", "", file[grep("# 2: ", file)])
  score <- as.numeric(gsub("# Score: ", "", file[grep("# Score: ", file)]))

  # split into list of reads
  pos <- grep("#=======================================", file)
  starts <- pos[c(FALSE, TRUE)] + 1
  ends <- pos[c(TRUE, FALSE)] - 1
  ends <- ends[-1]
  ends <- c(ends, length(file))

  headers_pos <- seq2(starts - 1, ends - 1, by = 2)
  split_v <- vector(mode = "character", length = length(file))
  uniq_ids <- seq_len(length(headers_pos))
  split_v[unlist(headers_pos, use.names = FALSE)] <-
    rep(uniq_ids, times = lengths(headers_pos))
  alns <- split(file, factor(split_v, levels = uniq_ids))

  # for every sequence parse
  seq1_s <- substr(seq1, 1, 13)
  seq2_s <- substr(seq2, 1, 13)
  alns <- unlist(mapply(function(aln, s1, s2) {
    query <- paste0(gsub(paste0(" |[0-9]|", s1), "",
                         aln[grep(s1, aln)]), collapse = "")
    subj <- paste0(gsub(paste0(" |[0-9]|", s2), "",
                        aln[grep(s2, aln)]), collapse = "")
    c(query, subj)
  }, alns, seq1_s, seq2_s, SIMPLIFY = FALSE, USE.NAMES = FALSE),
  use.names = FALSE)

  # parsing into GRanges
  events <- getEvents(
    alns[c(TRUE, FALSE)], alns[c(FALSE, TRUE)], scores = score,
    ID = ID, ampl_shift = 1, ampl_start = 1, strand_info = strand_info)
  events$counts <- 1 # assume reads are not collapsed into unique
  events
}
