seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))


#' Reverse and complement given string or list of strings
#'
#' @param x (string or vector of strings)
#' @return (string or vector of strings) reverse complemented input
#' @importFrom Biostrings DNAStringSet reverseComplement
#'
revComp <- function(x) {
  return(
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  )
}


#' Creates equal label spacing.
#' Used to calculate x label ticks.
#'
#' @param box (IRanges) specifies where predcited cut site is
#' @param ampl_len (numeric) Length of the amplicon
#' @param spacing (numeric) desired spacing between ticks
#' @return numeric vector
#' @importFrom IRanges start end
#'
xlabels_spacing <- function(box, ampl_len, spacing) {

  if (length(box) >= 1) {
    seq(-IRanges::start(box[1]) + 1, ampl_len - IRanges::start(box[1]), spacing)
  } else {
    seq(1, ampl_len, spacing)
  }
}


#' Filter Events Overlapping Primers. Message user when many of the events are
#' filtered.
#'
#' @param idRanges (data.frame) Contains events.
#' @param frPrimer (character) forward primer location
#' @param rvPrimer (character) reverse primer location
#' @return (data.frame) filtered data frame of events
#' @importFrom stringr str_locate
#'
filterEOP <- function(idRanges, frPrimer, rvPrimer) {

  totalSum <- sum(idRanges$count)
  # events starting before end of forward primer
  idRanges <- idRanges[!idRanges$start < frPrimer[2],]
  # events ending after start of reverse primer
  idRanges <- idRanges[!idRanges$end > rvPrimer[1],]

  totalFiltered <- totalSum - sum(idRanges$count)

  # message user when more than 50% of the reads filtered
  if (totalFiltered/totalSum >= 0.5) {
    warning(paste0(toString(totalFiltered),
                   " events out of ",
                   toString(totalSum),
                   " were filtered due to overlaping primers"))
  }

  return(idRanges)
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


#' Map events to their respective relative coordinates specified with
#' UPPER case.
#'
#' Translate coordinates of GRanges events so that they can be relative to the
#' amplicon. As point zero we assume first left sided UPPER case letter in the
#' amplicon.
#'
#' @param events (data.frame) List of events to map to the relative coordinates.
#' @param config (data.frame) config table
#' @return (GRanges) Same as events, but the coordinates are relative to the
#' @importFrom GenomicRanges GRanges shift
#' @export
#' @examples
#' # example config
#' config <- read.csv(system.file("extdata", "config.csv",
#'                    package = "amplican"))
#' # example events
#' events <- read.csv(system.file("extdata", "results",
#'                    "alignments_events.csv", package = "amplican"))
#' # make events relative to the UPPER case
#' map_to_relative(events, config)
#'
map_to_relative <- function(events, config) {
  events <- GRanges(events)
  no_upper <- FALSE

  for (id in unique(seqnames(events))) {
    amplicon <- get_amplicon(config, id)
    zero_point <- upperGroups(amplicon)
    if (length(zero_point) == 0) {
      no_upper <- TRUE
      events <- events[seqnames(events) != id]
      next()
    }
    events[seqnames(events) == id] <- shift(events[seqnames(events) == id],
                                            shift = -1 * start(zero_point)[1])
  }

  if (no_upper) {
    warning(paste0("Events for amplicons without UPPER case are filtered."))
  }

  return(events)
}
