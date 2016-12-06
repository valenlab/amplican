#' @include helpers_general.R
NULL

#' Reverse complement events that have amplicons with direction 1.
#'
#' @param idRanges (data.frame) Loaded events.
#' @param configTable (data.frame) Loaded configuration file.
#' @return (data.frame) Returns input idRanges, but events for amplicons with
#' direction 1 reverse complemented.
#'
flipRanges <- function(idRanges, configTable) {

  is_dir <- as.logical(configTable$Direction)
  to_flip <- configTable[is_dir, "ID"]
  to_flip <- idRanges$seqnames %in% to_flip

  if (any(to_flip)) {
    ampl_lengths <- nchar(as.character(configTable[is_dir, "Amplicon"]))
    ampl_ids <- as.character(configTable[is_dir, "ID"])

    ids_mapping <- match(idRanges[to_flip, "seqnames"], ampl_ids)
    ampl_lengths <- ampl_lengths[ids_mapping]

    idRanges[to_flip, "originally"] <- revComp(idRanges[to_flip, "originally"])
    idRanges[to_flip, "replacement"] <- revComp(idRanges[to_flip, "replacement"])

    old_starts <- idRanges[to_flip, "start"]
    idRanges[to_flip, "start"] <- ampl_lengths - idRanges[to_flip, "end"] + 1
    idRanges[to_flip, "end"] <- ampl_lengths - old_starts + 1
  }

  return(idRanges)
}


#' Creates equal label spacing.
#' Used to calculate x label ticks.
#'
#' @param direction (numeric) 0 or 1
#' @param box (IRanges) specifies where predcited cut site is
#' @param ampl_len (numeric) Length of the amplicon
#' @param spacing (numeric) desired spacing between ticks
#' @return numeric vector
#' @importFrom IRanges start end
#'
xlabels_spacing <- function(direction, box, ampl_len, spacing) {

  if (direction != 1) {
    if (length(box) >= 1) {
      seq(-IRanges::start(box[1]) + 1, ampl_len - IRanges::start(box[1]), spacing)
    } else {
      seq(1, ampl_len, spacing)
    }
  } else {
    if (length(box) >= 1) {
      rev(seq(IRanges::end(box[1]) - ampl_len + 1, IRanges::end(box[1]), spacing))
    } else {
      rev(seq(1, ampl_len, spacing))
    }
  }
}


#' Plots amplicon sequence using ggplot2.
#'
#' @param amplicon (character) Sequence of the amplicon to plot.
#' @return (amplicon plot) ggplot2 object of amplicon plot
#' @importFrom ggplot2 ggplot aes theme_bw theme
#' scale_colour_manual element_blank unit geom_text ylim
#' @importFrom ggbio xlim
#'
amplican_plot_amplicon <- function(amplicon) {

  ampl_df <- data.frame(seq(1, nchar(amplicon)), strsplit(amplicon, "")[[1]],
                        strsplit(toupper(amplicon), "")[[1]], 1)
  names(ampl_df) <- c("position", "nucleotide", "upper", "count")

  amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                       "#009E73", "#D55E00", "#F0E442", "#0072B2")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")

  nucleotide <- position <- count <- upper <- NULL
  p <- ggplot2::ggplot(ampl_df,
                  ggplot2::aes(x = position,
                               label = nucleotide,
                               colour = upper,
                               y = count)) +
    ggplot2::geom_text(size = I(4)) +
    ggplot2::ylim(0.7, 1.1) +
    ggbio::xlim(1, dim(ampl_df)[1]) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   line = element_blank(),
                   rect = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text = element_blank(),
                   axis.title.y = element_blank()) +
    ggplot2::scale_colour_manual(drop = FALSE, values = amplicon_colors)
  return(p)
}


#' Filter Events Overlapping Primers. Message user when many of the events are
#' filtered.
#'
#' @param idRanges (data.frame) Contains events.
#' @param frPrimer (character) forward primer
#' @param rwPrimer (character) reverse primer
#' @param amplicon (character) amplicon sequence
#' @return (data.frame) filtered data frame of events
#' @importFrom stringr str_locate
#'
filterEOP <- function(idRanges, frPrimer, rwPrimer, amplicon) {

  totalSum <- sum(idRanges$count)
  # events starting before end of forward primer
  idRanges <- idRanges[!idRanges$start < frPrimer[2],]
  # events ending after start of rewerse primer
  idRanges <- idRanges[!idRanges$end > rwPrimer[1],]

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


#' Plots mismatches using ggplot2 and ggbio.
#'
#' This function plots mismatches in relation to the amplicon.
#' Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs, if it is possible to group them. They have to
#' have the same amplicon, amplicons on the reverse strand will be reverse
#' complemented to match forward strand amplicons.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
#' @param filter (boolean) Whether deletions overlapping primers should be
#' removed. By deafult it set to TRUE.
#' @return (mismatches plot) ggplot2 object of mismatches plot
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle
#' scale_colour_manual geom_bar
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse
#' element_blank unit geom_text ylab ylim
#' @importFrom ggbio tracks xlim
#' @importFrom stats na.omit aggregate
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "config.csv", package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments_events.csv", package = "amplican")
#' alignments <- read.csv(alignments_file)
#' amplican_plot_mismatches(alignments, config, c('ID_1', 'ID_3'))
#'
amplican_plot_mismatches <- function(alignments,
                                     config,
                                     id,
                                     cut_buffer = 5,
                                     xlab_spacing = 4,
                                     filter = TRUE) {

  # check for amplicons being the same
  ampliconIntegrityCheck(config, id)

  idRanges <- alignments[alignments$seqnames %in% id, ]
  idRanges <- idRanges[idRanges$type == "mismatch", ]

  if (dim(idRanges)[1] == 0) {
    return("No mismatches to plot.")
  }

  # reverse map events when amplicons have Direction 1
  idRanges <- flipRanges(idRanges, config)

  amplicon <- as.character(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)
  amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                       "#009E73", "#D55E00", "#F0E442", "#0072B2")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(config[which(config$ID == id[1]), "Direction"],
                             box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)
  box <- box + cut_buffer

  frPrimer <- as.character(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(frPrimer)))
  rwPrimer <- as.character(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 revComp(rwPrimer)))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {
    frPrimer <- c(1, 1)
  }
  if (length(rwPrimer) == 0) {
    rwPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(idRanges,
                frPrimer,
                rwPrimer,
                amplicon)
  }

  if (dim(idRanges)[1] == 0) {
    return(message("No mismatches to plot."))
  }

  # variables to NULL first to negate CRAN check
  frequency <- replacement <- start <- strand <- NULL
  freqAgr <- stats::aggregate(
    cbind(count, frequency) ~ replacement + start + strand,
    idRanges,
    sum)
  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  if (dim(freqAgrPlus)[1] == 0) {
    freqAgrPlus <- rbind(
      freqAgrPlus, data.frame(replacement = "G",
                              start = 0,
                              strand = "+",
                              count = 0,
                              frequency = 0))
  }

  if (dim(freqAgrMinus)[1] == 0) {
    freqAgrMinus <- rbind(
      freqAgrMinus, data.frame(replacement = "G",
                               start = 0,
                               strand = "+",
                               count = 0,
                               frequency = 0))
  }

  mut_fr <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = freqAgrPlus,
                      ggplot2::aes(x = as.numeric(start) + 1,
                                   y = frequency,
                                   fill = replacement),
                      stat = "identity", position = "identity") +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue")

  mut_re <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = freqAgrMinus,
                      ggplot2::aes(x = as.numeric(start) + 1,
                                   y = frequency,
                                   fill = replacement),
                      stat = "identity",
                      position = "identity") +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::scale_y_reverse()

  p <- ggbio::tracks(mut_fr,
                     amplican_plot_amplicon(amplicon),
                     mut_re,
                     heights = c(0.5, 0.03, 0.53),
                     padding = -1,
                     xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots deletions using ggplot2 and ggbio.
#'
#' This function plots deletions in relation to the amplicon.
#' Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from
#' alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
#' @param filter (boolean) Whether deletions overlapping primers should be
#' removed. By deafult it set to TRUE.
#' @return (deletions plot) ggplot2 object of deletions plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle
#' scale_colour_manual
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse
#' element_blank unit geom_text ylab ylim scale_colour_gradientn
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks geom_arch xlim
#' @importFrom stats na.omit
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "config.csv", package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments_events.csv", package = "amplican")
#' alignments <- read.csv(alignments_file)
#' amplican_plot_deletions(alignments, config, c('ID_1','ID_3'), 5)
#'
amplican_plot_deletions <- function(alignments,
                                    config,
                                    id,
                                    cut_buffer = 5,
                                    xlab_spacing = 4,
                                    filter = TRUE) {

  # check for amplicons being the same
  ampliconIntegrityCheck(config, id)

  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]

  if (dim(archRanges)[1] == 0) {
    return(print("No deletions to plot."))
  }

  archRanges <- flipRanges(archRanges, config)

  amplicon <- as.character(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  archRanges <- stats::aggregate(
    cbind(count, frequency, cut) ~ strand + start + end, archRanges, sum)
  archRanges$cut <- archRanges$cut > 0

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(config[which(config$ID == id[1]), "Direction"],
                             box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)
  box <- box + cut_buffer

  frPrimer <- as.character(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(frPrimer)))
  rwPrimer <- as.character(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 revComp(rwPrimer)))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {
    frPrimer <- c(1, 1)
  }
  if (length(rwPrimer) == 0) {
    rwPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    archRanges <-
      filterEOP(archRanges,
                frPrimer,
                rwPrimer,
                amplicon)
  }

  if (dim(archRanges)[1] == 0) {
    return(print("No deletions to plot."))
  }

  frequency <- cut <- NULL
  isCutPalette <- c("#CC79A7", "#0072B2")
  names(isCutPalette) <- c("FALSE", "TRUE")

  arch_plot_fr <- ggplot2::ggplot() +
    ggbio::geom_arch(data = archRanges[archRanges$strand == "+", ],
                     ggplot2::aes(alpha = frequency,
                                  colour = cut,
                                  size = frequency,
                                  x = start,
                                  xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::scale_colour_manual(values = isCutPalette)

  arch_plot_re <- ggplot2::ggplot() +
    ggbio::geom_arch(data = archRanges[archRanges$strand == "-", ],
                     ggplot2::aes(alpha = frequency,
                                  colour = cut,
                                  size = frequency,
                                  x = start,
                                  xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_colour_manual(values = isCutPalette)

  p <- ggbio::tracks(arch_plot_fr,
                     amplican_plot_amplicon(amplicon),
                     arch_plot_re,
                     heights = c(0.5, 0.03, 0.53),
                     padding = -1,
                     xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots insertions using ggplot2 and ggbio.
#'
#' This function plots insertions in relation to the amplicon. Top plot is for
#' the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
#' @param filter (boolean) Whether deletions overlapping primers should be
#' removed. By deafult it set to TRUE.
#' @return (insertions plot) ggplot2 object of insertions plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle
#' scale_colour_manual geom_polygon scale_fill_manual scale_x_continuous
#' geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks xlim
#' @importFrom stats na.omit aggregate
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "config.csv", package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments_events.csv", package = "amplican")
#' alignments <- read.csv(alignments_file)
#' amplican_plot_insertions(alignments, config, c('ID_1','ID_3'), 5)
#'
amplican_plot_insertions <- function(alignments,
                                     config,
                                     id,
                                     cut_buffer = 5,
                                     xlab_spacing = 4,
                                     filter = TRUE) {

  idRanges <- alignments[alignments$seqnames %in% id &
                           alignments$type == "insertion", ]

  amplicon <- as.character(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  if (dim(idRanges)[1] == 0) {
    return(print("No insertions to plot."))
  }

  idRanges <- flipRanges(idRanges, config)

  frPrimer <- as.character(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(frPrimer)))
  rwPrimer <- as.character(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 revComp(rwPrimer)))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {
    frPrimer <- c(1, 1)
  }
  if (length(rwPrimer) == 0) {
    rwPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(idRanges,
                frPrimer,
                rwPrimer,
                amplicon)
  }

  if (dim(idRanges)[1] == 0) {
    return(print("No insertions to plot."))
  }

  # reduce
  idRangesReduced <- stats::aggregate(
    frequency ~ strand + start + end, idRanges, sum)
  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+", ]

  if (dim(idRangesFr)[1] != 0) {

    idRangesFrgroup = rep(1:dim(idRangesFr)[1], each = 3)
    idRangesFrFrequency = rep(idRangesFr$frequency, each = 3)
    idRangesFrFrequency[c(TRUE, FALSE, FALSE)] <- 0
    idRangesFrX <- as.vector(rbind(idRangesFr$start,
                                   idRangesFr$start,
                                   idRangesFr$end))
    traingleFr <- data.frame(frequency = idRangesFrFrequency,
                             position = idRangesFrX,
                             group = idRangesFrgroup)
  } else {
    idRangesFrX <- c()
    traingleFr <- data.frame(frequency = c(), position = c(), group = c())
  }

  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-", ]

  if (dim(idRangesRe)[1] != 0) {

    idRangesRegroup = rep(1:dim(idRangesRe)[1], each = 3)
    idRangesReFrequency = rep(idRangesRe$frequency, each = 3)
    idRangesReFrequency[c(TRUE, FALSE, FALSE)] <- 0
    idRangesReX <- as.vector(rbind(idRangesRe$start,
                                   idRangesRe$start,
                                   idRangesRe$end))
    traingleRe <- data.frame(frequency = idRangesReFrequency,
                             position = idRangesReX,
                             group = idRangesRegroup)
  } else {
    idRangesReX <- c()
    traingleRe <- data.frame(frequency = c(), position = c(), group = c())
  }

  if (dim(idRangesRe)[1] != 0 | dim(idRangesFr)[1] != 0) {
    ampl_len <- max(c(idRangesFrX, idRangesReX, ampl_len))
  }

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(config[which(config$ID == id[1]), "Direction"],
                             box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)
  box <- box + cut_buffer

  position <- frequency <- group <- NULL
  ins_fr <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = traingleFr,
                          ggplot2::aes(x = position,
                                       y = frequency,
                                       group = group,
                                       alpha = frequency,
                                       size = frequency),
                          fill = "red") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue")

  ins_re <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = traingleRe,
                          ggplot2::aes(x = position,
                                       y = frequency,
                                       group = group,
                                       alpha = frequency,
                                       size = frequency),
                          fill = "red") +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::scale_y_reverse()

  p <- ggbio::tracks(ins_fr,
                     amplican_plot_amplicon(amplicon),
                     ins_re,
                     heights = c(0.5, 0.03, 0.53),
                     padding = -1,
                     xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots cuts using ggplot2 and ggbio.
#'
#' This function plots cuts in relation to the amplicon with distinction for
#' each ID. Top plot is meta-plot of all cuts combined with frequencies for all
#' reads in the  amplicon group. Bottom plot shows cuts with facets for each ID.
#'
#' @param alignments (GRanges object) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is possible to group them. First amplicon
#' will be used as the basis for plot.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
#' @param filter (boolean) Whether deletions overlapping primers should be
#' removed. By deafult it set to TRUE.
#' @return (cuts plot) ggplot2 object of cuts plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle
#' scale_colour_manual scale_fill_manual scale_x_continuous geom_vline
#' scale_y_reverse element_blank unit geom_text ylab ylim scale_color_manual
#' facet_grid element_text
#' @importFrom stringr str_locate
#' @importFrom ggbio geom_arch xlim
#' @importFrom stats na.omit
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "config.csv", package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments_events.csv", package = "amplican")
#' alignments <- read.csv(alignments_file)
#' amplican_plot_cuts(alignments, config, c('ID_1','ID_3'))
#'
amplican_plot_cuts <- function(alignments,
                               config,
                               id,
                               xlab_spacing = 4,
                               filter = TRUE) {

  ampliconIntegrityCheck(config, id)

  archRanges <- alignments[alignments$seqnames %in% id & alignments$type == "deletion", ]
  archRanges <- archRanges[archRanges$cut == TRUE, ]

  if (length(unique(archRanges$strand)) == 2) {
    archRanges <- archRanges[archRanges$strand == "+", ]
  }

  amplicon <- as.character(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  archRanges <- flipRanges(archRanges, config)

  if (dim(archRanges)[1] == 0) {
    return(print("No cuts to plot."))
  }

  archRanges <- stats::aggregate(
    cbind(count, frequency, cut) ~ strand + start + end + seqnames,
    archRanges,
    sum)

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(config[which(config$ID == id[1]), "Direction"],
                             box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)

  frPrimer <- as.character(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(frPrimer)))
  rwPrimer <- as.character(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 revComp(rwPrimer)))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {
    frPrimer <- c(1, 1)
  }
  if (length(rwPrimer) == 0) {
    rwPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(archRanges,
                frPrimer,
                rwPrimer,
                amplicon)
  }

  if (dim(archRanges)[1] == 0) {
    return(print("No cuts to plot."))
  }

  amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                       "#009E73", "#D55E00", "#F0E442", "#0072B2")
  amplicon_letters <- c("A", "C", "G", "T", "a", "c", "g", "t")
  ampl_df <- data.frame(seq(1, ampl_len),
                        strsplit(amplicon, "")[[1]],
                        strsplit(toupper(amplicon), "")[[1]], 1)
  names(ampl_df) <- c("position", "nucleotide", "upper", "count")
  ampl_df$colour <- amplicon_colors[match(ampl_df$upper, amplicon_letters)]

  frequency <- seqnames <- NULL
  p <- ggplot2::ggplot() +
    ggbio::geom_arch(data = archRanges,
                     ggplot2::aes(alpha = frequency,
                                  size = frequency,
                                  x = start,
                                  xend = end,
                                  colour = seqnames)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::guides(size=FALSE, alpha = FALSE) +
    ggplot2::theme(legend.position = c(1, 1),
                   legend.justification = c(1.01, 1.01),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(y = "Frequency [%]",
                  colour = "Experiments",
                  x = "Nucleotide Position Relative to PAM") +
    ggplot2::scale_x_continuous(labels = xlabels,
                                breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::annotate("text",
                      x = ampl_df$position,
                      label = ampl_df$nucleotide,
                      y = -0.1,
                      colour = ampl_df$colour)
  return(p)
}
