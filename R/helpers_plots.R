#' @include helpers_general.R
NULL

# alignments <- amplican::map_to_relative(config, alignments)
# alignments <- alignments[seqnames(alignments) %in% config$ID[config$Barcode == "B1"]]
# alignment <- alignments[alignments$type == "deletion"]
# alignment$test <- seqnames(alignment)
# seqlevels(alignment) <- c("B1", seqlevels(alignment))
# seqnames(alignment)[seq_along(alignment)] <- rep("B1", length(alignment))
#
# ggplot() +
#   geom_arch(data = alignment,
#             ggplot2::aes(alpha = frequency,
#                          colour = cut,
#                          size = frequency,
#                          height = frequency,
#                          x = start,
#                          xend = end))

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
    ggplot2::theme(legend.position = "none",
                   legend.spacing = unit(0, "cm"),
                   line = element_blank(),
                   rect = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text = element_blank(),
                   axis.title.y = element_blank()) +
    ggplot2::scale_colour_manual(drop = FALSE, values = amplicon_colors)
  return(p)
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

  idRanges <- alignments[alignments$seqnames %in% id, ]
  idRanges <- idRanges[idRanges$type == "mismatch", ]

  if (dim(idRanges)[1] == 0) {
    return("No mismatches to plot.")
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)
  amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                       "#009E73", "#D55E00", "#F0E442", "#0072B2")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)
  box <- box + cut_buffer

  leftPrimer <- get_left_primer(config, id)
  leftPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                   toupper(leftPrimer)))
  rightPrimer <- get_right_primer(config, id)
  rightPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                    toupper(rightPrimer)))
  primers <- c(leftPrimer, rightPrimer)
  if (length(leftPrimer) == 0) {
    leftPrimer <- c(1, 1)
  }
  if (length(rightPrimer) == 0) {
    rightPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(idRanges,
                leftPrimer,
                rightPrimer,
                amplicon)
  }

  if (dim(idRanges)[1] == 0) {
    return("No mismatches to plot.")
  }

  # variables to NULL first to negate CRAN check
  frequency <- replacement <- start <- strand <- NULL
  freqAgr <- stats::aggregate(
    cbind(count, frequency) ~ replacement + start + strand,
    idRanges,
    sum)
  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  # sometimes barplot gets confused so we add mock data
  mock <- data.frame(replacement = rep("G", ampl_len),
                     start = 1:ampl_len,
                     strand = rep("+", ampl_len),
                     count = rep(0, ampl_len),
                     frequency = rep(0, ampl_len))
  freqAgrPlus <- rbind(freqAgrPlus, mock)
  freqAgrMinus <- rbind(freqAgrMinus, mock)

  mut_fr <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = freqAgrPlus,
                      ggplot2::aes(x = as.numeric(start),
                                   y = frequency,
                                   fill = replacement),
                      stat = "identity", position = "identity") +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                      ggplot2::aes(x = as.numeric(start),
                                   y = frequency,
                                   fill = replacement),
                      stat = "identity",
                      position = "identity") +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                     xlab = "Relative Nucleotide Position")
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

  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]

  if (dim(archRanges)[1] == 0) {
    return("No deletions to plot.")
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  archRanges <- stats::aggregate(
    cbind(count, frequency, cut) ~ strand + start + end, archRanges, sum)
  archRanges$cut <- archRanges$cut > 0

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)
  box <- box + cut_buffer

  leftPrimer <- get_left_primer(config, id)
  leftPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                   toupper(leftPrimer)))
  rightPrimer <- get_right_primer(config, id)
  rightPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                    toupper(rightPrimer)))
  primers <- c(leftPrimer, rightPrimer)
  if (length(leftPrimer) == 0) {
    leftPrimer <- c(1, 1)
  }
  if (length(rightPrimer) == 0) {
    rightPrimer <- c(ampl_len, ampl_len)
  }
  if (length(leftPrimer) == 0) {
    leftPrimer <- c(1, 1)
  }
  if (length(rightPrimer) == 0) {
    rightPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    archRanges <-
      filterEOP(archRanges,
                leftPrimer,
                rightPrimer,
                amplicon)
  }

  if (dim(archRanges)[1] == 0) {
    return("No deletions to plot.")
  }

  frequency <- cut <- NULL
  isCutPalette <- c("#CC79A7", "#0072B2")
  names(isCutPalette) <- c("FALSE", "TRUE")

  arch_plot_fr <- ggplot2::ggplot() +
    ggbio::geom_arch(data = archRanges[archRanges$strand == "+", ],
                     ggplot2::aes(alpha = frequency,
                                  colour = cut,
                                  size = frequency,
                                  height = frequency,
                                  x = start,
                                  xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                                  height = frequency,
                                  x = start,
                                  xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                     xlab = "Relative Nucleotide Position")
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

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  if (dim(idRanges)[1] == 0) {
    return("No insertions to plot.")
  }

  leftPrimer <- get_left_primer(config, id)
  leftPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                   toupper(leftPrimer)))
  rightPrimer <- get_right_primer(config, id)
  rightPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                    toupper(rightPrimer)))
  primers <- c(leftPrimer, rightPrimer)
  if (length(leftPrimer) == 0) {
    leftPrimer <- c(1, 1)
  }
  if (length(rightPrimer) == 0) {
    rightPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(idRanges,
                leftPrimer,
                rightPrimer,
                amplicon)
  }

  if (dim(idRanges)[1] == 0) {
    return("No insertions to plot.")
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
  xlabels <- xlabels_spacing(box, ampl_len, xlab_spacing)
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
                          fill = "#FF0000") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                          fill = "#FF0000") +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
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
                     xlab = "Relative Nucleotide Position")
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

  archRanges <- alignments[alignments$seqnames %in% id & alignments$type == "deletion", ]
  archRanges <- archRanges[archRanges$cut == TRUE, ]

  if (length(unique(archRanges$strand)) == 2) {
    archRanges <- archRanges[archRanges$strand == "+", ]
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  if (dim(archRanges)[1] == 0) {
    return("No cuts to plot.")
  }

  archRanges <- stats::aggregate(
    cbind(count, frequency, cut) ~ strand + start + end + seqnames,
    archRanges,
    sum)

  box <- upperGroups(amplicon)
  xlabels <- xlabels_spacing(box, ampl_len, xlab_spacing)
  xbreaks <- seq(1, ampl_len, xlab_spacing)

  leftPrimer <- get_left_primer(config, id)
  leftPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                   toupper(leftPrimer)))
  rightPrimer <- get_right_primer(config, id)
  rightPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                    toupper(rightPrimer)))
  primers <- c(leftPrimer, rightPrimer)
  if (length(leftPrimer) == 0) {
    leftPrimer <- c(1, 1)
  }
  if (length(rightPrimer) == 0) {
    rightPrimer <- c(ampl_len, ampl_len)
  }

  if (filter) {
    idRanges <-
      filterEOP(archRanges,
                leftPrimer,
                rightPrimer,
                amplicon)
  }

  if (dim(archRanges)[1] == 0) {
    return("No cuts to plot.")
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
                                  height = frequency,
                                  x = start,
                                  xend = end,
                                  colour = seqnames)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::guides(size=FALSE, alpha = FALSE) +
    ggplot2::theme(legend.position = c(1, 1),
                   legend.justification = c(1.01, 1.01)) +
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
                      y = 0,
                      colour = ampl_df$colour)
  return(p)
}
