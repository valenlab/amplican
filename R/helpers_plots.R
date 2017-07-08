#' @include helpers_general.R
NULL


amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                     "#009E73", "#D55E00", "#F0E442", "#0072B2")
names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")
is_cut_colors <- c("#CC79A7", "#0072B2")
names(is_cut_colors) <- c("FALSE", "TRUE")


amplicon_primers <- function(config, id, amplicon) {
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
    rightPrimer <- c(nchar(amplicon), nchar(amplicon))
  }

  list(leftPrimer = leftPrimer, rightPrimer = rightPrimer, primers = primers)
}


amplican_style <- function(p) {
  p +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = unit(0, "cm")) +
    ggplot2::ylab("Frequency [%]")
}


amplican_xlim <- function(p, xlabels, box, primers) {
  p +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xlabels) +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = primers,
                        linetype = "dotdash",
                        colour = "blue")
}


ggplot_mismatches <- function(xData) {
  # variables to NULL first to negate CRAN check
  frequency <- replacement <- start <- strand <- NULL
  ggplot2::ggplot() +
    ggplot2::geom_bar(data = xData,
                      ggplot2::aes(x = as.numeric(start),
                                   y = frequency,
                                   fill = replacement),
                      stat = "identity", position = "identity") +
    ggplot2::scale_fill_manual(values = amplicon_colors)
}


ggplot_deletions <- function(xData) {
  frequency <- overlaps <- NULL
  ggplot2::ggplot() +
    ggbio::geom_arch(data = xData,
                     ggplot2::aes(alpha = frequency,
                                  colour = overlaps,
                                  size = frequency,
                                  height = frequency,
                                  x = start,
                                  xend = end)) +
    ggplot2::scale_colour_manual(values = is_cut_colors)
}


ggplot_insertions <- function(xData) {
  position <- frequency <- group <- NULL
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = xData,
                          ggplot2::aes(x = position,
                                       y = frequency,
                                       group = group,
                                       alpha = frequency,
                                       size = frequency),
                          fill = "#FF0000")
}


triangulate_ranges <- function(xRanges) {
  if (dim(xRanges)[1] != 0) {
    idRangesFrequency = rep(xRanges$frequency, each = 3)
    idRangesFrequency[c(TRUE, FALSE, FALSE)] <- 0
    data.frame(frequency = idRangesFrequency,
               position = as.vector(rbind(xRanges$start,
                                          xRanges$start,
                                          xRanges$end)),
               group = rep(1:dim(xRanges)[1], each = 3))
  } else {
    data.frame(frequency = c(), position = c(), group = c())
  }
}


group_to_selection <- function(alnmt, config, group, selection) {
  alnmt <- alnmt[alnmt$seqnames %in%
                   unique(config$ID[config[[group]] %in% selection]), ]
  alnmt$ID <- alnmt$seqnames
  alnmt$seqnames <- alnmt[[group]]
  return(alnmt)
}


mock_mm_df <- function(ampl_max, ampl_min = 0) {
  how_many <- length(ampl_min:ampl_max)
  data.frame(replacement = rep("G", how_many),
             start = ampl_min:ampl_max,
             strand = rep("+", how_many),
             frequency = rep(0, how_many))
}


return_metaplot <- function(freqAgr, plot_fr, plot_re) {
  if (any(freqAgr$strand == "+") & any(freqAgr$strand == "-")) {
    return(ggbio::tracks(plot_fr,
                         plot_re +
                           ggplot2::scale_y_reverse(
                             limits = c(max(freqAgr$frequency, na.rm = TRUE),
                                        0)),
                         heights = c(0.5, 0.5),
                         padding = -1))
  } else if (all(freqAgr$strand == "+")) {
    return(plot_fr)
  } else {
    return(plot_re)
  }
}


annotate_with_amplicon <- function(p, amplicon, from, to) {
  ampl_len <- nchar(amplicon)
  amplicon <- strsplit(amplicon, "")[[1]]
  p +
    ggplot2::annotate(
      "text",
      x = seq(from, to),
      label = amplicon,
      y = 0,
      colour = amplicon_colors[match(toupper(amplicon),
                                     names(amplicon_colors))])
}


return_plot <- function(freqAgr, amplicon, from, to, plot_fr, plot_re) {
  if (any(freqAgr$strand == "+") & any(freqAgr$strand == "-")) {
    return(ggbio::tracks(plot_fr,
                         plot_amplicon(amplicon, from, to),
                         plot_re,
                         heights = c(0.5, 0.06, 0.5),
                         padding = -1,
                         xlim = from:to,
                         xlab = "Relative Nucleotide Position"))
  } else if (all(freqAgr$strand == "+")) {
    return(annotate_with_amplicon(plot_fr, amplicon))
  } else {
    return(annotate_with_amplicon(plot_re, amplicon))
  }
}


#' MetaPlots mismatches using ggplot2 and ggbio.
#'
#' Plots mismatches in relation to the amplicons for given
#' selection vector that groups values by given config group. All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{map_to_relative}}.
#' For zero position on new coordinates is the most left UPPER case letter of
#' the respective amplicon. This function filters out all alignment events
#' that have amplicons without UPPER case defined.
#' Top plot is for the forward reads and bottom plot is for reverse reads.
#'
#' @param alnmt (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param group (string) Name of the column from the config file to use for
#' grouping. Events are subselected based on this column and values from
#' selection.
#' @param selection (string or vector of strings) Values from config column
#' specified in group argument.
#' @return (mismatches metaplot) ggplot2 object of mismatches metaplot
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' metaplot_mismatches(alignments, config, "Group", "Betty")
#'
metaplot_mismatches <- function(alnmt, config, group,
                                selection) {
  alnmt <- alnmt[alnmt$type == "mismatch",]
  if (length(alnmt) == 0) {
    return("No mismatches to plot.")
  }

  alnmt[[group]] <- config[[group]][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)
  if (dim(alnmt)[1] == 0) {
    return("No mismatches to plot.")
  }
  if (!any(colnames(alnmt) == "frequency")) {
    alnmt$frequency <- alnmt$counts /
      config$Reads_noPD[match(alnmt$ID, config$ID)]
  }

  freqAgr <- stats::aggregate(frequency ~ replacement + start + strand + ID,
                              alnmt, sum)
  freqAgr <- stats::aggregate(frequency ~ replacement + start + strand,
                              freqAgr, mean)

  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  # sometimes barplot gets confused so we add mock data
  mock <- mock_mm_df(max(freqAgr$start, na.rm = TRUE),
                     min(freqAgr$start, na.rm = TRUE))
  freqAgrPlus <- rbind(freqAgrPlus, mock)
  freqAgrMinus <- rbind(freqAgrMinus, mock)

  mut_fr <- ggplot_mismatches(freqAgrPlus) +
    ggplot2::ylim(0, max(freqAgr$frequency, na.rm = TRUE))
  mut_fr <- amplican_style(mut_fr)

  mut_re <- ggplot_mismatches(freqAgrMinus)
  mut_re <- amplican_style(mut_re)

  return_metaplot(freqAgr, mut_fr, mut_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' MetaPlots deletions using ggplot2 and ggbio.
#'
#' This function plots deletions in relation to the amplicons for given
#' selection vector that groups values by given config group. All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{map_to_relative}}.
#' Top plot is for the forward reads and bottom plot is for reverse reads.
#'
#' @param alnmt (data.frame) Loaded alignment information from
#' events_filtered_shifted_normalized.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param group (string) Name of the column from the config file to use for
#' grouping. Events are subselected based on this column and values from
#' selection.
#' @param selection (string or vector of strings) Values from config column
#' specified in group argument.
#' @return (deletions metaplot) ggplot2 object of deletions metaplot
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' metaplot_deletions(alignments, config, "Group", "Tom")
#'
metaplot_deletions <- function(alnmt, config, group,
                              selection) {
  alnmt <- alnmt[alnmt$type == "deletion",]
  alnmt[[group]] <- config[[group]][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)

  if (dim(alnmt)[1] == 0) {
    return("No deletions to plot.")
  }
  if (!any(colnames(alnmt) == "frequency")) {
    alnmt$frequency <- alnmt$counts /
      config$Reads_noPD[match(alnmt$ID, config$ID)]
  }

  archRanges <- stats::aggregate(
    cbind(frequency, overlaps) ~ strand + start + end + ID, alnmt, sum)
  archRanges$overlaps <- archRanges$overlaps > 0
  archRanges <- stats::aggregate(
    frequency ~ overlaps + strand + start + end, archRanges, mean)
  archRanges$overlaps <- archRanges$overlaps > 0

  arch_plot_fr <- ggplot_deletions(archRanges[archRanges$strand == "+", ])
  arch_plot_fr <- amplican_style(arch_plot_fr) +
    ggplot2::ylim(0, max(archRanges$frequency, na.rm = TRUE))

  arch_plot_re <- ggplot_deletions(archRanges[archRanges$strand == "-", ])
  arch_plot_re <- amplican_style(arch_plot_re)

  return_metaplot(archRanges, arch_plot_fr, arch_plot_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' MetaPlots insertions using ggplot2 and ggbio.
#'
#' This function plots insertions in relation to the amplicons for given
#' selection vector that groups values by given config group. All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{map_to_relative}}.
#' Top plot is for the forward reads and bottom plot is for reverse reads.
#'
#' @param alnmt (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param group (string) Name of the column from the config file to use for
#' grouping. Events are subselected based on this column and values from
#' selection.
#' @param selection (string or vector of strings) Values from config column
#' specified in group argument.
#' @return (insertions metaplot) ggplot2 object of insertions metaplot
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' metaplot_insertions(alignments, config, "Group", "Betty")
#'
metaplot_insertions <- function(alnmt, config, group,
                                selection) {
  alnmt <- alnmt[alnmt$type == "insertion",]
  alnmt[,group] <- config[,group][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)

  if (dim(alnmt)[1] == 0) {
    return("No insertions to plot.")
  }

  if (!any(colnames(alnmt) == "frequency")) {
    alnmt$frequency <- alnmt$counts /
      config$Reads_noPD[match(alnmt$ID, config$ID)]
  }

  # reduce
  idRangesReduced <- stats::aggregate(
    frequency ~ strand + start + end + ID, alnmt, sum)
  idRangesReduced <- stats::aggregate(
    frequency ~ strand + start + end, idRangesReduced, mean)
  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+", ]
  triangleFr <- triangulate_ranges(idRangesFr)
  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-", ]
  triangleRe <- triangulate_ranges(idRangesRe)

  ins_fr <- ggplot_insertions(triangleFr) +
    ggplot2::ylim(0, max(idRangesReduced$frequency, na.rm = TRUE))
  ins_fr <- amplican_style(ins_fr)

  ins_re <- ggplot_insertions(triangleRe)
  ins_re <- amplican_style(ins_re)

  return_metaplot(idRangesReduced, ins_fr, ins_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' Plots amplicon sequence using ggplot2.
#'
#' @param amplicon (character) Sequence of the amplicon to plot.
#' @param from (number) Minimum on x axis
#' @param to (number) Maximum on x axis
#' @return (amplicon plot) ggplot2 object of amplicon plot
#'
plot_amplicon <- function(amplicon, from, to) {

  ampl_df <- data.frame(position = seq(from, to),
                        nucleotide = strsplit(amplicon, "")[[1]],
                        upper = strsplit(toupper(amplicon), "")[[1]],
                        counts = 1)

  nucleotide <- position <- counts <- upper <- NULL
  p <- ggplot2::ggplot(ampl_df,
                       ggplot2::aes(x = position,
                                    label = nucleotide,
                                    colour = upper,
                                    y = counts)) +
    ggplot2::geom_text(size = I(4)) +
    ggplot2::ylim(0.7, 1.1) +
    ggbio::xlim(from, to) +
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
#' This function plots mismatches in relation to the amplicon, function assumes
#' your reads are relative to the respective amplicon sequences prediced cut
#' sites.
#' Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs, if it is possible to group them. They have to
#' have the same amplicon, amplicons on the reverse strand will be reverse
#' complemented to match forward strand amplicons.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
#' @return (mismatches plot) ggplot2 object of mismatches plot
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_mismatches(alignments, config, c('ID_1', 'ID_3'))
#'
plot_mismatches <- function(alignments,
                            config,
                            id,
                            cut_buffer = 5,
                            xlab_spacing = 4) {

  idRanges <- alignments[alignments$seqnames %in% id, ]
  idRanges <- idRanges[idRanges$type == "mismatch", ]

  if (dim(idRanges)[1] == 0) {
    return("No mismatches to plot.")
  }

  if (!any(colnames(idRanges) == "frequency")) {
    idRanges$frequency <- idRanges$counts /
      config$Reads_noPD[match(idRanges$seqnames, config$ID)]
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)
  box <- upperGroups(amplicon)

  from <- if (length(box) >= 1) -IRanges::start(box[1]) + 1 else 1
  to <- if (length(box) >= 1) ampl_len - IRanges::start(box[1]) else ampl_len
  xlabels <- seq(from, to, xlab_spacing)
  pr <- amplicon_primers(config, id, amplicon)

  if (length(box) >= 1) {
    box_shift <- IRanges::start(box)[1]
    box <- IRanges::shift(box, -1 * box_shift)
    pr$primers <- pr$primers - box_shift
  }
  box <- box + cut_buffer

  freqAgr <- stats::aggregate(
    frequency ~ replacement + start + strand + seqnames, idRanges, sum)
  freqAgr <- stats::aggregate(
    frequency ~ replacement + start + strand, freqAgr, mean)
  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  # sometimes barplot gets confused so we add mock data
  mock <- mock_mm_df(ampl_len)
  freqAgrPlus <- rbind(freqAgrPlus, mock)
  freqAgrMinus <- rbind(freqAgrMinus, mock)

  mut_fr <- ggplot_mismatches(freqAgrPlus) +
    ggplot2::ylim(0, max(freqAgr$frequency, na.rm = TRUE))
  mut_fr <- amplican_xlim(mut_fr, xlabels, box, pr$primers)
  mut_fr <- amplican_style(mut_fr)

  mut_re <- ggplot_mismatches(freqAgrMinus)
  mut_re <- amplican_xlim(mut_re, xlabels, box, pr$primers)
  mut_re <- amplican_style(mut_re) +
    ggplot2::scale_y_reverse(
      limits = c(max(freqAgr$frequency, na.rm = TRUE), 0))

  return_plot(freqAgr, amplicon, from, to, mut_fr, mut_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' Plots deletions using ggplot2 and ggbio.
#'
#' This function plots deletions in relation to the amplicon, assumes events
#' are relative to the expected cut site.
#' Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
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
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_deletions(alignments, config, c('ID_1','ID_3'), 5)
#'
plot_deletions <- function(alignments,
                           config,
                           id,
                           cut_buffer = 5,
                           xlab_spacing = 4) {

  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]

  if (dim(archRanges)[1] == 0) {
    return("No deletions to plot.")
  }

  if (!any(colnames(archRanges) == "frequency")) {
    archRanges$frequency <- archRanges$counts /
      config$Reads_noPD[match(archRanges$seqnames, config$ID)]
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)
  box <- upperGroups(amplicon)

  from <- if (length(box) >= 1) -IRanges::start(box[1]) + 1 else 1
  to <- if (length(box) >= 1) ampl_len - IRanges::start(box[1]) else ampl_len
  xlabels <- seq(from, to, xlab_spacing)
  pr <- amplicon_primers(config, id, amplicon)

  if (length(box) >= 1) {
    box_shift <- IRanges::start(box)[1]
    box <- IRanges::shift(box, -1 * box_shift)
    pr$primers <- pr$primers - box_shift
  }
  box <- box + cut_buffer

  archRanges <- stats::aggregate(
    cbind(frequency, overlaps) ~ strand + start + end + seqnames, archRanges,
    sum)
  archRanges$overlaps <- archRanges$overlaps > 0
  archRanges <- stats::aggregate(
    cbind(frequency, overlaps) ~ strand + start + end, archRanges, mean)
  archRanges$overlaps <- archRanges$overlaps > 0

  if (dim(archRanges)[1] == 0) {
    return("No deletions to plot.")
  }

  arch_plot_fr <- ggplot_deletions(archRanges[archRanges$strand == "+", ]) +
    ggplot2::ylim(0, max(archRanges$frequency, na.rm = TRUE))
  arch_plot_fr <- amplican_xlim(arch_plot_fr, xlabels, box, pr$primers)
  arch_plot_fr <- amplican_style(arch_plot_fr)

  arch_plot_re <- ggplot_deletions(archRanges[archRanges$strand == "-", ])
  arch_plot_re <- amplican_xlim(arch_plot_re, xlabels, box, pr$primers)
  arch_plot_re <- amplican_style(arch_plot_re) +
    ggplot2::scale_y_reverse(
      limits = c(max(archRanges$frequency, na.rm = TRUE), 0))

  return_plot(archRanges, amplicon, from, to, arch_plot_fr, arch_plot_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' Plots insertions using ggplot2 and ggbio.
#'
#' This function plots insertions in relation to the amplicon. Top plot is for
#' the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
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
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_insertions(alignments, config, c('ID_1','ID_3'), 5)
#'
plot_insertions <- function(alignments,
                            config,
                            id,
                            cut_buffer = 5,
                            xlab_spacing = 4) {

  idRanges <- alignments[alignments$seqnames %in% id &
                           alignments$type == "insertion", ]

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  if (dim(idRanges)[1] == 0) {
    return("No insertions to plot.")
  }

  if (!any(colnames(idRanges) == "frequency")) {
    idRanges$frequency <- idRanges$counts /
      config$Reads_noPD[match(idRanges$seqnames, config$ID)]
  }

  # reduce
  idRangesReduced <- stats::aggregate(
    frequency ~ strand + start + end + seqnames, idRanges, sum)
  idRangesReduced <- stats::aggregate(
    frequency ~ strand + start + end, idRangesReduced, mean)

  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+", ]
  triangleFr <- triangulate_ranges(idRangesFr)
  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-", ]
  triangleRe <- triangulate_ranges(idRangesRe)

  if (dim(idRangesRe)[1] != 0 | dim(idRangesFr)[1] != 0) {
    ampl_len <- max(c(max(triangleFr$position, na.rm = TRUE),
                      max(triangleRe$position, na.rm = TRUE),
                      ampl_len), na.rm = TRUE)
  }

  box <- upperGroups(amplicon)
  from <- if (length(box) >= 1) -IRanges::start(box[1]) + 1 else 1
  to <- if (length(box) >= 1) ampl_len - IRanges::start(box[1]) else ampl_len
  xlabels <- seq(from, to, xlab_spacing)
  pr <- amplicon_primers(config, id, amplicon)

  if (length(box) >= 1) {
    box_shift <- IRanges::start(box)[1]
    box <- IRanges::shift(box, -1 * box_shift)
    pr$primers <- pr$primers - box_shift
  }
  box <- box + cut_buffer

  ins_fr <- ggplot_insertions(triangleFr) +
    ggplot2::ylim(0, max(idRangesReduced$frequency, na.rm = TRUE))
  ins_fr <- amplican_xlim(ins_fr, xlabels, box, pr$primers)
  ins_fr <- amplican_style(ins_fr)

  ins_re <- ggplot_insertions(triangleRe)
  ins_re <- amplican_xlim(ins_re, xlabels, box, pr$primers)
  ins_re <- amplican_style(ins_re) +
    ggplot2::scale_y_reverse(
      limits = c(max(idRangesReduced$frequency, na.rm = TRUE), 0))

  return_plot(idRangesReduced, amplicon, from, to, ins_fr, ins_re) +
    ggplot2::xlab("Relative Nucleotide Position")
}


#' Plots cuts using ggplot2 and ggbio.
#'
#' This function plots cuts in relation to the amplicon with distinction for
#' each ID. Top plot is meta-plot of all cuts combined with frequencies for all
#' reads in the  amplicon group. Bottom plot shows cuts with facets for each ID.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is possible to group them. First amplicon
#' will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param xlab_spacing (numeric) Spacing of the x axis labels. Default is 4.
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
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_cuts(alignments, config, c('ID_1','ID_3'))
#'
plot_cuts <- function(alignments,
                      config,
                      id,
                      cut_buffer = 5,
                      xlab_spacing = 4) {

  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]
  archRanges <- archRanges[archRanges$overlaps, ]

  if (length(unique(archRanges$strand)) == 2) {
    archRanges <- archRanges[archRanges$strand == "+", ]
  }

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  if (dim(archRanges)[1] == 0) {
    return("No cuts to plot.")
  }

  if (!any(colnames(archRanges) == "frequency")) {
    archRanges$frequency <- archRanges$counts /
      config$Reads_noPD[match(archRanges$seqnames, config$ID)]
  }

  archRanges <- stats::aggregate(
    frequency ~ strand + start + end + seqnames,
    archRanges, sum)

  box <- upperGroups(amplicon)
  from <- if (length(box) >= 1) -IRanges::start(box[1]) + 1 else 1
  to <- if (length(box) >= 1) ampl_len - IRanges::start(box[1]) else ampl_len
  xlabels <- seq(from, to, xlab_spacing)
  pr <- amplicon_primers(config, id, amplicon)

  if (length(box) >= 1) {
    box_shift <- IRanges::start(box)[1]
    box <- IRanges::shift(box, -1 * box_shift)
    pr$primers <- pr$primers - box_shift
  }
  box <- box + cut_buffer

  if (dim(archRanges)[1] == 0) {
    return("No cuts to plot.")
  }

  amplicon <- strsplit(amplicon, "")[[1]]
  frequency <- seqnames <- NULL
  p <- ggplot2::ggplot() +
    ggbio::geom_arch(data = archRanges,
                     ggplot2::aes(alpha = frequency,
                                  size = frequency,
                                  height = frequency,
                                  x = start,
                                  xend = end,
                                  colour = seqnames)) +
    ggbio::xlim(from, to) +
    ggplot2::theme_bw() +
    ggplot2::guides(size = FALSE, alpha = FALSE) +
    ggplot2::theme(legend.position = c(1, 1),
                   legend.justification = c(1.01, 1.01)) +
    ggplot2::labs(y = "Frequency [%]",
                  colour = "Experiments",
                  x = "Nucleotide Position Relative to PAM") +
    ggplot2::scale_x_continuous(labels = xlabels,
                                breaks = xlabels) +
    ggplot2::geom_vline(xintercept = pr$primers,
                        linetype = "dotdash",
                        colour = "blue") +
    ggplot2::geom_vline(xintercept = c(IRanges::start(box), IRanges::end(box)),
                        linetype = "longdash",
                        colour = "black") +
    ggplot2::annotate("text",
                      x = seq(from, to),
                      label = amplicon,
                      y = 0,
                      colour = amplicon_colors[match(toupper(amplicon),
                                                     names(amplicon_colors))])
  return(p)
}


#' Plots heterogeneity of the reads using ggplot2 and ggbio.
#'
#' This function creates stacked barplot explaining
#' reads heterogeneity. It groups reads
#' by user defined levels and measures how unique are reads in
#' this level. Uniqueness of reads is simplified to the bins and
#' colored according to the color gradient. Default color black
#' indicates very high heterogeneity of the reads. The more yellow
#' the more similar are reads and less heterogenic.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param level (string) Name of the column from config
#' file specifying levels to group by.
#' @param colors (html colors vector) Two colours for gradient,
#' eg. c('#000000', '#F0E442').
#' @param bins (numeric vector) Numeric vector from 0 to 100 specyfying bins eg.
#' c(0, 5, seq(10, 100, 10)).
#' @return (heterogeneity plot) ggplot2 object of heterogeneity plot
#' @importFrom ggplot2 ggplot aes_string theme_bw theme geom_label ggtitle
#' scale_colour_manual scale_fill_manual scale_x_continuous geom_vline xlab
#' coord_flip
#' scale_y_reverse element_blank unit geom_text ylab ylim scale_color_manual
#' facet_grid element_text
#' @importFrom ggbio geom_arch xlim
#' @importFrom stats na.omit formula
#' @importFrom grDevices colorRampPalette
#' @export
#' @family specialized plots
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_heterogeneity(alignments, config)
#'
plot_heterogeneity <- function(alignments,
                               config,
                               level = "ID",
                               colors = c('#000000', '#F0E442'),
                               bins = c(0, 5, seq(10, 100, 10))) {

  alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)
  uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id),
                                c('seqnames', 'read_id', 'counts')]
  if (level != "ID") {
    by = level
    howManyTimes <- aggregate(read_id ~ seqnames,
                              data = uniqueReadsByID, length)
    howManyTimes <- howManyTimes[
      order(match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames))),]

    uniqueReadsByID[[by]] <- rep(
      config[order(match(howManyTimes$seqnames, unique(config$ID))), by],
      times = howManyTimes$read_id)
  } else {
    by = "seqnames"
  }
  c_ord <- order(uniqueReadsByID[[by]], uniqueReadsByID$counts, decreasing = T)
  uniqueReadsByID <- uniqueReadsByID[c_ord, ]

  cumsum_list <- tapply(uniqueReadsByID$counts,
                        uniqueReadsByID[[by]], FUN = cumsum)
  cumsum_list <- cumsum_list[unique(uniqueReadsByID[[by]])]
  uniqueReadsByID$cumsum <- unlist(cumsum_list)

  howManyTimes <- table(uniqueReadsByID[[by]])
  howManyTimes <- howManyTimes[
    order(match(names(howManyTimes), unique(uniqueReadsByID[[by]])))]
  dimHMT <- if (is.null(dim(howManyTimes)[1])) {
    length(howManyTimes)
  } else {
    dim(howManyTimes)[1]
  }
  uniqueReadsByID$read_number <-
    as.vector(unlist(seq2(from = rep(1, dimHMT), to = howManyTimes, by = 1)))

  ids_with_reads <- tapply(
    uniqueReadsByID$cumsum, uniqueReadsByID[[by]], FUN = max)
  ids_with_reads <- ids_with_reads[
    order(match(names(ids_with_reads), unique(uniqueReadsByID[[by]])))]
  toDivide <- rep(ids_with_reads, times = howManyTimes)
  uniqueReadsByID$read_share_percentage_normal <-
    uniqueReadsByID$counts * 100 / toDivide

  # divide into bins for colour
  uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal,
                              bins)
  # reduce number of reads in 0-5 group - faster plots without artifacts
  uniqueReadsByID <- aggregate(
    stats::formula(paste0("read_share_percentage_normal ~ ", by, " + bins")),
    data = uniqueReadsByID, sum)
  colorPalette <- colorRampPalette(colors)(length(levels(uniqueReadsByID$bins)))
  names(colorPalette) <- levels(uniqueReadsByID$bins)
  ggplot(data = uniqueReadsByID,
         aes_string(x = paste0("as.factor(", by, ")"),
                    y = "read_share_percentage_normal",
                    fill = "bins",
                    order =  paste0("as.factor(", by, ")"))) +
    ggplot2::geom_bar(position='stack', stat='identity') +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = 'bold'),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.title = element_blank()) +
    ylab('Unique reads percentage of shares') +
    xlab(level) +
    scale_fill_manual(values = colorPalette) +
    coord_flip()
}
