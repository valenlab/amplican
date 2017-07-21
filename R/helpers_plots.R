#' @include helpers_general.R
NULL

# orange #E15D44
# red #BC243C
# blue #98B4D4
# magenta #D65076
# green #009B77
# grey #DFCFBE
amplicon_colors <- c("#009E73", "#D55E00", "#F0E442", "#0072B2",
                     "#009E73", "#D55E00", "#F0E442", "#0072B2", "#FFFFFF")
amplicon_colors <- c("#34A853", "#EA4335","#FBBC05", "#4285F4",
                     "#34A853", "#EA4335","#FBBC05", "#4285F4", "#FFFFFF")
names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t", "-")
codon_colors <- c("#E15D44", "#98B4D4", "#D65076", "#BC243C", "#009B77",
                  "#D65076", "#BC243C", "#E15D44", "#D65076", "#009B77",
                  "#009B77", "#98B4D4", "#009B77", "#009B77", "#009B77",
                  "#E15D44", "#E15D44", "#009B77", "#009B77", "#009B77",
                  rep("#DFCFBE", 10), "#FFFFFF")
names(codon_colors) <- c(Biostrings::AA_ALPHABET, "")
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
    ggplot2::theme(panel.spacing = grid::unit(0, "cm"),
                   legend.position = "none",
                   legend.spacing = grid::unit(0, "cm")) +
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
#' metaplot_deletions(alignments, config, "Group", "Betty")
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
                   legend.spacing = grid::unit(0, "cm"),
                   line = ggplot2::element_blank(),
                   rect = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
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
    howManyTimes <- stats::aggregate(read_id ~ seqnames,
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
  uniqueReadsByID <- stats::aggregate(
    stats::formula(paste0("read_share_percentage_normal ~ ", by, " + bins")),
    data = uniqueReadsByID, sum)
  colorPalette <- grDevices::colorRampPalette(colors)(
    length(levels(uniqueReadsByID$bins)))
  names(colorPalette) <- levels(uniqueReadsByID$bins)
  ggplot2::ggplot(data = uniqueReadsByID,
                  ggplot2::aes_string(x = paste0("as.factor(", by, ")"),
                                      y = "read_share_percentage_normal",
                                      fill = "bins",
                                      order =  paste0("as.factor(", by, ")"))) +
    ggplot2::geom_bar(position='stack', stat='identity') +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 14, face = 'bold'),
                   legend.position = 'top',
                   legend.direction = 'horizontal',
                   legend.title = ggplot2::element_blank()) +
    ggplot2::ylab('Unique reads percentage of shares') +
    ggplot2::xlab(level) +
    ggplot2::scale_fill_manual(values = colorPalette) +
    ggplot2::coord_flip()
}


aa_frame <- function(amplicon, sense = TRUE, frame = 1, ymin = 0, ymax = 1) {
  amplicon_pos <- seq_len(length(amplicon))

  if (sense) {
    xmax <- seq(-1 + frame, length(amplicon), 3)[-1]
    xmin <- xmax - 3
  } else {
    xmin <- rev(seq(length(amplicon) - frame + 1, 0, -3)[-1])
    xmax <- xmin + 3
  }
  f_t <- min(xmin + 1):max(xmax)
  f_e <- amplicon_pos[!amplicon_pos %in% f_t]
  s3_f <- decode(if (sense) amplicon[f_t] else revComp(amplicon[f_t]))

  data.frame(xmin = c(xmin, f_e - 1),
                      xmax = c(xmax, f_e),
                      ymin = ymin,
                      ymax = ymax,
                      codon = c(s3_f, rep("", length(f_e))))
}
range01 <- function(x) (x - min(x))/diff(range(x))
cRamp <- function(x){
  cols <- grDevices::colorRamp(c("#FFFFFF", "#98DDDE"))(range01(x))
  apply(cols, 1, function (xt) grDevices::rgb(xt[1], xt[2], xt[3],
                                              maxColorValue = 255))
}
cRampF <- function(x) {
  greenF <- rep("#79C753", length(x))
  greenF[x %% 3 == 0] <- "#FFFFFF"
  greenF
}


#' Plots most frequent variants using ggplot2 and ggbio.
#'
#' This function plots variants in relation to the amplicon. Shows sequences of
#' top mutants without agregating on deletions, insertions and mismatches.
#'
#' Top plot shows all six possible frames for given amplicon. Amino acids are
#' colored as follows:\cr
#' \tabular{rrrrr}{
#' Small nonpolar \tab G, A, S, T \tab Orange\cr
#' Hydrophobic \tab C, V, I, L, P, F, Y, M, W	\tab Green\cr
#' Polar \tab N, Q, H \tab Magenta\cr
#' Negatively charged \tab D, E \tab Red\cr
#' Positively charged \tab K, R \tab Blue\cr
#' Other \tab eg. *, U, + \tab Grey
#' }
#' Variant plot shows amplicon reference, UPPER letters which were the basis for
#' window selection are highlighted with dashed white box (guideRNA). Black
#' triangles are reflecting insertion points. Dashed letters indicate deletions.
#' Table associated with variant plot represents:
#' \itemize{
#' \item{Freq - }{Frequency of given read in experiment. Variants are ordered by
#' frequency value.}
#' \item{Count - }{Represents raw count of this variant reads in experiment.}
#' \item{Score - }{Alignment score from alignment to the amplicon sequence.}
#' \item{F - }{Sum of deletion and insertion widths over presented window. Green
#' background indicates that frameshift is aparent in presented window. Be weary
#' that this frameshift does not takes into account possible insetions and
#' deletions that happen outside of presented window.}}
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config
#' file or name of multiple IDs if it is possible to group them. First amplicon
#' will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as
#' used in the analysis.
#' @param top (numeric) Specify number of most frequent reads to plot. By
#' default it is 10. Check \code{\link{plot_heterogeneity}} to see how many
#' reads will be enough to give good overview of your variants.
#' @return (variant plot) ggplot2 object of variants plot
#' @export
#' @family specialized plots
#' @note
#' This function is inspired by \code{\link[CrispRVariants]{plotAlignments}}.
#' @examples
#' #example config
#' config <- read.csv(system.file("extdata", "results", "config_summary.csv",
#'                                package = "amplican"))
#' #example alignments results
#' alignments_file <- system.file("extdata", "results", "alignments",
#'                                "events_filtered_shifted_normalized.csv",
#'                                package = "amplican")
#' alignments <- read.csv(alignments_file)
#' plot_variants(alignments, config, c('ID_1','ID_3'))
#'
plot_variants <- function(alignments, config, id,
                          cut_buffer = 5, top = 10) {

  archRanges <- alignments[alignments$seqnames %in% id, ]
  archRanges <- archRanges[archRanges$overlaps, ]

  if (length(unique(archRanges$strand)) == 2) {
    archRanges <- archRanges[archRanges$strand == "+", ]
  }

  if (dim(archRanges)[1] == 0) {
    return("No variants to plot.")
  }

  if (!any(colnames(archRanges) == "frequency")) {
    archRanges$frequency <- archRanges$counts /
      config$Reads_noPD[match(archRanges$seqnames, config$ID)]
  }

  archRanges$read_names <- paste0(archRanges$seqnames, ":", archRanges$read_id)
  archRanges <- archRanges[order(-archRanges$frequency, archRanges$read_names),]

  if (length(unique(archRanges$read_names)) < top) {
    top <- length(unique(archRanges$read_names))
  }

  amplicon <- get_amplicon(config, id)
  box <- upperGroups(amplicon)[1]
  if (length(box) == 1) {
    box_shift <- IRanges::start(box)[1]
    upperBox <- IRanges::start(box):IRanges::end(box) - box_shift
    box <- box + cut_buffer
  } else {
    box_shift <- 0
    upperBox <- 1:nchar(amplicon)
    box <- IRanges::IRanges(upperBox)
  }
  amplicon <- strsplit(amplicon, "")[[1]][IRanges::start(box):IRanges::end(box)]
  box <- IRanges::shift(box, -1 * box_shift)

  xaxis <- IRanges::start(box[1]):IRanges::end(box[1])
  yaxis <- seq_len(top + 1) # + amplicon reference
  yaxis_names <- c("amplicon", unique(archRanges$read_names)[seq_len(top)])
  archRanges <- archRanges[archRanges$read_names %in% yaxis_names, ]
  variants <- matrix(toupper(amplicon),
                     nrow = length(yaxis),
                     ncol = length(amplicon), byrow = TRUE)

  archRanges <- GenomicRanges::restrict(GenomicRanges::GRanges(archRanges),
                                        start = xaxis[1],
                                        end = xaxis[length(xaxis)])
  archRanges <- data.frame(archRanges)
  # deletions and mismatches
  for (i in seq_len(dim(archRanges)[1])) {
    if (archRanges[i, "type"] == "insertion") next()
    if (archRanges[i, "type"] == "mismatch") {
      variants[which(archRanges[i, "read_names"] == yaxis_names),
               which(archRanges[i, "start"] == xaxis)] <-
        as.character(archRanges[i, "replacement"])
    }
    if (archRanges[i, "type"] == "deletion") {
      variants[which(archRanges[i, "read_names"] == yaxis_names),
               which(archRanges[i, "start"] == xaxis):
                 which(archRanges[i, "end"] == xaxis)] <- "-"
    }
  }
  colnames(variants) <- xaxis
  rownames(variants) <- length(yaxis):1
  variants_melt <- reshape2::melt(variants)
  variants_melt$ymin <- variants_melt$Var1 - 1
  variants_melt$ymax <- variants_melt$Var1
  variants_melt$xmin <- variants_melt$Var2 - 0.5
  variants_melt$xmax <- variants_melt$Var2 + 0.5

  insertion_melt <- archRanges[archRanges$type == "insertion", ]
  insertion_melt$y <- match(insertion_melt$read_names, rev(yaxis_names))
  insertion_melt$x <- insertion_melt$start + 0.5

  x<-xlab<-xmax<-xmin<-y<-ylab<-ymax<-ymin<-value<-codon<-NULL
  vplot <- ggplot2::ggplot(variants_melt,
                           ggplot2::aes((xmin + xmax) / 2,
                                        (ymin + ymax) / 2)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax,
                                    fill = value)) +
    ggplot2::geom_hline(yintercept  = 1:(length(yaxis) - 1),
                        colour = "white", linetype = "dashed") +
    ggplot2::geom_rect(data = data.frame(xmax = max(upperBox) + 0.5,
                                         xmin = min(upperBox) - 0.5,
                                         ymax = length(yaxis),
                                         ymin = length(yaxis) - 1),
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax),
                       colour = "black", alpha = 0) +
    ggplot2::geom_point(data = insertion_melt,
                        ggplot2::aes(x = x, y = y), shape = 25,
                        size = 5, fill = "black") +
    ggplot2::geom_text(ggplot2::aes(label = value)) +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none",
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(breaks = (length(yaxis)-1):0 + 0.5,
                                labels = yaxis_names,
                                expand = c(0,0)) +
    ggplot2::labs(y = "Variant",
                  x = "Nucleotide Position Relative to PAM")

  codon_melt <- rbind(aa_frame(amplicon, TRUE, 1, 0, 1),
                      aa_frame(amplicon, TRUE, 2, 1, 2),
                      aa_frame(amplicon, TRUE, 3, 2, 3),
                      aa_frame(amplicon, FALSE, 1, 3, 4),
                      aa_frame(amplicon, FALSE, 2, 4, 5),
                      aa_frame(amplicon, FALSE, 3, 5, 6))
  fnames <- c("1st, 5' <- 3'", "2nd, 5' <- 3'", "3rd, 5' <- 3'",
              "1st, 3' -> 5'", "2nd, 3' -> 5'", "3rd, 3' -> 5'")
  cplot <- ggplot2::ggplot(codon_melt,
                           ggplot2::aes((xmin + xmax) / 2,
                                        (ymin + ymax) / 2)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax,
                                    fill = codon), colour = "#FFFFFF") +
    ggplot2::geom_text(ggplot2::aes(label = codon)) +
    ggplot2::scale_y_continuous(breaks = 0:5 + 0.5,
                                labels = fnames) +
    ggplot2::scale_fill_manual(values = codon_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(y = "Frame")

  vtable <- archRanges[, c("read_names", "frequency", "counts", "score")]
  widthT <- archRanges[archRanges$type != "mismatch", ]
  widthT <- stats::aggregate(width ~ read_names, widthT, sum)
  vtable$Frameshift <- widthT$width[match(vtable$read_names, widthT$read_names)]
  vtable <- vtable[!duplicated(vtable), ]
  vtable <- vtable[, -1]
  vtable$frequency <- round(vtable$frequency, 2)
  colnames(vtable) <- c("Freq", "Count", "Score", "F")

  tgb <- gridExtra::tableGrob(
    vtable, theme = gridExtra::ttheme_minimal(core = list(
      bg_params = list(fill = c(cRamp(vtable$Freq), cRamp(vtable$Count),
                                cRamp(vtable$Score), cRampF(vtable$F)),
                       col = NA))), rows = NULL)
  # tgb <- gtable::gtable_add_grob(tgb,
  #                      grobs = grid::rectGrob(gp = grid::gpar(fill = NA)),
  #                      t = 2, b = nrow(tgb), l = 1, r = ncol(tgb))
  # tgb <- gtable::gtable_add_grob(tgb,
  #                      grobs = grid::rectGrob(gp = grid::gpar(fill = NA)),
  #                      t = 1, l = 1, r = ncol(tgb))
  separators <- replicate(ncol(tgb) - 1,
                          grid::segmentsGrob(x1 = grid::unit(0, "npc"),
                                             gp=grid::gpar(lty = 2)),
                          simplify=FALSE)
  tgb <- gtable::gtable_add_grob(tgb, grobs = separators,
                               t = 1, b = nrow(tgb), l = 2:4)
  separators <- replicate(nrow(tgb) - 1,
                          grid::segmentsGrob(y1 = grid::unit(0, "npc"),
                                             gp=grid::gpar(lty = 2)),
                          simplify=FALSE)
  tgb <- gtable::gtable_add_grob(tgb, grobs = separators,
                                 t = 1:(nrow(tgb) - 1), l = 1, r = 4)
  cgb <- ggplot2::ggplotGrob(cplot)
  vgb <- ggplot2::ggplotGrob(vplot)
  egb <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(1,1), colour="white") +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
  egb <- ggplot2::ggplotGrob(egb)

  maxWidth <- grid::unit.pmax(vgb$widths, cgb$widths)
  vgb$widths <- as.list(maxWidth)
  cgb$widths <- as.list(maxWidth)
  tgb$heights <- grid::unit(rep(1/(nrow(tgb)), nrow(tgb)), "npc")

  bot <- gtable::gtable_add_cols(vgb, sum(tgb$widths))
  bot <- gtable::gtable_add_grob(bot, grobs = tgb,
                                 t = 6, l = ncol(bot), b = 6, r = ncol(bot))
  top <- gtable::gtable_add_cols(cgb, sum(tgb$widths))
  top <- gtable::gtable_add_grob(top, grobs = egb,
                                 t = 6, l = ncol(top), b = 6, r = ncol(top))
  fin <- gridExtra::grid.arrange(
    top, bot,
    heights = grid::unit.c(grid::unit(8, "char"),
                           grid::unit(1, "npc") - grid::unit(8, "char")),
    nrow = 2)
  grid::grid.newpage()
  grid::grid.draw(fin)
}
