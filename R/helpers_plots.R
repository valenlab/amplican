#' @include helpers_general.R
NULL

# scaling for the bezier archs, euation to calculate height relative to y
# simplifies to y = 4/3 * h
sh <- 1.33

# orange #E15D44
# red #BC243C
# blue #98B4D4
# magenta #D65076
# green #009B77
# grey #DFCFBE

# A "#EA4335" red
# C "#4285F4" blue
# G "#FBBC05"  yellow
# T "#34A853" green
amplicon_colors <- c("#EA4335", "#4285F4", "#FBBC05", "#34A853",
                     "#EA4335", "#4285F4", "#FBBC05", "#34A853", "#FFFFFF")
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


amplican_xlim <- function(p, xlabels, box, primers, limits) {
  p +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xlabels,
                                limits = limits) +
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
  frequency <- overlaps <- x <- y <- group <- NULL
  xr <- nrow(xData)
  xData <- data.frame(x = c(xData$start, xData$start, xData$end, xData$end),
                      y = c(rep(0, xr), xData$frequency*sh,
                            xData$frequency*sh, rep(0, xr)),
                      group = rep(seq_len(xr), 4),
                      overlaps = rep(xData$overlaps, 4),
                      frequency = rep(xData$frequency, 4))
  ggplot2::ggplot() +
    geom_bezier(ggplot2::aes(x = x, y = y, group = group,
                                      alpha = frequency,
                                      colour = overlaps,
                                      size = frequency),
                         data = xData) +
    ggplot2::scale_colour_manual(values = is_cut_colors) +
    ggplot2::xlab("Relative Nucleotide Position")
}


ggplot_insertions <- function(xData) {
  position <- frequency <- group <- frequencyReal <- NULL
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = xData,
                          ggplot2::aes(x = position,
                                       y = frequency,
                                       group = group,
                                       alpha = frequencyReal,
                                       size = frequency),
                          fill = "#FF0000")
}


triangulate_ranges <- function(xRanges) {
  if (dim(xRanges)[1] != 0) {
    ifR <- ifr <- rep(xRanges$frequency, each = 3)
    ifr[c(TRUE, FALSE, FALSE)] <- 0
    data.frame(frequency = ifr,
               frequencyReal = ifR,
               position = as.vector(rbind(xRanges$start - 0.5,
                                          xRanges$start - 0.5,
                                          xRanges$end + 0.5)), # ins start 108
               group = rep(1:dim(xRanges)[1], each = 3))       # means it's ins
  } else {                                                     # between 107/108
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
             counts = rep(0, how_many),
             frequency = rep(0, how_many))
}


scale_freq <- function(p, freq) {
  p +
    ggplot2::scale_alpha(limits = c(
      0, max(freq, na.rm = TRUE))) +
    ggplot2::scale_size(limits = c(
      0, max(freq, na.rm = TRUE))) +
    ggplot2::coord_cartesian(
      ylim = c(0, max(freq, na.rm = TRUE)))
}


return_metaplot <- function(freqAgr, plot_fr, plot_re) {
  if (any(freqAgr$strand == "+") & any(freqAgr$strand == "-")) {
    return(gridExtra::grid.arrange(
      scale_freq(plot_fr, freqAgr$frequency) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank()),
      scale_freq(plot_re, freqAgr$frequency) +
        ggplot2::scale_y_reverse() +
        ggplot2::xlab("Relative Nucleotide Position"),
      ncol = 1,
      heights = c(0.5, 0.5),
      padding = -1))
  } else if (all(freqAgr$strand == "+")) {
    return(plot_fr +
             ggplot2::xlab("Relative Nucleotide Position"))
  } else {
    return(plot_re + + ggplot2::scale_y_reverse() +
             ggplot2::xlab("Relative Nucleotide Position"))
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
      colour = amplicon_colors[
        match(toupper(amplicon), names(amplicon_colors))])
}

return_plot <- function(freqAgr, amplicon, from, to, plot_fr, plot_re) {
  if (any(freqAgr$strand == "+") & any(freqAgr$strand == "-")) {
    plot_fr <- scale_freq(plot_fr, freqAgr$frequency) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    plot_re <- scale_freq(plot_re, freqAgr$frequency) +
      ggplot2::scale_y_reverse() +
      ggplot2::xlab("Relative Nucleotide Position")
    amplicon <- plot_amplicon(amplicon, from, to)
    plot_fr <- ggplot2::ggplotGrob(plot_fr)
    amplicon <- ggplot2::ggplotGrob(amplicon)
    plot_re <- ggplot2::ggplotGrob(plot_re)
    g <- gridExtra::gtable_rbind(plot_fr, amplicon, plot_re)
    len_t <- length(plot_fr$heights)
    len_a <- length(amplicon$heights)
    g$heights[c(len_t, len_t + len_a, len_t + len_a + 1)] <-
      grid::unit(0, "cm")
    g$heights[len_t + which(as.character(amplicon$heights) == "1null")] <-
      grid::unit(1.5, "char")
  } else if (all(freqAgr$strand == "+")) {
    g <- annotate_with_amplicon(plot_fr, amplicon, from, to)
    g <- ggplot2::ggplotGrob(g)
  } else {
    g <- annotate_with_amplicon(
      plot_re + ggplot2::scale_y_reverse(), amplicon, from, to)
    g <- ggplot2::ggplotGrob(g)
  }
  grid::grid.newpage()
  grid::grid.draw(g)
  return(g)
}


#' MetaPlots mismatches using ggplot2.
#'
#' Plots mismatches in relation to the amplicons for given
#' \code{selection} vector that groups values by given config \code{group}.
#' All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{amplicanMap}}.
#' Zero position on new coordinates is the most left UPPER case letter of
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
#' metaplot_mismatches(alignments,
#'                     config, "Group", "Betty")
#'
metaplot_mismatches <- function(alnmt, config, group, selection) {
  frequency <- NULL
  alnmt <- alnmt[alnmt$type == "mismatch", ]
  if (length(alnmt) == 0) return("No mismatches to plot.")
  alnmt[[group]] <- config[[group]][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)
  if (dim(alnmt)[1] == 0) return("No mismatches to plot.")

  data.table::setDT(alnmt)
  freqAgr <- alnmt[, list(counts = sum(counts)),
                   by = c("replacement", "start", "strand")]
  freqAgr$frequency <-
    freqAgr$counts/sum(config$Reads_Filtered[config[[group]] %in% selection])

  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  # sometimes barplot gets confused so we add mock data
  mock <- mock_mm_df(max(freqAgr$start, na.rm = TRUE),
                     min(freqAgr$start, na.rm = TRUE))
  freqAgrPlus <- rbind(freqAgrPlus, mock)
  freqAgrMinus <- rbind(freqAgrMinus, mock)

  mut_fr <- ggplot_mismatches(freqAgrPlus)
  mut_fr <- amplican_style(mut_fr)

  mut_re <- ggplot_mismatches(freqAgrMinus)
  mut_re <- amplican_style(mut_re)

  return_metaplot(freqAgr, mut_fr, mut_re)
}


#' MetaPlots deletions using ggplot2.
#'
#' This function plots deletions in relation to the amplicons for given
#' \code{selection} vector that groups values by given config \code{group}.
#' All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{amplicanMap}}.
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
#' @param over (string) Specify which column contains overlaps with
#' expected cut sites generated by \code{\link{amplicanOverlap}}
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
#' metaplot_deletions(alignments[alignments$consensus, ],
#'                    config, "Group", "Betty")
#'
metaplot_deletions <- function(alnmt, config, group,
                              selection, over = "overlaps") {
  frequency <- NULL
  alnmt <- alnmt[alnmt$type == "deletion",]
  alnmt[[group]] <- config[[group]][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)

  if (dim(alnmt)[1] == 0) return("No deletions to plot.")
  data.table::setDT(alnmt)
  archRanges <- alnmt[, list(counts = sum(counts),
                             overlaps = sum(get(over)) > 0),
                      by = c("strand", "start", "end")]
  archRanges$frequency <-
    archRanges$counts/sum(config$Reads_Filtered[config[[group]] %in% selection])

  arch_plot_fr <- ggplot_deletions(archRanges[archRanges$strand == "+", ])
  arch_plot_fr <- amplican_style(arch_plot_fr)

  arch_plot_re <- ggplot_deletions(archRanges[archRanges$strand == "-", ])
  arch_plot_re <- amplican_style(arch_plot_re)

  return_metaplot(archRanges, arch_plot_fr, arch_plot_re)
}


#' MetaPlots insertions using ggplot2.
#'
#' This function plots insertions in relation to the amplicons for given
#' \code{selection} vector that groups values by given config \code{group}.
#' All reads should
#' already be converted to their relative position to their respective amplicon
#' using \code{\link{amplicanMap}}.
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
#' metaplot_insertions(alignments[alignments$consensus, ], config,
#'                     "Group", "Betty")
#'
metaplot_insertions <- function(alnmt, config, group, selection) {
  frequency <- NULL
  alnmt <- alnmt[alnmt$type == "insertion",]
  alnmt[,group] <- config[[group]][match(alnmt$seqnames, config$ID)]
  alnmt <- group_to_selection(alnmt, config, group, selection)

  if (dim(alnmt)[1] == 0) return("No insertions to plot.")
  data.table::setDT(alnmt)
  idRangesReduced <- alnmt[, list(counts = sum(counts)),
                           by = c("strand", "start", "end")]
  idRangesReduced$frequency <-
    idRangesReduced$counts/sum(
      config$Reads_Filtered[config[[group]] %in% selection])

  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+", ]
  triangleFr <- triangulate_ranges(idRangesFr)
  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-", ]
  triangleRe <- triangulate_ranges(idRangesRe)

  ins_fr <- ggplot_insertions(triangleFr)
  ins_fr <- amplican_style(ins_fr)

  ins_re <- ggplot_insertions(triangleRe)
  ins_re <- amplican_style(ins_re)

  return_metaplot(idRangesReduced, ins_fr, ins_re)
}


#' Plots amplicon sequence using ggplot2.
#'
#' @keywords internal
#' @param amplicon (character) Sequence of the amplicon to plot.
#' @param from (number) Minimum on x axis - start of the amplicon
#' @param to (number) Maximum on x axis - not necessarily end of the amplicon
#' @return (amplicon plot) ggplot2 object of amplicon plot
#'
plot_amplicon <- function(amplicon, from, to) {

  ampl_df <- data.frame(position = seq(from, by = 1,
                                       length.out = nchar(amplicon)),
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
    ggplot2::xlim(from, to) + # check me
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


#' Plots mismatches using ggplot2.
#'
#' Plots mismatches in relation to the amplicon, assumes
#' your reads are relative to the respective amplicon sequences predicted cut
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
#' @return (mismatches plot) gtable object of mismatches plot
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
#' p <- plot_mismatches(alignments, config, c('ID_1', 'ID_3'))
#'
plot_mismatches <- function(alignments,
                            config,
                            id,
                            cut_buffer = 5,
                            xlab_spacing = 4) {

  idRanges <- alignments[alignments$seqnames %in% id, ]
  idRanges <- idRanges[idRanges$type == "mismatch", ]
  if (dim(idRanges)[1] == 0) return("No mismatches to plot.")

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

  data.table::setDT(idRanges)
  frequency <- NULL
  freqAgr <- idRanges[, list(counts = sum(counts)),
                      by = c("replacement", "start", "strand")]
  freqAgr$frequency <- freqAgr$counts/sum(config$Reads_Filtered[config$ID %in% id])
  freqAgrPlus <- freqAgr[freqAgr$strand == "+", ]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-", ]

  # sometimes barplot gets confused so we add mock data
  mock <- mock_mm_df(to)
  freqAgrPlus <- rbind(freqAgrPlus, mock)
  freqAgrMinus <- rbind(freqAgrMinus, mock)

  mut_fr <- ggplot_mismatches(freqAgrPlus)
  mut_fr <- amplican_xlim(mut_fr, xlabels, box, pr$primers, c(from, to))
  mut_fr <- amplican_style(mut_fr)

  mut_re <- ggplot_mismatches(freqAgrMinus)
  mut_re <- amplican_xlim(mut_re, xlabels, box, pr$primers, c(from, to))
  mut_re <- amplican_style(mut_re)

  return_plot(freqAgr, amplicon, from, to, mut_fr, mut_re)
}


#' Plots deletions using ggplot2.
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
#' @param over (string) Specify which columns contains overlaps with
#' expected cut sites generated by \code{\link{amplicanOverlap}}
#' @return (deletions plot) gtable object of deletions plot
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
#' p <- plot_deletions(alignments[alignments$consensus, ],
#'                     config, c('ID_1','ID_3'))
#'
plot_deletions <- function(alignments,
                           config,
                           id,
                           cut_buffer = 5,
                           xlab_spacing = 4,
                           over = "overlaps") {
  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]
  if (dim(archRanges)[1] == 0) return("No deletions to plot.")

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
  data.table::setDT(archRanges)
  frequency <- NULL
  archRanges <- archRanges[, list(counts = sum(counts),
                                  overlaps = sum(get(over)) > 0),
                           by = c("strand", "start", "end")]
  archRanges$frequency <-
    archRanges$counts/sum(config$Reads_Filtered[config$ID %in% id])

  arch_plot_fr <- ggplot_deletions(archRanges[archRanges$strand == "+", ])
  arch_plot_fr <- amplican_xlim(arch_plot_fr, xlabels, box,
                                pr$primers, c(from, to))
  arch_plot_fr <- amplican_style(arch_plot_fr)

  arch_plot_re <- ggplot_deletions(archRanges[archRanges$strand == "-", ])
  arch_plot_re <- amplican_xlim(arch_plot_re, xlabels, box,
                                pr$primers, c(from, to))
  arch_plot_re <- amplican_style(arch_plot_re)

  return_plot(archRanges, amplicon, from, to, arch_plot_fr, arch_plot_re)
}


#' Plots insertions using ggplot2.
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
#' @return (insertions plot) gtable object of insertions plot
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
#' p <- plot_insertions(alignments, config, c('ID_1','ID_3'))
#'
plot_insertions <- function(alignments,
                            config,
                            id,
                            cut_buffer = 5,
                            xlab_spacing = 4) {

  idRanges <- alignments[alignments$seqnames %in% id &
                           alignments$type == "insertion", ]
  if (dim(idRanges)[1] == 0) return("No insertions to plot.")

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)
  data.table::setDT(idRanges)
  idRangesReduced <- idRanges[, list(counts = sum(counts)),
                              by = c("strand", "start", "end")]
  idRangesReduced$frequency <-
    idRangesReduced$counts/sum(config$Reads_Filtered[config$ID %in% id])

  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+", ]
  triangleFr <- triangulate_ranges(idRangesFr)
  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-", ]
  triangleRe <- triangulate_ranges(idRangesRe)

  if (dim(idRangesRe)[1] != 0 | dim(idRangesFr)[1] != 0) {
    ampl_len <- max(c(max(c(0, triangleFr$position), na.rm = TRUE),
                      max(c(0, triangleRe$position), na.rm = TRUE),
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

  ins_fr <- ggplot_insertions(triangleFr)
  ins_fr <- amplican_xlim(ins_fr, xlabels, box, pr$primers, c(from, to))
  ins_fr <- amplican_style(ins_fr)

  ins_re <- ggplot_insertions(triangleRe)
  ins_re <- amplican_xlim(ins_re, xlabels, box, pr$primers, c(from, to))
  ins_re <- amplican_style(ins_re)

  return_plot(idRangesReduced, amplicon, from, to, ins_fr, ins_re)
}


#' Plots cuts using ggplot2.
#'
#' This function plots cuts in relation to the amplicon with distinction for
#' each ID.
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
#' @return (cuts plot) gtable object of cuts plot
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
#' plot_cuts(alignments[alignments$consensus & alignments$overlaps, ],
#'           config, c('ID_1','ID_3'))
#'
plot_cuts <- function(alignments,
                      config,
                      id,
                      cut_buffer = 5,
                      xlab_spacing = 4) {

  archRanges <- alignments[alignments$seqnames %in% id &
                             alignments$type == "deletion", ]
  archRanges <- archRanges[archRanges$strand == "+", ]

  amplicon <- get_amplicon(config, id)
  ampl_len <- nchar(amplicon)

  if (dim(archRanges)[1] == 0) return("No cuts to plot.")
  data.table::setDT(archRanges)
  frequency <- NULL
  archRanges <- archRanges[, list(counts = sum(counts)),
                           by = c("strand", "start", "end", "seqnames")]
  archRanges$frequency <- archRanges$counts /
    config$Reads_Filtered[match(archRanges$seqnames, config$ID)]
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
  if (dim(archRanges)[1] == 0) return("No cuts to plot.")

  amplicon <- strsplit(amplicon, "")[[1]]
  frequency <- seqnames <- x <- y <- group <- NULL

  xr <- nrow(archRanges)
  archRanges <- data.frame(x = c(archRanges$start, archRanges$start,
                                 archRanges$end, archRanges$end),
                           y = c(rep(0, xr), archRanges$frequency*sh,
                                 archRanges$frequency*sh, rep(0, xr)),
                           group = rep(seq_len(xr), 4),
                           seqnames = rep(archRanges$seqnames, 4),
                           frequency = rep(archRanges$frequency, 4))
  p <- ggplot2::ggplot() +
    geom_bezier(ggplot2::aes(x= x, y = y, group = group,
                                      alpha = frequency,
                                      colour = seqnames,
                                      size = frequency),
                         data = archRanges) +
    ggplot2::theme_bw() +
    ggplot2::guides(size = FALSE, alpha = FALSE) +
    ggplot2::theme(legend.position = c(1, 1),
                   legend.justification = c(1.01, 1.01)) +
    ggplot2::labs(y = "Frequency [%]",
                  colour = "Experiments",
                  x = "Relative Nucleotide Position") +
    ggplot2::scale_x_continuous(labels = xlabels,
                                breaks = xlabels,
                                limits = c(from, to)) +
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


#' Plots heterogeneity of the reads using ggplot2.
#'
#' This function creates stacked barplot explaining
#' reads heterogeneity. It groups reads
#' by user defined levels and measures how unique are reads in
#' this level. Uniqueness of reads is simplified to the bins and
#' colored according to the color gradient. Default color black
#' indicates very high heterogeneity of the reads. The more yellow (default)
#' the more similar are reads and less heterogeneous.
#'
#' @param alignments (data.frame) Loaded alignment information from
#' alignments_events.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param level (string) Name of the column from config
#' file specifying levels to group by.
#' @param colors (html colors vector) Two colours for gradient,
#' eg. c('#000000', '#F0E442').
#' @param bins (numeric vector) Numeric vector from 0 to 100 specifying bins eg.
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
#' plot_heterogeneity(alignments[alignments$consensus, ], config)
#'
plot_heterogeneity <- function(alignments,
                               config,
                               level = "ID",
                               colors = c('#000000', '#F0E442'),
                               bins = c(0, 5, seq(10, 100, 10))) {
  seqnames <- read_id <- replacement <- read_shares <- NULL
  #long into wide
  data.table::setDT(alignments)
  alignments <- alignments[
    order(seqnames, read_id, type, start, end, replacement), ]
  alignments <- alignments[, `:=` (event_id = as.numeric(seq_len(.N))),
                           by = c("seqnames", "read_id")]
  alignments <- data.table::dcast(alignments,
                                  seqnames + read_id + counts ~ event_id,
                                  value.var = c(
                                    "start", "end", "type", "replacement"))
  alignments[[level]] <- config[[level]][match(alignments$seqnames, config$ID)]
  # collapse reads
  cols <- colnames(alignments)[!colnames(alignments) %in%
                                 c("seqnames", "read_id", "counts")]
  alignments <- alignments[, list(counts = sum(counts)), by = cols]
  alignments <- alignments[, `:=` (read_id = as.numeric(seq_len(.N))),
                           by = level]
  alignments <- alignments[, c(level, "counts", "read_id"), with = FALSE]

  c_ord <- order(alignments[[level]], alignments$counts, decreasing = TRUE)
  alignments <- alignments[c_ord, ]
  alignments[, read_shares:= counts*100/sum(counts),
             by = level]

  # divide into bins for colour
  alignments$bins <- cut(alignments$read_shares, bins)
  # reduce number of reads in 0-5 group - faster plots without artifacts
  alignments <- alignments[, list(read_shares = sum(read_shares)),
                                  by = c(level, "bins")]
  colorPalette <- grDevices::colorRampPalette(colors)(
    length(levels(alignments$bins)))
  names(colorPalette) <- levels(alignments$bins)
  ggplot2::ggplot(data = alignments,
                  ggplot2::aes_string(
                    x = paste0("as.factor(", level, ")"),
                    y = "read_shares",
                    fill = "bins",
                    order =  paste0("as.factor(", level, ")"))) +
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
range01 <- function(x) {
  nx <- (x - min(x))/diff(range(x))
  nx[!is.finite(nx)] <- 0
  nx
}
cRamp <- function(x){
  cols <- grDevices::colorRamp(c("#FFFFFF", "#98DDDE"))(range01(x))
  apply(cols, 1, function (xt) {
    grDevices::rgb(xt[1], xt[2], xt[3], maxColorValue = 255)
  })
}
cRampF <- function(x) {
  greenF <- rep("#79C753", length(x))
  greenF[x %% 3 == 0] <- "#FFFFFF"
  greenF
}
merge_4_grobs <- function(cgb, vgb, tgb, sgb) {
  maxWidth <- grid::unit.pmax(vgb$widths, cgb$widths)
  vgb$widths <- as.list(maxWidth)
  cgb$widths <- as.list(maxWidth)
  tgb$heights <- grid::unit(rep(1/(nrow(tgb)), nrow(tgb)), "npc")

  bot <- gtable::gtable_add_cols(vgb, sum(tgb$widths))
  panel_id <- bot$layout[bot$layout$name == "panel", c("t","b")]
  bot <- gtable::gtable_add_grob(bot, grobs = tgb,
                                 t = panel_id$t, l = ncol(bot),
                                 b = panel_id$b, r = ncol(bot))
  top <- gtable::gtable_add_cols(cgb, sum(tgb$widths))
  top <- gtable::gtable_add_grob(top, grobs = sgb,
                                 t = panel_id$t, l = ncol(top),
                                 b = panel_id$b, r = ncol(top))
  fin <- gridExtra::grid.arrange(
    top, bot,
    heights = grid::unit.c(grid::unit(8, "char"),
                           grid::unit(1, "npc") - grid::unit(8, "char")),
    nrow = 2)
  fin
}
merge_3_grobs <- function(cgb, vgb, tgb) {
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
  panel_id <- bot$layout[bot$layout$name == "panel", c("t","b")]
  bot <- gtable::gtable_add_grob(
    bot, grobs = tgb, t = panel_id$t, l = ncol(bot),
    b = panel_id$b, r = ncol(bot))
  top <- gtable::gtable_add_cols(cgb, sum(tgb$widths))
  top <- gtable::gtable_add_grob(
    top, grobs = egb,
    t = panel_id$t, l = ncol(top),
    b = panel_id$b, r = ncol(top))
  fin <- gridExtra::grid.arrange(
    top, bot,
    heights = grid::unit.c(grid::unit(8, "char"),
                           grid::unit(1, "npc") - grid::unit(8, "char")),
    nrow = 2)
  fin
}
merge_2_grobs <- function(vgb, tgb) {
  tgb$heights <- grid::unit(rep(1/(nrow(tgb)), nrow(tgb)), "npc")
  bot <- gtable::gtable_add_cols(vgb, sum(tgb$widths))
  panel_id <- bot$layout[bot$layout$name == "panel", c("t","b")]
  bot <- gtable::gtable_add_grob(bot, grobs = tgb,
                                 t = panel_id$t, l = ncol(bot), b = panel_id$b,
                                 r = ncol(bot))
  grid::grid.newpage()
  grid::grid.draw(bot)
  bot
}


#' Plots most frequent variants using ggplot2.
#'
#' This function plots variants in relation to the amplicon. Shows sequences of
#' top mutants without aggregating on deletions, insertions and mismatches.
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
#' \item{F - }{Sum of deletion and insertion widths of events overlapping
#' presented window. Green background indicates frameshift. }}
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
#' @param annot ("codon" or NA) What to display for
#' annotation top plot. When NA will not display anything, also not display
#' total summary.
#' @param summary_plot (boolean) Whether small summary plot in the upper right
#' corner should be displayed. Top bar summarizes total reads with
#' frameshift (F), reads with Edits without Frameshift (Edits) and reads
#' without Edits (Match).
#' \preformatted{
#' annot            on    | off
#  summary_plot  on | off | off
#' }
#' @return (variant plot) gtable object of variants plot
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
#' p <- plot_variants(alignments[alignments$consensus & alignments$overlaps, ],
#'                    config, c('ID_1','ID_3'))
#'
plot_variants <- function(alignments, config, id,
                          cut_buffer = 5, top = 10,
                          annot = "codon", summary_plot = TRUE) {
  seqnames <- read_id <- replacement <- NULL

  archRanges <- alignments[alignments$seqnames %in% id, ]
  if (dim(archRanges)[1] == 0) return("No variants to plot.")
  archRanges$strand <- "*"

  amplicon <- get_amplicon(config, id)
  box <- upperGroups(amplicon)[1]
  if (length(box) == 1) {
    box_shift <- IRanges::start(box)[1]
    upperBox <- IRanges::start(box):IRanges::end(box) - box_shift
    box <- box + cut_buffer
    box <- IRanges::restrict(box, start = 1,
                             end = nchar(amplicon))
  } else {
    box_shift <- 0
    upperBox <- 1:nchar(amplicon)
    box <- IRanges::IRanges(upperBox)
  }
  amplicon <- strsplit(amplicon, "")[[1]][IRanges::start(box):IRanges::end(box)]
  box <- IRanges::shift(box, -1 * box_shift)
  xaxis <- IRanges::start(box[1]):IRanges::end(box[1])

  # calculate frameshifts beforehand
  data.table::setDT(archRanges)
  widthT <- archRanges
  widthT[type == "mismatch", width := 0]
  widthT[type == "deletion", width := width * -1L] # by reference
  widthT <- widthT[, list(width = sum(width)), c("seqnames", "read_id")]

  # cast reads from long to wide format
  archRanges <- archRanges[
    order(seqnames, read_id, type, start, end, replacement), ]
  archRanges <- archRanges[, `:=` (event_id = as.numeric(seq_len(.N))),
                           by = c("seqnames", "read_id")]
  archRanges <- data.table::dcast(archRanges,
                                  seqnames + read_id + counts ~ event_id,
                                  value.var = c(
                                    "start", "end", "type", "replacement"))
  # add frameshift info
  data.table::setkeyv(widthT, c("seqnames", "read_id"))
  data.table::setkeyv(archRanges, c("seqnames", "read_id"))
  archRanges <- merge(archRanges, widthT, all=TRUE)
  # collapse reads
  cols <- colnames(archRanges)[!colnames(archRanges) %in%
                                 c("seqnames", "read_id", "counts")]
  archRanges <- archRanges[, list(counts = sum(counts)), by = cols]
  archRanges$frequency <-
    archRanges$counts/sum(config$Reads_Filtered[config$ID %in% id])

  # restrict ranges
  max_event <- max(as.numeric(gsub("start_", "", cols[grepl("start_", cols)])))
  for (i in seq_len(max_event)) {
    na_cases <- is.na(archRanges[[paste0("start_", i)]])
    if (all(na_cases)) next()
    shftGR <- IRanges::restrict(IRanges::IRanges(
      start = archRanges[[paste0("start_", i)]][!na_cases],
      end = archRanges[[paste0("end_", i)]][!na_cases]),
      start = xaxis[1],
      end = xaxis[length(xaxis)], keep.all.ranges = TRUE)
    archRanges[[paste0("start_", i)]][!na_cases] <- IRanges::start(shftGR)
    archRanges[[paste0("end_", i)]][!na_cases] <- IRanges::end(shftGR)
    cols <- paste0(c("start_", "end_", "type_", "replacement_"), i)
    data.table::set(archRanges,
                    i = which(!na_cases)[IRanges::width(shftGR) == 0],
                    j = cols, value = NA)
  }

  if (dim(archRanges)[1] == 0) return("No variants to plot.")
  archRanges <- archRanges[order(-archRanges$frequency), ]
  ampl_freq <- 1 - sum(archRanges$frequency)
  ampl_count <-
    sum(config$Reads_Filtered[config$ID %in% id]) - sum(archRanges$counts)

  if (dim(archRanges)[1] < top) top <- dim(archRanges)[1]
  yaxis <- seq_len(top + 2) # + amplicon reference + empty (column header)
  yaxis_names <- c("", "amplicon", seq_len(top))
  archRanges <- archRanges[seq_len(top), ]
  variants <- matrix(toupper(amplicon),
                     nrow = length(yaxis),
                     ncol = length(amplicon), byrow = TRUE)
  variants[1, ] <- "" # header empty

  insertion_melt <- data.frame()
  # deletions and mismatches
  for (i in seq_len(top)) {
    for (j in seq_len(max_event)) {
      if (is.na(archRanges[[paste0("type_", j)]][i])) next()
      if (archRanges[[paste0("type_", j)]][i] == "insertion") {
        insertion_melt <- rbind(insertion_melt,
                                c(top - i + 1,
                                  archRanges[[paste0("start_", j)]][i] - 0.5))
      }
      if (archRanges[[paste0("type_", j)]][i] == "mismatch") {
        variants[i + 2, which(archRanges[[paste0("start_", j)]][i] == xaxis)] <-
          as.character(archRanges[[paste0("replacement_", j)]][i])
      }
      if (archRanges[[paste0("type_", j)]][i] == "deletion") {
        variants[i + 2, which(archRanges[[paste0("start_", j)]][i] == xaxis):
                   which(archRanges[[paste0("end_", j)]][i] == xaxis)] <- "-"
      }
    }
  }
  colnames(variants) <- xaxis
  rownames(variants) <- length(yaxis):1
  variants_melt <- data.table::melt(variants)
  variants_melt$ymin <- variants_melt$Var1 - 1
  variants_melt$ymax <- variants_melt$Var1
  variants_melt$xmin <- variants_melt$Var2 - 0.5
  variants_melt$xmax <- variants_melt$Var2 + 0.5
  if (dim(insertion_melt)[1] > 0) colnames(insertion_melt) <- c("y", "x")

  x<-xlab<-xmax<-xmin<-y<-ylab<-ymax<-ymin<-value<-codon<-NULL
  vplot <- ggplot2::ggplot(variants_melt,
                           ggplot2::aes((xmin + xmax) / 2,
                                        (ymin + ymax) / 2)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax,
                                    fill = value)) +
    ggplot2::geom_hline(yintercept  = 1:(length(yaxis) - 2),
                        colour = "white", linetype = "dashed") +
    ggplot2::geom_rect(data = data.frame(xmax = max(upperBox) + 0.5,
                                         xmin = min(upperBox) - 0.5,
                                         ymax = length(yaxis) - 1,
                                         ymin = length(yaxis) - 2),
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax),
                       colour = "black", alpha = 0) +
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
                                expand = c(0, 0)) +
    ggplot2::labs(y = trimws(paste0(
      if (length(id) == 1) id else "", collapse = "")),
      x = "Relative Nucleotide Position")

  if (dim(insertion_melt)[1] > 0) {
    vplot <- vplot + ggplot2::geom_point(data = data.frame(insertion_melt),
                                         ggplot2::aes(x = x, y = y), shape = 25,
                                         size = 4,
                                         fill = "black")
  }
  if (!is.na(annot) & annot == "codon") {
    codon_melt <- rbind(aa_frame(amplicon, TRUE, 1, 5, 6),
                        aa_frame(amplicon, TRUE, 2, 4, 5),
                        aa_frame(amplicon, TRUE, 3, 3, 4),
                        aa_frame(amplicon, FALSE, 1, 2, 3),
                        aa_frame(amplicon, FALSE, 2, 1, 2),
                        aa_frame(amplicon, FALSE, 3, 0, 1))
    fnames <- rev(c("1st, 5' -> 3'", "2nd, 5' -> 3'", "3rd, 5' -> 3'",
                    "1st, 3' <- 5'", "2nd, 3' <- 5'", "3rd, 3' <- 5'"))
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
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.margin = grid::unit(
                       c(0, 0, -2, 0), "char")) +
      ggplot2::labs(y = "Frame")
    cgb <- ggplot2::ggplotGrob(cplot)
  }

  vtable <- archRanges[, c("frequency", "counts", "width")]
  vtable <- rbind(
    data.frame(frequency = ampl_freq, counts = ampl_count, width = 0), vtable)
  vtable$frequency <- round(vtable$frequency, 2)
  vtable[is.na(vtable) | !apply(vtable, 2, is.finite)] <- 0
  colnames(vtable) <- c("Freq", "Count", "F")
  tgb <- gridExtra::tableGrob(
    vtable, theme = gridExtra::ttheme_minimal(core = list(
      bg_params = list(
        fill = c(cRamp(vtable$Freq), cRamp(vtable$Count),
                 cRampF(vtable[["F"]])),
        col = NA))), rows = NULL)

  separators <- replicate(ncol(tgb) - 1,
                          grid::segmentsGrob(x1 = grid::unit(0, "npc"),
                                             gp=grid::gpar(lty = 2)),
                          simplify=FALSE)
  tgb <- gtable::gtable_add_grob(tgb, grobs = separators,
                               t = 1, b = nrow(tgb), l = 2:3)
  separators <- replicate(nrow(tgb) - 1,
                          grid::segmentsGrob(y1 = grid::unit(0, "npc"),
                                             gp=grid::gpar(lty = 2)),
                          simplify=FALSE)
  tgb <- gtable::gtable_add_grob(tgb, grobs = separators,
                                 t = 1:(nrow(tgb) - 1), l = 1, r = 3)
  vgb <- ggplot2::ggplotGrob(vplot)

  if (!is.na(annot) & annot == "codon" & summary_plot) {
    bnames <- c("Edited", "Match", "F")
    cfgS <- config[config$ID %in% id,
                   c("Reads_Filtered", "Reads_Edited", "Reads_Frameshifted")]
    cfgS <- colSums(cfgS)
    uT <- c(cfgS[2] - cfgS[3], cfgS[1] - cfgS[2], cfgS[3])/cfgS[1]
    uTo <- c(3, 1, 2)
    cfgS <- data.frame(
      x = round(c(uT * 100)),
      group = factor(bnames,
                     levels = bnames[uTo],
                     ordered = TRUE))
    scols <- c("#7FCDCD", "#898E8C", "#79C753")
    names(scols) <- cfgS$group

    hjust_x <- c(-0.5, -0.5, -0.5)
    hjust_x[uT > 0.8] <- 1.5
    group <- NULL #Check
    splot <- ggplot2::ggplot(cfgS,
                             ggplot2::aes(x = group, y = x, fill = group)) +
      ggplot2::geom_bar(stat='identity') +
      ggplot2::geom_text(ggplot2::aes(label = x), hjust = hjust_x,
                         size = 8*5/15) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(size = 8),
                     axis.title = ggplot2::element_text(size = 8),
                     text = ggplot2::element_text(size = 8),
                     plot.margin = grid::unit(
                       c(1, 1, 2, 1), "char")) +
      ggplot2::scale_y_continuous(position = "right", name = "[ % ]",
                                  limits = c(0, 100)) +
      ggplot2::scale_fill_manual(values = scols) +
      ggplot2::coord_flip()

    sgb <- ggplot2::ggplotGrob(splot)
  }

  if (is.na(annot)) {
    merge_2_grobs(vgb, tgb)
  } else if (summary_plot) {
    merge_4_grobs(cgb, vgb, tgb, sgb)
  } else {
    merge_3_grobs(cgb, vgb, tgb)
  }
}
