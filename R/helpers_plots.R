#' Plots mismatches using ggplot2 and ggbio.
#'
#' This function plots mismatches in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as used in the analysis.
#' @return (mismatches plot) ggplot2 object of mismatches plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual geom_bar
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom ggbio tracks xlim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom stats na.omit aggregate
#' @export
#'
amplican_plot_mismatches <- function(alignments, config, id, cut_buffer = 5) {

  idRanges <- alignments[alignments$seqnames %in% id,]
  idRanges <- idRanges[idRanges$type == "mismatch",]

  if (dim(idRanges)[1] == 0) {return(print("No mismatches to plot."))}

  amplicon <- toString(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)
  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id[1]), "Direction"] != 1) {
    if (length(box) >= 1) {
      seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
    } else {
      seq(1, ampl_len, 4)
    }
  } else {
    if (length(box) >= 1) {
      rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
    } else {
      rev(seq(1, ampl_len, 4))
    }
  }
  xbreaks <- seq(1, ampl_len, 4)
  box <- box + cut_buffer

  frPrimer <- toString(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon), toupper(frPrimer)))
  rwPrimer <- toString(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer)))))))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {frPrimer <- c(1, 1)}
  if (length(rwPrimer) == 0) {rwPrimer <- c(ampl_len, ampl_len)}

  #FILTER EVENTS
  # forward before primer after amplicon length
  idRanges <- idRanges[!(idRanges$strand == "+" & idRanges$start <= frPrimer[1]),]
  idRanges <- idRanges[!(idRanges$strand == "+" & idRanges$end >= ampl_len),]
  # reverse before primer and after amplicon start
  idRanges <- idRanges[!(idRanges$strand == "-" & idRanges$start >= rwPrimer[1]),]
  idRanges <- idRanges[!(idRanges$strand == "-" & idRanges$start == 1),]
  if (dim(idRanges)[1] == 0) {return(print("No mismatches to plot."))}

  amplicon_colors <- c("green3", "red", "gold2", "blue", "green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")

  nucleotide <- position <- count <- upper <- NULL
  amp_plot <- ggplot(ampl_df, aes(x = position, label = nucleotide, colour = upper, y = count)) +
    geom_text(size = I(4)) + ggplot2::ylim(0.7, 1.1) + ggbio::xlim(1, ampl_len) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position="none", legend.margin = unit(0, "cm"), line = element_blank(),
          rect = element_blank(), axis.title.x = element_blank(), axis.text = element_blank(),
          axis.title.y = element_blank()) +
    scale_colour_manual(drop = F, values = amplicon_colors)

  #reduce mismatches
  # Setting the variables to NULL first for retarded CRAN check
  frequency <- mm_replacement <- start <- strand <- NULL
  freqAgr <- stats::aggregate(cbind(count, frequency) ~ mm_replacement + start + strand, idRanges, sum)
  freqAgrPlus <- freqAgr[freqAgr$strand == "+",]
  freqAgrMinus <- freqAgr[freqAgr$strand == "-",]

  if (dim(freqAgrPlus)[1] == 0) {
    freqAgrPlus <- rbind(freqAgrPlus, data.frame(mm_replacement = "G", start = 0, strand = "+", count = 0, frequency = 0))
  }

  if (dim(freqAgrMinus)[1] == 0) {
    freqAgrMinus <- rbind(freqAgrMinus, data.frame(mm_replacement = "G", start = 0, strand = "+", count = 0, frequency = 0))
  }

  mut_fr <- ggplot() +
    geom_bar(data = freqAgrPlus,
             aes(x = as.numeric(start) + 1,
                 y = frequency,
                 width = 1,
                 order = mm_replacement,
                 fill = mm_replacement),
             stat = "identity", position = "identity")  +
    scale_fill_manual(values = amplicon_colors) +
    xlim(1, ampl_len) +
    theme_bw() +
    theme(panel.margin = unit(0, "cm"),
          legend.position="none", legend.margin = unit(0, "cm"),
          panel.grid=element_blank()) +
    ylab("Frequency [%]") +
    scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue")

  mut_re <- ggplot() +
    geom_bar(data = freqAgrMinus,
             aes(x = as.numeric(start) + 1,
                 y = frequency,
                 width = 1,
                 order = mm_replacement,
                 fill = mm_replacement),
             stat = "identity", position = "identity")  +
    scale_fill_manual(values = amplicon_colors) +
    xlim(1, ampl_len) +
    theme_bw() +
    theme(panel.margin = unit(0, "cm"),
          legend.position="none", legend.margin = unit(0, "cm"),
          panel.grid=element_blank()) +
    ylab("Frequency [%]") +
    scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue") +
    scale_y_reverse()

  p <- ggbio::tracks(mut_fr, amp_plot, mut_re,  heights = c(0.5, 0.03, 0.53), padding = -1, xlim = 1:ampl_len,
              xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots deletions using ggplot2 and ggbio.
#'
#' This function plots deletions in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as used in the analysis.
#' @return (deletions plot) ggplot2 object of deletions plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks geom_arch xlim
#' @importFrom stats na.omit
#' @export
#'
amplican_plot_deletions <- function(alignments, config, id, cut_buffer = 5) {

  archRanges <- alignments[alignments$seqnames %in% id & alignments$type == "deletion",]

  amplicon <- toString(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  if (dim(archRanges)[1] == 0) {return(print("No deletions to plot."))}
  archRanges <- stats::aggregate(cbind(count, frequency, cut) ~ strand + start + end, archRanges, sum)

  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id[1]), "Direction"] != 1) {
    if (length(box) >= 1) {
      seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
    } else {
      seq(1, ampl_len, 4)
    }
  } else {
    if (length(box) >= 1) {
      rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
    } else {
      rev(seq(1, ampl_len, 4))
    }
  }
  xbreaks <- seq(1, ampl_len, 4)
  box <- box + cut_buffer

  frPrimer <- toString(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon), toupper(frPrimer)))
  rwPrimer <- toString(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer)))))))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {frPrimer <- c(1, 1)}
  if (length(rwPrimer) == 0) {rwPrimer <- c(ampl_len, ampl_len)}

  #FILTER EVENTS
  # forward before primer after amplicon length
  archRanges <- archRanges[!(archRanges$strand == "+" & archRanges$start <= frPrimer[1]),]
  archRanges <- archRanges[!(archRanges$strand == "+" & archRanges$end >= ampl_len),]
  # reverse before primer and after amplicon start
  archRanges <- archRanges[!(archRanges$strand == "-" & archRanges$start >= rwPrimer[1]),]
  archRanges <- archRanges[!(archRanges$strand == "-" & archRanges$start == 1),]
  if (dim(archRanges)[1] == 0) {return(print("No deletions to plot."))}

  frequency <- cut <- NULL
  arch_plot_fr <- ggplot2::ggplot() + ggbio::geom_arch(data = archRanges[archRanges$strand == "+",],
                                                       ggplot2::aes(height = frequency,
                                                                    alpha = frequency,
                                                                    colour = cut,
                                                                    order = frequency,
                                                                    size = frequency,
                                                                    x = start,
                                                                    xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue")

  arch_plot_re <- ggplot2::ggplot() + ggbio::geom_arch(data = archRanges[archRanges$strand == "-",],
                                                       ggplot2::aes(height = frequency,
                                                                    alpha = frequency,
                                                                    colour = cut,
                                                                    order = frequency,
                                                                    size = frequency,
                                                                    x = start,
                                                                    xend = end)) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue") +
    ggplot2::scale_y_reverse()

  amplicon_colors <- c("green3", "red", "gold2", "blue", "green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")
  position <- nucleotide <- upper <- count <- NULL
  amp_plot <- ggplot2::ggplot(ampl_df, ggplot2::aes(x = position, label = nucleotide, colour = upper, y = count)) +
    ggplot2::geom_text(size = I(4)) + ggplot2::ylim(0.7, 1.1) + ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   legend.position="none", legend.margin = unit(0, "cm"), line = element_blank(),
                   rect = element_blank(), axis.title.x = element_blank(), axis.text = element_blank(),
                   axis.title.y = element_blank()) +
    ggplot2::scale_colour_manual(drop = F, values = amplicon_colors)

  p <- ggbio::tracks(arch_plot_fr, amp_plot, arch_plot_re,
                     heights = c(0.5, 0.03, 0.53), padding = -1, xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots insertions using ggplot2 and ggbio.
#'
#' This function plots insertions in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignments (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @param cut_buffer (numeric) Default is 5, you should specify the same as used in the analysis.
#' @return (insertions plot) ggplot2 object of insertions plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual
#' geom_polygon scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks xlim
#' @importFrom stats na.omit aggregate
#' @export
#'
amplican_plot_insertions <- function(alignments, config, id, cut_buffer = 5) {

  idRanges <- alignments[alignments$seqnames %in% id & alignments$type == "insertion",]

  amplicon <- toString(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  if (dim(idRanges)[1] == 0) { return(print("No insertions to plot.")) }

  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  frPrimer <- toString(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon), toupper(frPrimer)))
  rwPrimer <- toString(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer)))))))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {frPrimer <- c(1, 1)}
  if (length(rwPrimer) == 0) {rwPrimer <- c(ampl_len, ampl_len)}

  idRanges <- idRanges[!(idRanges$strand == "+" & idRanges$start <= frPrimer[1]),]
  idRanges <- idRanges[!(idRanges$strand == "+" & idRanges$end >= ampl_len),]
  # reverse before primer and after amplicon start
  idRanges <- idRanges[!(idRanges$strand == "-" & idRanges$start >= rwPrimer[1]),]
  idRanges <- idRanges[!(idRanges$strand == "-" & idRanges$start == 1),]
  if (dim(idRanges)[1] == 0) { return(print("No insertions to plot.")) }

  #reduce
  idRangesReduced <- stats::aggregate(frequency ~ strand + start + end, idRanges, sum)
  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+",]
  if (dim(idRangesFr)[1] != 0) {
    idRangesFrgroup = rep(1:dim(idRangesFr)[1], each = 3)
    idRangesFrFrequency = rep(idRangesFr$frequency, each = 3)
    idRangesFrFrequency[c(T, F, F)] <- 0
    idRangesFrX <- as.vector(rbind(idRangesFr$start, idRangesFr$start, idRangesFr$end))
    traingleFr <- data.frame(frequency = idRangesFrFrequency, position = idRangesFrX, group = idRangesFrgroup)
  } else {
    idRangesFrX <- c()
    traingleFr <- data.frame(frequency = c(), position = c(), group = c())
  }

  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-",]
  if (dim(idRangesRe)[1] != 0) {
    idRangesRegroup = rep(1:dim(idRangesRe)[1], each = 3)
    idRangesReFrequency = rep(idRangesRe$frequency, each = 3)
    idRangesReFrequency[c(T, F, F)] <- 0
    idRangesReX <- as.vector(rbind(idRangesRe$start, idRangesRe$start, idRangesRe$end))
    traingleRe <- data.frame(frequency = idRangesReFrequency, position = idRangesReX, group = idRangesRegroup)
  } else {
    idRangesReX <- c()
    traingleRe <- data.frame(frequency = c(), position = c(), group = c())
  }

  if (dim(idRangesRe)[1] != 0 | dim(idRangesFr)[1] != 0) {
    ampl_len <- max(c(idRangesFrX, idRangesReX, ampl_len))
  }

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id[1]), "Direction"] != 1) {
    if (length(box) >= 1) {
      seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
    } else {
      seq(1, ampl_len, 4)
    }
  } else {
    if (length(box) >= 1) {
      rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
    } else {
      rev(seq(1, ampl_len, 4))
    }
  }
  xbreaks <- seq(1, ampl_len, 4)
  box <- box + cut_buffer

  amplicon_colors <- c("green3", "red", "gold2", "blue", "green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")
  position <- nucleotide <- upper <- count <- NULL
  amp_plot <- ggplot2::ggplot(ampl_df, aes(x = position, label = nucleotide, colour = upper, y = count)) +
    ggplot2::geom_text(size = I(4)) + ggplot2::ylim(0.7, 1.1) + ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   legend.position="none", legend.margin = unit(0, "cm"), line = element_blank(),
                   rect = element_blank(), axis.title.x = element_blank(), axis.text = element_blank(),
                   axis.title.y = element_blank()) +
    ggplot2::scale_colour_manual(drop = F, values = amplicon_colors)

  frequency <- group <- NULL
  ins_fr <- ggplot2::ggplot() + ggplot2::geom_polygon(data = traingleFr,
                                                      aes(x = position,
                                                          y = frequency,
                                                          group = group,
                                                          alpha = frequency,
                                                          order = frequency,
                                                          size = frequency))  +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    #ylim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue")

  ins_re <- ggplot2::ggplot() + ggplot2::geom_polygon(data = traingleRe,
                                                      ggplot2::aes(x = position,
                                                                   y = frequency,
                                                                   group = group,
                                                                   alpha = frequency,
                                                                   order = frequency,
                                                                   size = frequency))  +
    ggplot2::scale_fill_manual(values = amplicon_colors) +
    ggbio::xlim(1, ampl_len) +
    #ylim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue") +
    ggplot2::scale_y_reverse()

  p <- ggbio::tracks(ins_fr, amp_plot, ins_re,  heights = c(0.5, 0.03, 0.53), padding = -1, xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}


#' Plots cuts using ggplot2 and ggbio.
#'
#' This function plots cuts in relation to the amplicon with distinction for each ID.
#' Top plot is metaplot of all cuts combined with frequee to all reads in the amplicon group.
#' Bottom plot shows cuts with facets for each ID.
#'
#' @param alignments (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string or vector of strings) Name of the ID column from config file or name of multiple IDs if it is
#' possible to group them. First amplicon will be used as the basis for plot.
#' @return (cuts plot) ggplot2 object of cuts plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' scale_color_manual facet_grid element_text
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks geom_arch xlim
#' @importFrom stats na.omit
#' @export
#'
amplican_plot_cuts <- function(alignments, config, id) {

  archRanges <- alignments[alignments$seqnames %in% id & alignments$type == "deletion",]
  archRanges <- archRanges[archRanges$cut == TRUE,]

  if (length(unique(archRanges$strand)) == 2) { archRanges <- archRanges[archRanges$strand == "+",] }

  amplicon <- toString(config[which(config$ID == id[1]), "Amplicon"])
  ampl_len <- nchar(amplicon)

  if (dim(archRanges)[1] == 0) {return(print("No cuts to plot."))}
  archRanges <- stats::aggregate(cbind(count, frequency, cut) ~ strand + start + end + seqnames, archRanges, sum)

  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id[1]), "Direction"] != 1) {
    if (length(box) >= 1) {
      seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
    } else {
      seq(1, ampl_len, 4)
    }
  } else {
    if (length(box) >= 1) {
      rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
    } else {
      rev(seq(1, ampl_len, 4))
    }
  }
  xbreaks <- seq(1, ampl_len, 4)

  frPrimer <- toString(config[which(config$ID == id[1]), "Forward_Primer"])
  frPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon), toupper(frPrimer)))
  rwPrimer <- toString(config[which(config$ID == id[1]), "Reverse_Primer"])
  rwPrimer <- stats::na.omit(stringr::str_locate(toupper(amplicon),
                                                 toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer)))))))
  primers <- c(frPrimer, rwPrimer)
  if (length(frPrimer) == 0) {frPrimer <- c(1, 1)}
  if (length(rwPrimer) == 0) {rwPrimer <- c(ampl_len, ampl_len)}

  #FILTER EVENTS
  # forward before primer after amplicon length
  archRanges <- archRanges[!(archRanges$strand == "+" & archRanges$start <= frPrimer[1]),]
  archRanges <- archRanges[!(archRanges$strand == "+" & archRanges$end >= ampl_len),]
  # reverse before primer and after amplicon start
  archRanges <- archRanges[!(archRanges$strand == "-" & archRanges$start >= rwPrimer[1]),]
  archRanges <- archRanges[!(archRanges$strand == "-" & archRanges$start == 1),]
  if (dim(archRanges)[1] == 0) {return(print("No cuts to plot."))}

  frequency <- cut <- NULL
  arch_plot_top <- ggplot2::ggplot() + ggbio::geom_arch(data = archRanges,
                                                        ggplot2::aes(height = frequency,
                                                                     alpha = frequency,
                                                                     order = frequency,
                                                                     size = frequency,
                                                                     x = start,
                                                                     xend = end)) +
    scale_color_manual(values=c("dogerblue2")) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank()) +
    ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue")

  arch_plot_bot <- ggplot2::ggplot() + ggbio::geom_arch(data = archRanges,
                                                       ggplot2::aes(height = frequency,
                                                                    alpha = frequency,
                                                                    order = frequency,
                                                                    size = frequency,
                                                                    x = start,
                                                                    xend = end)) +
    facet_grid(seqnames ~ .) +
    ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.margin = unit(0, "cm"),
                   legend.position="none", legend.margin = unit(0, "cm"),
                   panel.grid=element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   strip.text.y = element_text(angle = 0)) +
    #ggplot2::ylab("Frequency [%]") +
    ggplot2::scale_x_continuous(labels = xlabels, breaks = xbreaks) +
    #ggplot2::geom_vline(xintercept = c(start(box), end(box)), linetype = "longdash", colour = "black") +
    ggplot2::geom_vline(xintercept = primers, linetype = "dotdash", colour = "blue")

  amplicon_colors <- c("green3", "red", "gold2", "blue", "green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T", "a", "c", "g", "t")
  position <- nucleotide <- upper <- count <- NULL
  amp_plot <- ggplot2::ggplot(ampl_df, ggplot2::aes(x = position, label = nucleotide, colour = upper, y = count)) +
    ggplot2::geom_text(size = I(4)) + ggplot2::ylim(0.7, 1.1) + ggbio::xlim(1, ampl_len) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   legend.position="none", legend.margin = unit(0, "cm"), line = element_blank(),
                   rect = element_blank(), axis.title.x = element_blank(), axis.text = element_blank(),
                   axis.title.y = element_blank()) +
    ggplot2::scale_colour_manual(drop = F, values = amplicon_colors)

  p <- ggbio::tracks(arch_plot_top, amp_plot, arch_plot_bot,
                     heights = c(0.3, 0.03, 0.53), padding = -1, xlim = 1:ampl_len,
                     xlab = "Nucleotide Position Relative to PAM")
  return(p)
}
