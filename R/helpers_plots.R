#' This function plots mutations in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignmentsGRanges (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string) Name of the ID column from config file.
#' @return (mutation plot) ggplot2 object of mutation plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual geom_bar
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom ggbio tracks xlim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom stats na.omit aggregate
#' @export
#'
amplican_mutations_plot <- function(alignmentsGRanges, config, id) {

  idRanges <- alignmentsGRanges[alignmentsGRanges$seqnames == id,]
  idRanges <- idRanges[idRanges$type == "mismatch",]

  amplicon <- toString(config[which(config$ID == id), "Amplicon"])
  ampl_len <- nchar(amplicon)
  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id), "Strand"] != 1) {
    seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
  } else {
    rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
  }
  xbreaks <- seq(1, nchar(amplicon), 4)
  box <- box + 5

  frPrimer <- toString(config[which(config$ID == id), "Forward_Primer"])
  frPrimer <- stringr::str_locate(toupper(amplicon), toupper(frPrimer))
  rwPrimer <- toString(config[which(config$ID == id), "Reverse_Primer"])
  rwPrimer <- stringr::str_locate(toupper(amplicon), toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer))))))
  primers <- stats::na.omit(c(frPrimer, rwPrimer))

  amplicon_colors <- c("green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T")

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

  mut_fr <- ggplot() +
    geom_bar(data = freqAgr[freqAgr$strand == "+",],
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
    geom_bar(data = freqAgr[freqAgr$strand == "-",],
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


#' This function plots deletions in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignmentsGRanges (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string) Name of the ID column from config file.
#' @return (deletions plot) ggplot2 object of mutation plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual
#' scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks geom_arch xlim
#' @importFrom stats na.omit
#' @export
#'
amplican_deletions_plot <- function(alignmentsGRanges, config, id) {

  archRanges <- alignmentsGRanges[alignmentsGRanges$seqnames == id,]
  archRanges <- archRanges[archRanges$type == "deletion",]
  archRanges <- stats::aggregate(cbind(count, frequency, cut) ~ strand + start + end, archRanges, sum)

  amplicon <- toString(config[which(config$ID == id), "Amplicon"])
  ampl_len <- nchar(amplicon)
  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id), "Strand"] != 1) {
    seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
  } else {
    rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
  }
  xbreaks <- seq(1, nchar(amplicon), 4)
  box <- box + 5

  frPrimer <- toString(config[which(config$ID == id), "Forward_Primer"])
  frPrimer <- stringr::str_locate(toupper(amplicon), toupper(frPrimer))
  rwPrimer <- toString(config[which(config$ID == id), "Reverse_Primer"])
  rwPrimer <- stringr::str_locate(toupper(amplicon), toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer))))))
  primers <- stats::na.omit(c(frPrimer, rwPrimer))

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

  amplicon_colors <- c("green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T")
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


#' This function plots insertions in relation to the amplicon. Top plot is for the forward reads, middle one shows
#' amplicon sequence, and bottom plot is for reverse reads.
#'
#' @param alignmentsGRanges (GRanges object) Loaded alignment information from alignments.csv file.
#' @param config (data.frame) Loaded table from config_summary.csv file.
#' @param id (string) Name of the ID column from config file.
#' @return (insertions plot) ggplot2 object of mutation plot
#' @import GenomicRanges
#' @importFrom ggplot2 ggplot aes theme_bw theme geom_label ggtitle scale_colour_manual
#' geom_polygon scale_fill_manual scale_x_continuous geom_vline scale_y_reverse element_blank unit geom_text ylab ylim
#' @importFrom seqinr s2c comp
#' @importFrom stringr str_locate
#' @importFrom ggbio tracks xlim
#' @importFrom stats na.omit aggregate
#' @export
#'
amplican_insertions_plot <- function(alignmentsGRanges, config, id) {
  idRanges <- alignmentsGRanges[alignmentsGRanges$seqnames == id,]
  idRanges <- idRanges[idRanges$type == "insertion",]

  amplicon <- toString(config[which(config$ID == id), "Amplicon"])
  ampl_len <- nchar(amplicon)
  selcolumns <- c("position", "nucleotide", "upper", "count")
  ampl_df <- data.frame(seq(1, ampl_len), seqinr::s2c(amplicon), seqinr::s2c(toupper(amplicon)), 1)
  names(ampl_df) <- selcolumns

  frPrimer <- toString(config[which(config$ID == id), "Forward_Primer"])
  frPrimer <- stringr::str_locate(toupper(amplicon), toupper(frPrimer))
  rwPrimer <- toString(config[which(config$ID == id), "Reverse_Primer"])
  rwPrimer <- stringr::str_locate(toupper(amplicon), toupper(seqinr::c2s(seqinr::comp(rev(seqinr::s2c(rwPrimer))))))
  primers <- stats::na.omit(c(frPrimer, rwPrimer))

  #reduce mismatches
  idRangesReduced <- stats::aggregate(frequency ~ strand + start + end, idRanges, sum)
  idRangesFr <- idRangesReduced[idRangesReduced$strand == "+",]
  idRangesFrgroup = rep(1:dim(idRangesFr)[1], each = 3)
  idRangesFrFrequency = rep(idRangesFr$frequency, each = 3)
  idRangesFrFrequency[c(T, F, F)] <- 0
  idRangesFrX <- as.vector(rbind(idRangesFr$start, idRangesFr$start, idRangesFr$end))
  traingleFr <- data.frame(frequency = idRangesFrFrequency, position = idRangesFrX, group = idRangesFrgroup)

  idRangesRe <- idRangesReduced[idRangesReduced$strand == "-",]
  idRangesRegroup = rep(1:dim(idRangesRe)[1], each = 3)
  idRangesReFrequency = rep(idRangesRe$frequency, each = 3)
  idRangesReFrequency[c(T, F, F)] <- 0
  idRangesReX <- as.vector(rbind(idRangesRe$start, idRangesRe$start, idRangesRe$end))
  traingleRe <- data.frame(frequency = idRangesReFrequency, position = idRangesReX, group = idRangesRegroup)
  ampl_len <- max(idRangesReX, idRangesFrX) #new x scale max

  box <- upperGroups(amplicon)
  xlabels <- if (config[which(config$ID == id), "Strand"] != 1) {
    seq(-start(box[1]) + 1, ampl_len-start(box[1]), 4)
  } else {
    rev(seq(end(box[1]) - ampl_len + 1, end(box[1]), 4))
  }
  xbreaks <- seq(1, ampl_len, 4)
  box <- box + 5

  amplicon_colors <- c("green3", "red", "gold2", "blue")
  names(amplicon_colors) <- c("A", "C", "G", "T")
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
