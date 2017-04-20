#' Pretty print forward and reverse reads aligned to each other.
#'
#' Usefull and needed for barcode reports.
#'
#' @param forward (character or vector of characters) Forward reads.
#' @param reverse (character or vector of characters) Will be reverse
#' complemented before alignment.
#' @return Vector with alignments ready to be printed.
#' @import Biostrings
#' @importFrom utils capture.output
#' @include helpers_general.R
#' @export
#'
amplican_print_reads <- function(forward, reverse) {

  alignForwardReverse <- Biostrings::pairwiseAlignment(forward,
                                                       revComp(reverse))
  wPA <- capture.output(Biostrings::writePairwiseAlignments(alignForwardReverse))
  wPA <- wPA[!grepl("#", wPA)]

  wPA
}


#' Function to calculate figure height
#' based on number of elements to plot.
#'
#' @param x (numeric) number of elements to fit onto height axis
#' @return (numeric)
#'
plot_height <- function(x) {
  if (x < 10) {
    x + 1
  } else {
    ceiling(9 * log(x))
  }
}

#' Writes head of the .Rmd report files.
#'
#' @param title (string) title to include in header
#' @return string of the header for rmd file
#'
write_head <- function(title) {
  return(c("---",
           paste0("title: '", title, "'"),
           "author: amplican",
           "date: '`r format(Sys.time(), \"%d %B %Y\")`'",
           "output:",
           "  html_document:",
           "    toc: true",
           "    theme: paper",
           "    toc_float: true",
           "    number_sections: true",
           "---\n"))
}

#' Writes explanation to the .Rmd report files.
#'
#' @return string of the explanation for rmd file
#'
write_explanation <- function() {
  return(c("**Read distribution plot** - plot shows number of reads assigned during read grouping  ",
           "**PRIMER DIMER filter** - plot shows percentage of assigned reads that have been recognized as PRIMER DIMERS  ",
           "**Cutting rates** - plot gives overview of percentage of reads (not filtered as PRIMER DIMER) that have cut  ",
           "**Frameshift** - plot shows what percentage of reads that have frameshift  ",
           paste0("**Read heterogeneity plot** - shows what is the share of each of the unique reads in total count of all reads.",
                  " The more blue each row, the less heterogeneity in the reads, more green means reads don't repeat often and are unique  ")))
}


#' Writes alignment plots - aggregate on the amplicon: deletion, insertion and mismatch.
#'
#' @param ID (string or vector of strings)
#' @param title (string) title for the alignment plots
#' @param cut_buffer (numeric) cut detection box will be extended left and right by box_buffer
#' @return string to include in rmd file
#'
write_alignment_plots <- function(ID, title, cut_buffer = 5) {

  id_pass <- paste0("c('", paste0(ID, collapse = "','"), "')")
  return(c(paste0("## ", title, "  \n"),
           "### Deletions  \n",
           paste0("```{r deletions ", title, ", echo=F, fig.height=14, fig.width=30, message=F, warning=F}"),
           paste0("amplican_plot_deletions(alignments, config, ", id_pass, ", ", cut_buffer, ")"),
           "```  \n",
           "### Insertions  \n",
           paste0("```{r insertions ", title, ", echo=F, fig.height=14, fig.width=30, message=F, warning=F}"),
           paste0("amplican_plot_insertions(alignments, config, ", id_pass, ", ", cut_buffer, ")"),
           "```  \n",
           "### Mismatches  \n",
           paste0("```{r mismatches ", title, ", echo=F, fig.height=14, fig.width=30, message=F, warning=F}"),
           paste0("amplican_plot_mismatches(alignments, config, ", id_pass, ", ", cut_buffer, ")"),
           "```  \n"))
}


#' Writes unassigned reads
#'
#' @param barcode (string) CSV name of the barcode unassigned reads file.
#' @param top (numeric) How many from the most abundant unassigned reads write.
#' @return string to include in rmd file
#'
write_unassigned_reads <- function(barcode, top = 5) {
  pure_bname <- gsub("_unassigned_reads.csv", "", barcode)
  return(c(paste0("## ", pure_bname, "  \n"),
           paste0("```{r plot unassigned reads ", pure_bname, ", echo=FALSE, message=F, warning=FALSE, comment = ''}"),
           paste0("unassigned_reads <- read.csv(file.path(results_folder, 'alignments', 'unassigned_sequences', '", barcode, "'), stringsAsFactors = FALSE)"),
           "unassigned_reads <- unassigned_reads[order(unassigned_reads$BarcodeFrequency, decreasing = TRUE), ]",
           "if (dim(unassigned_reads)[1] == 0) {",
           '  "No unassigned reads in this barcode. That\'s great!"',
           "} else {",
           paste0("  topN <- if (dim(unassigned_reads)[1] < ", top, ") dim(unassigned_reads)[1] else ", top),
           "  knitr::kable(data.frame(Forward = paste0('P', 1:topN),",
           "                          Reverse = paste0('S', 1:topN),",
           "                          Counts = unassigned_reads[1:topN, 'Total'],",
           "                          Frequency = unassigned_reads[1:topN, 'BarcodeFrequency']))",
           "",
           "  knitr::asis_output(cat(amplican_print_reads(unassigned_reads[1:topN, 'Forward'],",
           "                                              unassigned_reads[1:topN, 'Reverse']), sep = '\n'))",
           "}",
           "```\n"))
}


#' Make contents for .Rmd file aggregated on ID.
#'
#' @param results_folder (string) path to results
#' @param cut_buffer (numeric)
#' @return c() collection of strings to put into rmd file
#' @importFrom utils read.csv
#'
make_id_rmd <- function(results_folder, cut_buffer = 5) {

  config <- utils::read.csv(file.path(results_folder, "config_summary.csv"),
                            stringsAsFactors = FALSE)
  id_alignments_plots <- unlist(lapply(config$ID,
                                       function(x)
                                         write_alignment_plots(x, x, cut_buffer)))
  height <- plot_height(length(unique(config$ID)))

  return(c(write_head("Report breakdown by ID"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(file.path(results_folder, 'alignments_events.csv'), stringsAsFactors = FALSE)",
           "config <- read.csv(file.path(results_folder, 'config_summary.csv'), stringsAsFactors = FALSE)",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "**Deletions plot** - shows summary of deletions detected after alignments with distinction",
           "for forward (top plot) and reverse (bottom) reads, blue dotted lines represent primers as black",
           "dotted line represents cut site box, for deletions overlapping with cut site box there is distinction",
           "in color  ",
           "**Mismatches plot** - shows summary of mismatches detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, mismatches are colored in the same manner as amplicon  ",
           "**Insertions plot** - shows summary of insertions detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, insertion is shown as right-angled triangle where size of",
           "the insertion corresponds to the width of the triangle, size and transparency of triangle reflect on",
           "the frequency of the insertion\n",
           "***\n",
           "# ID Summary  \n",
           "***\n",
           "## Read distribution  \n",
           paste0("```{r plot total reads, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "ggplot(data = config, aes(x = ID, y = log10(Reads + 1), order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Number of reads + 1, log10 scaled')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip() +",
           "  geom_text(aes(x = ID, y = log10(Reads + 1), label = Reads), hjust = -1)",
           "```\n",
           "## PRIMER DIMER filter  \n",
           paste0("```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = PD_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip() +",
           "  geom_text(aes(x = ID, y = PD_percentage, label = PRIMER_DIMER), hjust = -1)",
           "```  \n",
           "## Cut rates  \n",
           paste0("```{r plot cut percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = cut_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip() +",
           "  geom_text(aes(x = ID, y = cut_percentage, label = Cut), hjust = -1)",
           "```  \n",
           "## Frameshift  \n",
           paste0("```{r plot frame shift percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = frameshift_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip() +",
           "  geom_text(aes(x = ID, y = frameshift_percentage, label = Frameshift), hjust = -1)",
           "```  \n",
           "## Heterogeneity of reads  \n",
           paste0("```{r plot read domination, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]  \n",
           "cumsum_list <- tapply(uniqueReadsByID$count, uniqueReadsByID$seqnames, FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$seqnames))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)\n",
           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)\n",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),]  \n",
           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "uniqueReadsByID$read_number <- as.vector(unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes$read_id, by = 1)))  \n",
           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, uniqueReadsByID$seqnames, FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$seqnames))]",
           "#config[match(howManyTimes$seqnames, config$ID),]",
           "toDivide <- rep(ids_with_reads, times = howManyTimes$read_id)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",
           "# divide into bins for colour",
           "uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal,",
           "                            c(0, 5, seq(10, 100, 10)))",
           "# reduce number of reads in 0-5 group - faster plots without artifacts",
           "uniqueReadsByID <- aggregate(read_share_percentage_normal ~ seqnames + bins,",
           "                             data = uniqueReadsByID, sum)",

           "colorPalette <- colorRampPalette(c('#009E73', '#0072B2'))(length(levels(uniqueReadsByID$bins)))",
           "names(colorPalette) <- levels(uniqueReadsByID$bins)",
           "ggplot(data = uniqueReadsByID, aes(x = seqnames, y = read_share_percentage_normal, fill = bins)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  theme(axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold'),",
           "        legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of shares') +",
           "  xlab('ID') +",
           "  scale_fill_manual(values = colorPalette) +",
           "  coord_flip()",
           "```  \n",
           "***\n",
           "# Alignments plots  \n",
           "***\n",
           id_alignments_plots))
}


#' Make contents for .Rmd file aggregated on amplicon.
#'
#' @param results_folder (string) path to results
#' @param cut_buffer (numeric)
#' @return c() collection of strings to put into rmd file
#' @importFrom utils read.csv
#'
make_amplicon_rmd <- function(results_folder, cut_buffer = 5) {

  config <- utils::read.csv(file.path(results_folder, "config_summary.csv"),
                            stringsAsFactors = FALSE)
  config$AmpliconUpper <- toupper(config$Amplicon)
  uniqueAmlicons <- unique(config$AmpliconUpper)
  height <- plot_height(length(uniqueAmlicons))

  ampl_alignments_plots <- c(mapply(function(x, i) {

    ID <- config$ID[config$AmpliconUpper == x]
    id_pass <- paste0("c('", paste0(ID, collapse = "','"), "')")

    c(write_alignment_plots(ID, paste0("Group ", i), cut_buffer),
      "### Cuts  \n",
      paste0("```{r cuts ", paste0("Group ", i), ", echo=F, fig.height=14, fig.width=30, message=F, warning=F}"),
      paste0("amplican_plot_cuts(alignments, config, ", id_pass, ")"),
      "```  \n")
  }, uniqueAmlicons, 1:length(uniqueAmlicons)))

  return(c(write_head("Report breakdown by Amplicon"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(file.path(results_folder, 'alignments_events.csv'), stringsAsFactors = FALSE)",
           "config <- read.csv(file.path(results_folder, 'config_summary.csv'), stringsAsFactors = FALSE)",
           "config$AmpliconUpper <- toupper(config$Amplicon)",
           "uniqueAmlicons <- unique(config$AmpliconUpper)",
           "labels <- sapply(uniqueAmlicons, function(x) {toString(config$ID[config$AmpliconUpper == x])})",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "**Deletions plot** - shows summary of deletions detected after alignments with distinction",
           "for forward (top plot) and reverse (bottom) reads, blue dotted lines represent primers as black",
           "dotted line represents cut site box, for deletions overlapping with cut site box there is distinction",
           "in color  ",
           "**Mismatches plot** - shows summary of mismatches detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, mismatches are colored in the same manner as amplicon  ",
           "**Insertions plot** - shows summary of insertions detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, insertion is shown as right-angled triangle where size of",
           "the insertion corresponds to the width of the triangle, size and transparency of triangle reflect on",
           "the frequency of the insertion  ",
           "**Cuts plot** - shows summary of cuts in the amplicon (insertions that are overlapping with ",
           "specified box of uppercase letters in the amplicon, and are supported by both forward and reverse reads)\n",
           "***\n",
           "# Amplicon Summary  \n",
           "***\n",

           "## Amplicon groups\n",

           paste0("```{r amplicon groups, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "library(knitr)",
           "ampliconDF <- data.frame(group = 1:length(uniqueAmlicons), IDs = unname(labels))",
           "kable(ampliconDF)",
           "ampliconDF$Amplicon <- uniqueAmlicons #for sorting",
           "config$group <- match(config$AmpliconUpper, ampliconDF$Amplicon)",
           "```\n",


           "## Read distribution  \n",

           paste0("```{r plot total reads, echo=FALSE, fig.height=",
                  height,", fig.width=14, message=F, warning=FALSE}"),
           "ampliconTable <- aggregate(cbind(Reads, PRIMER_DIMER, Cut, Frameshift) ~ group, data = config, sum)\n",

           "ggplot(data = ampliconTable, aes(x = group, y = log10(Reads + 1)), order = group) +",
           "  geom_bar(stat = 'identity') +",
           "  ylab('Number of reads + 1, log10 scaled')  +",
           "  xlab('Amplicon Group') +",
           "  scale_x_continuous(breaks = 1:length(uniqueAmlicons)) +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip() +",
           "  geom_text(aes(x = group, y = log10(Reads + 1), label = Reads), hjust = -1)",
           "```\n",

           "## PRIMER DIMER filter  \n",

           paste0("```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}\n"),

           "ampliconTable$PD_percentage <- ampliconTable$PRIMER_DIMER * 100/ampliconTable$Reads",
           "ampliconTable$PD_percentage[is.nan(ampliconTable$PD_percentage)] <- 0 \n",

           "ggplot(data = ampliconTable, aes(x = group, y = PD_percentage, order = group)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  xlab('Amplicon Group') +",
           "  scale_x_continuous(breaks = 1:length(uniqueAmlicons)) +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = group, y = PD_percentage, label = PRIMER_DIMER), hjust = -1)",
           "```  \n",

           "## Cut rates  \n",

           paste0("```{r plot cut rate, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "ampliconTable$Reads_noPD <- ampliconTable$Reads - ampliconTable$PRIMER_DIMER",
           "ampliconTable$cut_percentage <- ampliconTable$Cut * 100/ampliconTable$Reads_noPD",
           "ampliconTable$cut_percentage[is.nan(ampliconTable$cut_percentage)] <- 0  \n",

           "ggplot(data = ampliconTable, aes(x = group, y = cut_percentage, order = group)) +",
           "  geom_bar(stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  xlab('Amplicon Group') +",
           "  scale_x_continuous(breaks = 1:length(uniqueAmlicons)) +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0,100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = group, y = cut_percentage, label = Cut), hjust = -1)",
           "```  \n",

           "## Frameshift  \n",

           paste0("```{r plot frame shift percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "ampliconTable$frameshift_percentage <- ampliconTable$Frameshift * 100/ampliconTable$Reads_noPD",
           "ampliconTable$frameshift_percentage[is.nan(ampliconTable$frameshift_percentage)] <- 0 \n",

           "ggplot(data = ampliconTable, aes(x = group, y = frameshift_percentage, order = group)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  xlab('Amplicon Group') +",
           "  scale_x_continuous(breaks = 1:length(uniqueAmlicons)) +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = group, y = frameshift_percentage, label = Frameshift), hjust = -1)",
           "```  \n",

           "## Heterogeneity of reads  \n",

           paste0("```{r plot read heterogeneity, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$group <- rep(config$group[match(howManyTimes$seqnames, unique(config$ID))], ",
           "                             times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$group, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$group), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$group))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)\n",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$group))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$group))]",
           "dimHMT <- if (is.null(dim(howManyTimes)[1])) { length(howManyTimes) } else { dim(howManyTimes)[1] }",
           "uniqueReadsByID$read_number <- as.vector(unlist(seq2(from = rep(1, dimHMT), to = howManyTimes, by = 1)))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$group), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$group))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",
           "# divide into bins for colour",
           "uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal, ",
           "                            c(0, 5, seq(10, 100, 10)))",
           "# reduce number of reads in 0-5 group - faster plots without artifacts",
           "uniqueReadsByID <- aggregate(read_share_percentage_normal ~ group + bins,",
           "                             data = uniqueReadsByID, sum)",

           "colorPalette <- colorRampPalette(c('#009E73', '#0072B2'))(length(levels(uniqueReadsByID$bins)))",
           "names(colorPalette) <- levels(uniqueReadsByID$bins)",
           "ggplot(data = uniqueReadsByID, aes(x = as.factor(group), ",
           "                                   y = read_share_percentage_normal, ",
           "                                   fill = bins), ",
           "       order = group) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('Amplicon Group') +",
           "  scale_fill_manual(values = colorPalette) +",
           "  coord_flip()\n",

           "# fix frequency to be relative to amplicon - not ID",
           "alignments$group <- config$group[match(alignments$seqnames, config$ID)]",
           "alignments$frequency <- alignments$count / ampliconTable$Reads[alignments$group]",
           "``` \n",
           "# Alignments plots  \n",
           "***\n",
           ampl_alignments_plots))
}


#' Make contents for .Rmd file aggregated on Barcode
#'
#' @param results_folder path to results as string
#' @return c() collection of strings to put into rmd file
#'
make_barcode_rmd <- function(results_folder) {

  config <- read.csv(file.path(results_folder, 'config_summary.csv'),
                     stringsAsFactors = FALSE)
  height <- plot_height(length(unique(config$Barcode)))
  ub <- list.files(file.path(results_folder,
                             "alignments", "unassigned_sequences"))
  ub_reads <- unlist(lapply(ub, function(x) write_unassigned_reads(x)))

  return(c(write_head("Report breakdown by Barcode"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(file.path(results_folder, 'alignments_events.csv'), stringsAsFactors = FALSE)",
           "config <- read.csv(file.path(results_folder, 'config_summary.csv'), stringsAsFactors = FALSE)",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "**Top unassigned reads** - take a look at the alignment of most abundant ",
           "forward and reverse complemented reverse reads for each barcode, ",
           "if you find that there is many unassigned reads you can ivestigate here.  ",
           "\n***\n",
           "# Barcode Summary  \n",
           "***\n",
           "## Read distribution  \n",

           paste0("```{r plot total reads, echo=FALSE, fig.height=",
                  plot_height(length(unique(config$ID))),
                  ", fig.width=14, message=F, warning=FALSE}"),
           "ggplot(data = config, aes(x = Barcode, y = log10(Reads + 1), order = Barcode, fill = ID)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  ylab('Number of reads + 1, log10 scaled')  +",
           "  theme(legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold'),",
           "        legend.title = element_blank()) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           paste0("```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "barcodeTable <- aggregate(cbind(Reads, PRIMER_DIMER, Cut, Frameshift) ~ Barcode, data = config, sum)",
           "barcodeTable$PD_percentage <- barcodeTable$PRIMER_DIMER * 100/barcodeTable$Reads",
           "barcodeTable$PD_percentage[is.nan(barcodeTable$PD_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = PD_percentage, order = Barcode)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = Barcode, y = PD_percentage, label = PRIMER_DIMER), hjust = -1)",
           "``` \n",

           "## Cut rates  \n",

           paste0("```{r plot cut percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "barcodeTable$Reads_noPD <- barcodeTable$Reads - barcodeTable$PRIMER_DIMER",
           "barcodeTable$cut_percentage <- barcodeTable$Cut * 100/barcodeTable$Reads_noPD",
           "barcodeTable$cut_percentage[is.nan(barcodeTable$cut_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = cut_percentage, order = Barcode)) +",
           "  geom_bar(stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0,100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = Barcode, y = cut_percentage, label = Cut), hjust = -1)",
           "``` \n",

           "## Frameshift  \n",

           paste0("```{r plot frame shift percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "barcodeTable$frameshift_percentage <- barcodeTable$Frameshift * 100/barcodeTable$Reads_noPD",
           "barcodeTable$frameshift_percentage[is.nan(barcodeTable$frameshift_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = frameshift_percentage, order = Barcode)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip() +",
           "  geom_text(aes(x = Barcode, y = frameshift_percentage, label = Frameshift), hjust = -1)",
           "``` \n",

           "## Heterogeneity of reads  \n",

           paste0("```{r plot read heterogeneity, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$barcode <- rep(config$Barcode[match(howManyTimes$seqnames, unique(config$ID))], ",
           "                               times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$barcode, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$barcode), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$barcode))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)\n",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$barcode))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$barcode))]",
           "dimHMT <- if (is.null(dim(howManyTimes)[1])) { length(howManyTimes) } else { dim(howManyTimes)[1] }",
           "uniqueReadsByID$read_number <- as.vector(unlist(seq2(from = rep(1, dimHMT), to = howManyTimes, by = 1)))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$barcode), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$barcode))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "# divide into bins for colour",
           "uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal,",
           "                            c(0, 5, seq(10, 100, 10)))",
           "# reduce number of reads in 0-5 group - faster plots without artifacts",
           "uniqueReadsByID <- aggregate(read_share_percentage_normal ~ barcode + bins,",
           "                             data = uniqueReadsByID, sum)",

           "colorPalette <- colorRampPalette(c('#009E73', '#0072B2'))(length(levels(uniqueReadsByID$bins)))",
           "names(colorPalette) <- levels(uniqueReadsByID$bins)",

           "ggplot(data = uniqueReadsByID, aes(x = as.factor(barcode), y = read_share_percentage_normal, fill = bins)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  scale_fill_manual(values = colorPalette) +",
           "  xlab('Barcode') +",
           "  coord_flip()",
           "``` \n",
           "***\n",
           "# Top unassigned reads  \n",
           "***\n",
           ub_reads))
}


#' Make contents for .Rmd file aggregated on Type of Experiment
#'
#' @param results_folder path to results as string
#' @return c() collection of strings to put into rmd file
#'
make_group_rmd <- function(results_folder) {

  config <- read.csv(file.path(results_folder, 'config_summary.csv'),
                     stringsAsFactors = FALSE)
  height <- plot_height(length(unique(config$Group)))

  return(c(write_head("Report breakdown by Group"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(file.path(results_folder, 'alignments_events.csv'), stringsAsFactors = FALSE)",
           "config <- read.csv(file.path(results_folder, 'config_summary.csv'), stringsAsFactors = FALSE)",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "\n",
           "***\n",
           "# Group Summary  \n",
           "***\n",

           "## Read distribution  \n",

           paste0("```{r plot total reads, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "ggplot(data = config, aes(x = Group, y = log10(Reads + 1), order = Group, fill = Group)) +",
           "  geom_boxplot() +",
           "  ylab('Number of reads + 1, log10 scaled')  +",
           "  xlab('Group') + ",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           paste0("```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0\n",

           "ggplot(data = config, aes(x = Group, y = PD_percentage, order = Group, fill = Group)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'none') +",
           "  ylim(0, 100) +",
           "  coord_flip()",
           "``` \n",

           "## Cut rates  \n",

           paste0("```{r plot cut percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Group, y = cut_percentage, order = Group, fill = Group)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0,100) +",
           "  coord_flip()",
           "``` \n",

           "## Frameshift  \n",

           paste0("```{r plot frame shift percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Group, y = frameshift_percentage, order = Group, fill = Group)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0, 100) +",
           "  coord_flip()",
           "``` \n",

           "## Heterogeneity of reads  \n",

           paste0("```{r plot read heterogeneity, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$group <- rep(config$Group[match(howManyTimes$seqnames, unique(config$ID))],",
           "times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$group, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$group), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$group))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$group))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$group))]",
           "dimHMT <- if (is.null(dim(howManyTimes)[1])) { length(howManyTimes) } else { dim(howManyTimes)[1] }",
           "uniqueReadsByID$read_number <- as.vector(unlist(seq2(from = rep(1, dimHMT), to = howManyTimes, by = 1)))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$group), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$group))]",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "# divide into bins for colour",
           "uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal,",
           "                            c(0, 5, seq(10, 100, 10)))",
           "# reduce number of reads in 0-5 group - faster plots without artifacts",
           "uniqueReadsByID <- aggregate(read_share_percentage_normal ~ group + bins,",
           "                             data = uniqueReadsByID, sum)",

           "colorPalette <- colorRampPalette(c('#009E73', '#0072B2'))(length(levels(uniqueReadsByID$bins)))",
           "names(colorPalette) <- levels(uniqueReadsByID$bins)",

           "ggplot(data = uniqueReadsByID, aes(x = group, y = read_share_percentage_normal, fill = bins)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('Group') +",
           "  scale_fill_manual(values = colorPalette) +",
           "  coord_flip() +",
           "  coord_flip()",
           "``` \n"))
}


#' Make contents for .Rmd file aggregated on the guideRNA
#'
#' @param results_folder path to results as string
#' @return c() collection of strings to put into rmd file
#'
make_guide_rmd <- function(results_folder) {

  config <- read.csv(file.path(results_folder, 'config_summary.csv'),
                     stringsAsFactors = FALSE)
  height <- plot_height(length(unique(config$guideRNA)))

  return(c(write_head("Report breakdown by guideRNA"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(file.path(results_folder, 'alignments_events.csv'), stringsAsFactors = FALSE)",
           "config <- read.csv(file.path(results_folder, 'config_summary.csv'), stringsAsFactors = FALSE)",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "\n",
           "***\n",
           "# guideRNA Summary  \n",
           "***\n",

           "## Read distribution  \n",

           paste0("```{r plot total reads, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "ggplot(data = config, aes(x = guideRNA, y = log10(Reads + 1), order = guideRNA, fill = guideRNA)) +",
           "  geom_boxplot() +",
           "  ylab('Number of reads + 1, log10 scaled')  +",
           "  xlab('guideRNA') + ",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           paste0("```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0\n",

           "ggplot(data = config, aes(x = guideRNA, y = PD_percentage, order = guideRNA, fill = guideRNA)) +",
           "  geom_boxplot() +",
           "  xlab('guideRNA') + ",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'none') +",
           "  ylim(0, 100) + ",
           "  coord_flip()",
           "``` \n",

           "## Cut rates  \n",

           paste0("```{r plot cut percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = guideRNA, y = cut_percentage, order = guideRNA, fill = guideRNA)) +",
           "  geom_boxplot() +",
           "  xlab('guideRNA') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0,100) + ",
           "  coord_flip()",
           "``` \n",

           "## Frameshift  ",

           paste0("```{r plot frame shift percentage, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = guideRNA, y = frameshift_percentage, order = guideRNA, fill = guideRNA)) +",
           "  geom_boxplot() +",
           "  xlab('guideRNA') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0, 100) + ",
           "  coord_flip()",
           "``` \n",

           "## Heterogeneity of reads  \n",

           paste0("```{r plot read heterogeneity, echo=FALSE, fig.height=",
                  height, ", fig.width=14, message=F, warning=FALSE}"),
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$guideRNA <- rep(config$guideRNA[match(howManyTimes$seqnames, unique(config$ID))], \n",
           "                                times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$guideRNA, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$guideRNA), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$guideRNA))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$guideRNA))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$guideRNA))]",
           "dimHMT <- if (is.null(dim(howManyTimes)[1])) { length(howManyTimes) } else { dim(howManyTimes)[1] }",
           "uniqueReadsByID$read_number <- as.vector(unlist(seq2(from = rep(1, dimHMT), to = howManyTimes, by = 1)))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$guideRNA), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$guideRNA))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "# divide into bins for colour",
           "uniqueReadsByID$bins <- cut(uniqueReadsByID$read_share_percentage_normal,",
           "                            c(0, 5, seq(10, 100, 10)))",
           "# reduce number of reads in 0-5 group - faster plots without artifacts",
           "uniqueReadsByID <- aggregate(read_share_percentage_normal ~ guideRNA + bins,",
           "                             data = uniqueReadsByID, sum)",

           "colorPalette <- colorRampPalette(c('#009E73', '#0072B2'))(length(levels(uniqueReadsByID$bins)))",
           "names(colorPalette) <- levels(uniqueReadsByID$bins)",

           "ggplot(data = uniqueReadsByID, aes(x = guideRNA, y = read_share_percentage_normal, fill = bins)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('guideRNA') +",
           "  scale_fill_manual(values = colorPalette) +",
           "  coord_flip()",
           "``` \n"))
}


#' Prepare summary report as .Rmd file.
#'
#' @param results_folder (string) Folder containing results from the \code{\link{amplicanAlign}} function,
#' do not change names of the files.
#' @param links (string) A string containing already processed links to other reports.
#' @return (string) contents to write to files
#'
make_summary_rmd <- function(results_folder, links) {

  config <- utils::read.csv(file.path(results_folder, "config_summary.csv"), stringsAsFactors = FALSE)
  height <- plot_height(length(unique(config$Barcode)))

  return(c(write_head("Summary Read Report"),
           links,
           "\n***\n",

           "# Explanation of variables\n",

           "***\n",

           "**experiment_count** - how many IDs belongs to this barcode  ",
           "**read_count** - how many reads belongs to this barcode  ",
           "**bad_base_quality** - how many reads had base quality worse than specified (default is 0)  ",
           "**bad_average_quality** - how many reads had average base quality worse than specified (default is 0)  ",
           "**bad_alphabet** - how many reads had alphabet with bases other than A, C, G, T  ",
           "**filtered_read_count** - how many reads were left after filtering  ",
           "**unique_reads** - how many reads (forward and reverse together) for this barcode is unique  ",
           "**unassigned_reads/assigned_reads** - how many reads have been not assigned/assigned to any of the
           experiments  \n",

           "***\n",

           "# Total reads\n",

           "***\n",

           "## Read Quality\n",

           "```{r echo = F}",
           "library(ggplot2)",
           paste0("results_folder = '", results_folder, "'"),
           "summaryDF <- read.csv(file.path(results_folder, 'barcode_reads_filters.csv'),",
           "                      stringsAsFactors = FALSE)",
           "total_reads <- sum(summaryDF$read_count)",
           "reads_names <- c('Good', 'Bad',",
           "                 'Bad Base Quality', 'Bad Average Read Quality',",
           "                 'Bad Read Aplhabet', 'Unassigned Reads')",
           "total_bad_reads <- sum(sum(summaryDF[, c(4:6, 9)]))",
           "readDF <- data.frame(levels = c(rep('ALL READS', 2),",
           "                                rep('BAD READS', 4)),",
           "                     type = reads_names,",
           "                     percentage = c(",
           "                      sum(summaryDF$assigned_reads)/total_reads,",
           "                      1 - sum(summaryDF$assigned_reads)/total_reads,",
           "                      sum(summaryDF$bad_base_quality)/total_bad_reads,",
           "                      sum(summaryDF$bad_average_quality)/total_bad_reads,",
           "                      sum(summaryDF$bad_alphabet)/total_bad_reads,",
           "                      sum(summaryDF$unassigned_reads)/total_bad_reads))",
           "readDF$type <- factor(readDF$type, levels = reads_names)\n",

           "library(ggthemes)",
           "ggplot(readDF,",
           "       aes(levels, percentage*100, fill = type)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  scale_fill_colorblind() +",
           "  ylab('[ % ]') +",
           "  xlab('') +",
           "  theme(legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        legend.title = element_blank())",
           "```\n",

           "## Mutants\n",

           "```{r echo = F}",
           "configDF <- read.csv(file.path(results_folder, 'config_summary.csv'),",
           "                     stringsAsFactors = FALSE)",
           "total_reads <- sum(configDF$Reads)",
           "reads_names <- c('Cut', 'No Cut', 'Frameshift', 'No Frameshift')",
           "readDF <- data.frame(levels = c(rep('CUT RATE', 2),",
           "                                rep('FRAMESHIFT', 2)),",
           "type = reads_names,",
           "percentage = c(sum(configDF$Cut)/total_reads,",
           "               1 - sum(configDF$Cut)/total_reads,",
           "               sum(configDF$Frameshift)/total_reads,",
           "               1 - sum(configDF$Frameshift)/total_reads))",
           "readDF$type <- factor(readDF$type, levels = reads_names)",

           "library(ggthemes)",
           "ggplot(readDF,",
           "       aes(levels, percentage*100, fill = type)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  scale_fill_colorblind() +",
           "  ylab('[ % ]') +",
           "  xlab('') +",
           "  theme(legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        legend.title = element_blank())",
           "```\n",

           "***\n",

           "# Summary Table\n",

           "***\n",

           "```{r echo = F}",
           "library(knitr)",
           "kable(summaryDF)",
           "```\n",

           "Table 1. Reads distributed for each barcode\n",

           "***\n",

           "# Reads by barcode\n",

           "***\n",

           paste0("```{r fig.width=8, fig.height=", height, ", echo = F}"),
           "library(reshape2)",
           "summaryDFmelt = melt(summaryDF,",
           "                     id.vars = c('barcode'),",
           "                     measure.vars = c('bad_base_quality',",
           "                                      'bad_average_quality',",
           "                                      'bad_alphabet',",
           "                                      'unassigned_reads',",
           "                                      'assigned_reads'))",
           "ggplot(data = summaryDFmelt,",
           "       aes(x = as.factor(barcode),",
           "           y = value,",
           "           fill = variable)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  ylab('number of reads') +",
           "  xlab('Barcode') +",
           "  theme(legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        legend.title = element_blank()) +",
           "  coord_flip()",
           "```\n",

           "Plot 1. Reads distribution for each barcode.\n",

           "***\n"))
}
