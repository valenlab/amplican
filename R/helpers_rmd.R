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
           "**Read heterogeneity plot** - shows what is the share of each of the unique reads in total count of all reads  "))
}


#' Writes alignment plots - aggregate on the amplicon: deltion, insertion and mismatch.
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
           "### Mutations  \n",
           paste0("```{r mutations ", title, ", echo=F, fig.height=14, fig.width=30, message=F, warning=F}"),
           paste0("amplican_plot_mutations(alignments, config, ", id_pass, ", ", cut_buffer, ")"),
           "```  \n"))
}


#' Make contents for .Rmd file aggregated on ID.
#'
#' @param results_folder (string) path to results
#' @param cut_buffer (numeric)
#' @return c() collection of strings to put into rmd file
#' @importFrom utils read.csv
#'
make_id_rmd <- function(results_folder, cut_buffer = 5) {

  config <- utils::read.csv(paste0(results_folder, "/config_summary.csv"), sep = "\t")
  id_alignments_plots <- unlist(lapply(config$ID, function(x) write_alignment_plots(x, x, cut_buffer)))

  return(c(write_head("Report breakdown by ID"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(paste0(results_folder, '/alignments.csv'), sep = '\\t')",
           "config <- read.csv(paste0(results_folder, '/config_summary.csv'), sep = '\\t')",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "**Deletions plot** - shows summary of deletions detected after alignments with distinction",
           "for forward (top plot) and reverse (bottom) reads, blue dotted lines represent primers as black",
           "dotted line represents cut site box, for deletions overlapping with cut site box there is distinction",
           "in colour  ",
           "**Mismatches plot** - shows summary of mismatches detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, mismatches are colored in the same manner as amplicon  ",
           "**Insertions plot** - shows summary of insertions detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, insertion is shown as right-angled triangle where size of",
           "the insertion coresponds to the width of the triangle, size and transparency of traingle reflect on",
           "the frequency of the insertion\n",
           "***\n",
           "# ID Summary  \n",
           "***\n",
           "## Read distribution  \n",
           "```{r plot total reads, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "ggplot(data = config, aes(x = ID, y = log10(Reads + 1), order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Number of reads log10 scaled')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```\n",
           "## PRIMER DIMER filter  \n",
           "```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = PD_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```  \n",
           "## Cut rates  \n",
           "```{r plot cut percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = cut_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```  \n",
           "## Frameshift  \n",
           "```{r plot frame shift percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",
           "ggplot(data = config, aes(x = ID, y = frameshift_percentage, order = ID)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have framshift')  +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```  \n",
           "## Heterogeneity of reads  \n",
           "```{r plot read domination, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]  \n",
           "cumsum_list <- tapply(uniqueReadsByID$count, uniqueReadsByID$seqnames, FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$seqnames))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)\n",
           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)\n",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),]  \n",
           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "uniqueReadsByID$read_number <- unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes$read_id, by = 1))  \n",
           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, uniqueReadsByID$seqnames, FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$seqnames))]",
           "#config[match(howManyTimes$seqnames, config$ID),]",
           "toDivide <- rep(ids_with_reads, times = howManyTimes$read_id)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",
           "ggplot(data = uniqueReadsByID, aes(x = seqnames, y = read_share_percentage_normal, fill = read_share_percentage_normal)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  theme(axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold'),",
           "        legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of shares') +",
           "  xlab('ID') +",
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

  config <- utils::read.csv(paste0(results_folder, "/config_summary.csv"), sep = "\t")
  config$AmpliconUpper <- toupper(config$Amplicon)
  uniqueAmlicons <- unique(config$AmpliconUpper)

  ampl_alignments_plots <- c(mapply(function(x, i) {
    write_alignment_plots(config$ID[config$AmpliconUpper == x], paste0("Group ", i), cut_buffer)
  }, uniqueAmlicons, 1:length(uniqueAmlicons)))

  return(c(write_head("Report breakdown by Amplicon"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           "results_folder = '/home/ai/removemelater'",
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(paste0(results_folder, '/alignments.csv'), sep = '\\t')",
           "config <- read.csv(paste0(results_folder, '/config_summary.csv'), sep = '\\t')",
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
           "in colour  ",
           "**Mismatches plot** - shows summary of mismatches detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, mismatches are colored in the same manner as amplicon  ",
           "**Insertions plot** - shows summary of insertions detected after alignments split by forward",
           "(top plot) and reverse (bottom) reads, insertion is shown as right-angled triangle where size of",
           "the insertion coresponds to the width of the triangle, size and transparency of traingle reflect on",
           "the frequency of the insertion\n",
           "***\n",
           "# Amplicon Summary  \n",
           "***\n",

           "## Amplicon groups\n",

           "```{r amplicon groups, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "library(knitr)",
           "ampliconDF <- data.frame(group = 1:length(uniqueAmlicons), IDs = unname(labels))",
           "kable(ampliconDF)",
           "ampliconDF$Amplicon <- uniqueAmlicons #for sorting",
           "config$group <- match(config$AmpliconUpper, ampliconDF$Amplicon)",
           "```\n",


           "## Read distribution  \n",

           "```{r plot total reads, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "ampliconTable <- aggregate(cbind(Reads, PRIMER_DIMER, Cut, Frameshift) ~ group, data = config, sum)\n",

           "ggplot(data = ampliconTable, aes(x = group, y = log10(Reads + 1)), order = group) +",
           "  geom_bar(stat = 'identity') +",
           "  ylab('Number of reads log10 scaled')  +",
           "  xlab('Amplicon Group') +",
           "  scale_x_continuous(breaks = 1:length(uniqueAmlicons)) +",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           "```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}\n",

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
           "  coord_flip()",
           "```  \n",

           "## Cut rates  \n",

           "```{r plot cut rate, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
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
           "  coord_flip()",
           "```  \n",

           "## Frameshift  \n",

           "```{r plot frame shift percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
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
           "  coord_flip()",
           "```  \n",

           "## Heterogeneity of reads  \n",

           "```{r plot read heterogeneity, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
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
           "uniqueReadsByID$read_number <- unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes, by = 1))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$group), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$group))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "ggplot(data = uniqueReadsByID, aes(x = as.factor(group), ",
           "                                   y = read_share_percentage_normal, ",
           "                                   fill = read_share_percentage_normal), ",
           "       order = group) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('Amplicon Group') +",
           "  coord_flip()\n",

           "#fix frequency to be relative to amplicon - not ID",
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
  return(c(write_head("Report breakdown by Barcode"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(paste0(results_folder, '/alignments.csv'), sep = '\\t')",
           "config <- read.csv(paste0(results_folder, '/config_summary.csv'), sep = '\\t')",
           "```\n",
           "***\n",
           "# Description  \n",
           "***\n",
           write_explanation(),
           "\n",
           "***\n",
           "# Barcode Summary  \n",
           "***\n",
           "## Read distribution  \n",

           "```{r plot total reads, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "ggplot(data = config, aes(x = Barcode, y = log10(Reads + 1), order = Barcode, fill = ID)) +",
           "  geom_bar(position='stack', stat='identity') +",
           "  ylab('Number of reads on log10 scale')  +",
           "  theme(legend.position = 'top',",
           "        legend.direction = 'horizontal',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold'),",
           "        legend.title = element_blank()) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           "```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "barcodeTable <- aggregate(cbind(Reads, PRIMER_DIMER, Cut, Frameshift) ~ Barcode, data = config, sum)",
           "barcodeTable$PD_percentage <- barcodeTable$PRIMER_DIMER * 100/barcodeTable$Reads",
           "barcodeTable$PD_percentage[is.nan(barcodeTable$PD_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = PD_percentage, order = Barcode)) +",
           "  geom_bar(stat='identity') +",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip()",
           "``` \n",

           "## Cut rates  \n",

           "```{r plot cut percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "barcodeTable$Reads_noPD <- barcodeTable$Reads - barcodeTable$PRIMER_DIMER",
           "barcodeTable$cut_percentage <- barcodeTable$Cut * 100/barcodeTable$Reads_noPD",
           "barcodeTable$cut_percentage[is.nan(barcodeTable$cut_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = cut_percentage, order = Barcode)) +",
           "  geom_bar(stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold')) +",
           "  ylim(0,100) +",
           "  coord_flip()",
           "``` \n",

           "## Frameshift  \n",

           "```{r plot frame shift percentage, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
           "barcodeTable$frameshift_percentage <- barcodeTable$Frameshift * 100/barcodeTable$Reads_noPD",
           "barcodeTable$frameshift_percentage[is.nan(barcodeTable$frameshift_percentage)] <- 0  \n",

           "ggplot(data = barcodeTable, aes(x = Barcode, y = frameshift_percentage, order = Barcode)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold')) +",
           "  ylim(0, 100) +",
           "  coord_flip()",
           "``` \n",

           "## Heterogeneity of reads  \n",

           "```{r plot read heterogeneity, echo=FALSE, fig.height=30, fig.width=14, message=F, warning=FALSE}",
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
           "uniqueReadsByID$read_number <- unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes, by = 1))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$barcode), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$barcode))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "ggplot(data = uniqueReadsByID, aes(x = barcode, y = read_share_percentage_normal, fill = read_share_percentage_normal)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('Barcode') +",
           "  coord_flip()",
           "``` \n" ))
}


#' Make contents for .Rmd file aggregated on Type of Experiment
#'
#' @param results_folder path to results as string
#' @return c() collection of strings to put into rmd file
#'
make_group_rmd <- function(results_folder) {
  return(c(write_head("Report breakdown by Group"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(paste0(results_folder, '/alignments.csv'), sep = '\\t')",
           "config <- read.csv(paste0(results_folder, '/config_summary.csv'), sep = '\\t')",
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

           "```{r plot total reads, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "ggplot(data = config, aes(x = Experiment_Type, y = log10(Reads + 1), order = Experiment_Type, fill = Experiment_Type)) +",
           "  geom_boxplot() +",
           "  ylab('Number of reads on log10 scale')  +",
           "  xlab('Group') + ",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold'))",
           "```\n",

           "## PRIMER DIMER filter  \n",

           "```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0\n",

           "ggplot(data = config, aes(x = Experiment_Type, y = PD_percentage, order = Experiment_Type, fill = Experiment_Type)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads marked as PRIMER DIMERS')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'none') +",
           "  ylim(0, 100)",
           "``` \n",

           "## Cut rates  \n",

           "```{r plot cut percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Experiment_Type, y = cut_percentage, order = Experiment_Type, fill = Experiment_Type)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have cut site')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14, face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0,100)",
           "``` \n",

           "## Frameshift  \n",

           "```{r plot frame shift percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Experiment_Type, y = frameshift_percentage, order = Experiment_Type, fill = Experiment_Type)) +",
           "  geom_boxplot() +",
           "  xlab('Group') + ",
           "  ylab('Percentage of reads (not marked as PRIMER DIMERS) that have frameshift')  +",
           "  theme(axis.text = element_text(size=12),",
           "        axis.title = element_text(size=14,face = 'bold'),",
           "        legend.position = 'None') +",
           "  ylim(0, 100)",
           "``` \n",

           "## Heterogeneity of reads  \n",

           "```{r plot read heterogeneity, echo=FALSE, fig.height=15, fig.width=14, message=F, warning=FALSE}",
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$ex_type <- rep(config$Experiment_Type[match(howManyTimes$seqnames, unique(config$ID))],",
                                          "times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$ex_type, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$ex_type), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$ex_type))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$ex_type))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$ex_type))]",
           "uniqueReadsByID$read_number <- unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes, by = 1))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$ex_type), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$ex_type))]",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "ggplot(data = uniqueReadsByID, aes(x = ex_type, y = read_share_percentage_normal, fill = read_share_percentage_normal)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('Group') +",
           "  coord_flip()",
           "``` \n"))
}


#' Make contents for .Rmd file aggregated on the guideRNA
#'
#' @param results_folder path to results as string
#' @return c() collection of strings to put into rmd file
#'
make_guide_rmd <- function(results_folder) {
  return(c(write_head("Report breakdown by guideRNA"),
           "```{r load data, message=F, warning=FALSE, include=FALSE}",
           paste0("results_folder = '", results_folder, "'"),
           "library(amplican)",
           "library(ggplot2)",
           "alignments <- read.csv(paste0(results_folder, '/alignments.csv'), sep = '\\t')",
           "config <- read.csv(paste0(results_folder, '/config_summary.csv'), sep = '\\t')",
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

           "```{r plot total reads, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "ggplot(data = config, aes(x = Target_Primer, y = log10(Reads + 1), order = Target_Primer, fill = Target_Primer)) +",
           "  geom_boxplot() +",
           "  ylab('Number of reads on log10 scale')  +",
           "  xlab('guideRNA') + ",
           "  theme(legend.position = 'none',",
           "        axis.text = element_text(size = 12),",
           "        axis.title = element_text(size = 14, face = 'bold')) +",
           "  coord_flip()",
           "```\n",

           "## PRIMER DIMER filter  \n",

           "```{r plot PRIMER DIMER percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$PD_percentage <- config$PRIMER_DIMER * 100/config$Reads",
           "config$PD_percentage[is.nan(config$PD_percentage)] <- 0\n",

           "ggplot(data = config, aes(x = Target_Primer, y = PD_percentage, order = Target_Primer, fill = Target_Primer)) +",
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

           "```{r plot cut percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$Reads_noPD <- config$Reads - config$PRIMER_DIMER",
           "config$cut_percentage <- config$Cut * 100/config$Reads_noPD",
           "config$cut_percentage[is.nan(config$cut_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Target_Primer, y = cut_percentage, order = Target_Primer, fill = Target_Primer)) +",
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

           "```{r plot frame shift percentage, echo=FALSE, fig.height=14, fig.width=14, message=F, warning=FALSE}",
           "config$frameshift_percentage <- config$Frameshift * 100/config$Reads_noPD",
           "config$frameshift_percentage[is.nan(config$frameshift_percentage)] <- 0  \n",

           "ggplot(data = config, aes(x = Target_Primer, y = frameshift_percentage, order = Target_Primer, fill = Target_Primer)) +",
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

           "```{r plot read heterogeneity, echo=FALSE, fig.height=15, fig.width=14, message=F, warning=FALSE}",
           "alignments$ID_read_id <- paste0(alignments$seqnames, '_', alignments$read_id)",
           "uniqueReadsByID <- alignments[!duplicated(alignments$ID_read_id), c('seqnames', 'read_id', 'count')]",
           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$seqnames, uniqueReadsByID$count, decreasing = T),]\n",

           "howManyTimes <- aggregate(read_id ~ seqnames, data = uniqueReadsByID, length)",
           "howManyTimes <- howManyTimes[match(howManyTimes$seqnames, unique(uniqueReadsByID$seqnames)),] \n",

           "uniqueReadsByID$guideRNA <- rep(config$Target_Primer[match(howManyTimes$seqnames, unique(config$ID))], \n",
           "                                times = howManyTimes$read_id)\n",

           "uniqueReadsByID <- uniqueReadsByID[order(uniqueReadsByID$guideRNA, uniqueReadsByID$count, decreasing = T),]\n",

           "cumsum_list <- tapply(uniqueReadsByID$count, as.vector(uniqueReadsByID$guideRNA), FUN = cumsum)",
           "cumsum_list <- cumsum_list[match(names(cumsum_list), unique(uniqueReadsByID$guideRNA))]",
           "uniqueReadsByID$cumsum <- do.call(c, cumsum_list)",

           "seq2 <- Vectorize(seq.default, vectorize.args = c('from', 'to'))",
           "howManyTimes <- table(as.vector(uniqueReadsByID$guideRNA))",
           "howManyTimes <- howManyTimes[match(names(howManyTimes), unique(uniqueReadsByID$guideRNA))]",
           "uniqueReadsByID$read_number <- unlist(seq2(from = rep(1, dim(howManyTimes)[1]), to = howManyTimes, by = 1))\n",

           "ids_with_reads <- tapply(uniqueReadsByID$cumsum, as.vector(uniqueReadsByID$guideRNA), FUN = max)",
           "ids_with_reads <- ids_with_reads[match(names(ids_with_reads), unique(uniqueReadsByID$guideRNA))]\n",

           "toDivide <- rep(ids_with_reads, times = howManyTimes)",
           "uniqueReadsByID$read_share_percentage <- uniqueReadsByID$cumsum * 100 / toDivide",
           "uniqueReadsByID$read_share_percentage_normal <- uniqueReadsByID$count * 100 / toDivide  \n",

           "ggplot(data = uniqueReadsByID, aes(x = guideRNA, y = read_share_percentage_normal, fill = read_share_percentage_normal)) +",
           "  geom_bar(position = 'stack', stat = 'identity') +",
           "  theme(legend.position = 'top',",
           "        axis.text = element_text(size=12),",
           "        axis.title=element_text(size=14,face = 'bold'),",
           "        legend.direction = 'horizontal', ",
           "        legend.title = element_blank()) +",
           "  ylab('Unique reads percentage of contribution') +",
           "  xlab('guideRNA') +",
           "  coord_flip()",
           "``` \n"))
}
