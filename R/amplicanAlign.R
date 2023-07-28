#' Align reads to amplicons.
#'
#' amplicanAlign takes a configuration files, fastq reads and output
#' directory to prepare alignments and summary. It uses global Needleman-Wunsch
#' algorithm with parameters optimized for CRISPR experiment. After alignments,
#' object of \code{\link{AlignmentsExperimentSet}} is returned that allows for
#' coercion into GRanges (plus is for forward and minus for reverse reads).
#' It is also possible to output alignments in other, additional formats.
#' @inheritParams amplicanPipeline
#' @return (AlignmentsExperimentSet) Check \code{\link{AlignmentsExperimentSet}}
#' class for details. You can use \code{\link{lookupAlignment}} to examine
#' alignments visually.
#' @include helpers_alignment.R helpers_warnings.R helpers_directory.R
#' AlignmentsExperimentSet-class.R
#' @export
#' @family analysis steps
#' @examples
#' # path to example config file
#' config <- system.file("extdata", "config.csv", package = "amplican")
#' # path to example fastq files
#' fastq_folder <- system.file("extdata", package = "amplican")
#' aln <- amplicanAlign(config, fastq_folder)
#' aln
#'
amplicanAlign <- function(
  config,
  fastq_folder,
  use_parallel = FALSE,
  average_quality = 30,
  min_quality = 20,
  filter_n = FALSE,
  batch_size = 1e6,
  scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
    match = 5, mismatch = -4, baseOnly = FALSE, type = "DNA"),
  gap_opening = 25,
  gap_extension = 0,
  fastqfiles = 0.5,
  primer_mismatch = 0,
  donor_mismatch = 3,
  donor_strict = FALSE) {

  message("Checking configuration file...")
  cfgT <- data.frame(data.table::fread(config))
  cfgT[is.na(cfgT)] <- ""
  colnames(cfgT) <- c("ID", "Barcode", "Forward_Reads_File",
                      "Reverse_Reads_File", "Group", "Control", "guideRNA",
                      "Forward_Primer", "Reverse_Primer", "Direction",
                      "Amplicon", "Donor",
                      if (dim(cfgT)[2] > 12) colnames(cfgT)[13:dim(cfgT)[2]])
  ids <- cfgT$ID
  cfgT$Control <- as.logical(cfgT$Control)
  cfgT$Direction <- as.logical(cfgT$Direction)
  cfgT$Forward_Reads_File <-
    ifelse(cfgT$Forward_Reads_File == "",
           "",
           file.path(fastq_folder, cfgT$Forward_Reads_File))
  cfgT$Reverse_Reads_File <-
    ifelse(cfgT$Reverse_Reads_File == "",
           "",
           file.path(fastq_folder, cfgT$Reverse_Reads_File))

  if (sum(cfgT$Reverse_Reads_File == "") > 0 & fastqfiles != 1) {
    stop("Reverse_Reads_File has empty rows. ",
         "Change fastqfiles parameter to 1, ",
         "to operate only on forward reads.")
  }
  if (sum(cfgT$Forward_Reads_File == "" & fastqfiles != 2) > 0) {
    stop("Forward_Reads_File has empty rows. ",
         "Change fastqfiles parameter to 2, ",
         "to operate only on reverse reads.")
  }
  checkConfigFile(cfgT, fastq_folder)

  cfgT$Reverse_PrimerRC <- revComp(cfgT$Reverse_Primer)
  cfgT <- checkPrimers(cfgT, fastqfiles)

  cfgT$guideRNA[cfgT$Direction] <- revComp(cfgT$guideRNA[cfgT$Direction])
  cfgT$RguideRNA <- revComp(cfgT$guideRNA)
  cfgT$Found_Guide <- checkTarget(cfgT)
  cfgT$ampl_len <- nchar(cfgT$Amplicon)

  uBarcode <- unique(cfgT$Barcode)
  cfgT$Reads <- 0

  p <- if (!use_parallel) {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::bpparam()
  }
  message("Making alignments...")
  config_order <- cfgT$ID
  configSplit <- split(cfgT, f = cfgT$Barcode)
  finalAES <- BiocParallel::bplapply(configSplit, FUN = makeAlignment,
                                     average_quality,
                                     min_quality,
                                     filter_n,
                                     batch_size,
                                     scoring_matrix,
                                     gap_opening,
                                     gap_extension,
                                     fastqfiles,
                                     primer_mismatch,
                                     donor_mismatch,
                                     donor_strict, BPPARAM = p)
  finalAES <- Reduce(c, finalAES)

  # sort like at the entry point
  cfgT <- experimentData(finalAES)
  id_order <- match(ids, cfgT$ID)
  initialize(finalAES,
             fwdReads = fwdReads(finalAES)[id_order],
             rveReads = rveReads(finalAES)[id_order],
             fwdReadsType = fwdReadsType(finalAES)[id_order],
             rveReadsType = rveReadsType(finalAES)[id_order],
             readCounts = readCounts(finalAES)[id_order],
             experimentData = cfgT[id_order, ])
}
