#' Align reads to amplicons.
#'
#' amplicanAlign takes a configuration files, fastq reads and output
#' directory to prepare alignments and summary. It uses global Needleman-Wunsch
#' alghoritm with parameters optimized for CRISPR experiment. After alignments,
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
  total_processors = 1,
  average_quality = 30,
  min_quality = 20,
  scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(
    match = 5, mismatch = -4, baseOnly = TRUE, type = "DNA"),
  gap_opening = 50,
  gap_extension = 0,
  fastqfiles = 0) {

  message("Checking configuration file...")
  cfgT <- readr::read_csv(config, col_types = "ccccclccclc", na = "")
  cfgT[is.na(cfgT)] <- ""
  cfgT <- data.frame(cfgT)
  colnames(cfgT) <- c("ID", "Barcode", "Forward_Reads_File",
                      "Reverse_Reads_File", "Group", "Control", "guideRNA",
                      "Forward_Primer", "Reverse_Primer", "Direction",
                      "Amplicon")
  cfgT$Forward_Reads_File <-
    ifelse(cfgT$Forward_Reads_File == "",
           "",
           file.path(fastq_folder, cfgT$Forward_Reads_File))
  cfgT$Reverse_Reads_File <-
    ifelse(cfgT$Reverse_Reads_File == "",
           "",
           file.path(fastq_folder, cfgT$Reverse_Reads_File))

  if (sum(cfgT$Reverse_Reads_File == "") > 0) {
    message(paste0("Reverse_Reads_File has empty rows. ",
                   "Changing fastqfiles parameter to 1, ",
                   "operating only on forward reads."))
    fastqfiles <- 1
  }
  if (sum(cfgT$Forward_Reads_File == "") > 0) {
    message(paste0("Forward_Reads_File has empty rows. ",
                   "Changing fastqfiles parameter to 2, ",
                   "operating only on reverse reads."))
    fastqfiles <- 2
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

  if (total_processors > 1) {
    p <- BiocParallel::MulticoreParam(workers = total_processors)
    BiocParallel::register(p)
    configSplit <- split(cfgT, f = cfgT$Barcode)
    finalAES <- BiocParallel::bplapply(configSplit, FUN = makeAlignment,
                                       average_quality,
                                       min_quality,
                                       scoring_matrix,
                                       gap_opening,
                                       gap_extension,
                                       fastqfiles, BPPARAM=p)
  } else {
    finalAES <- vector("list", length(uBarcode))
    for (j in seq_along(uBarcode)) {
      finalAES[[j]] <- makeAlignment(cfgT[cfgT$Barcode == uBarcode[j], ],
                                     average_quality,
                                     min_quality,
                                     scoring_matrix,
                                     gap_opening,
                                     gap_extension,
                                     fastqfiles)
    }
  }
  BiocGenerics::Reduce(c, finalAES)
}
