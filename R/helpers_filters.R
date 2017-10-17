range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' Find Off-targets and Fragmented alignments from reads.
#'
#' Will try to detect off-targets and low quality alignments (outliers). It
#' tries k-means clustering on normalized number of events per read and read
#' alignment score. If there are 3 clusters (decided based on silhouette
#' criterion) cluster with high event count and low alignment score will be
#' marked for filtering. When there is less than 1000
#' scores in \code{aln} it will filter nothing.
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end, score.
#' @return (logical vector) where TRUE indicates events that are
#' potential off-targets or low quality alignments.
#' @export
#' @family filters
#' @seealso \code{\link{findPD}} \code{\link{findEOP}}
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' aln <- aln[seqnames == "ID_1"] # for first experiment
#' findLQR(aln)
#'
findLQR <- function(aln) {
  data.table::setDT(aln)
  if (dim(aln)[1] < 1000) return(logical(dim(aln)[1]))
  events <- NULL

  aln_n <- aln[, list(events = .N/max(end), score = max(score)),
               by = c("read_id", "strand", "seqnames")]
  aln_n <- aln_n[, list(events = events/length(unique(strand)),
                        score = score/length(unique(strand))),
                 by = c("read_id", "seqnames")]
  x <- cbind(range01(aln_n$score), range01(aln_n$events))

  k2 <- stats::kmeans(x, 2, iter.max = 1000, nstart = 1000)
  k2s <- clusterCrit::intCriteria(x, k2$cluster, "silhouette")$silhouette
  k3 <- stats::kmeans(x, 3, iter.max = 1000, nstart = 1000)
  k3s <- clusterCrit::intCriteria(x, k3$cluster, "silhouette")$silhouette
  if (k2s >= k3s) return(logical(dim(aln)[1])) else {
    # find top left center and filter it
    centers <- apply(k3$centers, 1,
                     function(x) sqrt((x[1] - 1) ^ 2 + x[2] ^ 2))
    bs <- aln_n[k3$cluster == which.max(centers)]
    return(aln$seqnames %in% bs$seqnames & aln$read_id %in% bs$read_id)
  }
}


#' Find Events Overlapping Primers.
#'
#' Very often alignments return deletions that are not real deletions, but
#' rather artifact of incomplete reads eg.: \cr
#' \preformatted{
#' ACTGAAAAA------- <- this "deletion" should be filtered
#' ACTG----ACTGACTG
#' }
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end.
#' @param cfgT (data.frame) Needs columns Forward_Primer, ReversePrimer and
#' Amplicon.
#' @return (logical vector) where TRUE indicates events that are overlapping
#' primers
#' @export
#' @family filters
#' @seealso \code{\link{findPD}} \code{\link{findLQR}}
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' findEOP(aln, cfgT)
#'
findEOP <- function(aln, cfgT) {
  mapID <- match(aln$seqnames, cfgT$ID)
  if (any(aln$start < 0 | aln$end < 0)) { # if events are relative
    for (i in seq_along(cfgT$ID)) {
      amplicon <- get_amplicon(cfgT, cfgT$ID[i])
      zero_point <- upperGroups(amplicon)
      if (length(zero_point) == 0) next()
      cfgT$fwdPrPosEnd[i] <- cfgT$fwdPrPosEnd[i] - 1 *
        GenomicRanges::start(zero_point)[1]
      cfgT$rvePrPos[i] <- cfgT$rvePrPos[i] - 1 *
        GenomicRanges::start(zero_point)[1]
    }
  }
  (aln$start < cfgT$fwdPrPosEnd[mapID]) |
    (aln$end > cfgT$rvePrPos[mapID] & aln$type != "insertion") |
    (aln$start > cfgT$rvePrPos[mapID] & aln$type == "insertion")
}


#' Find PRIMER DIMER reads.
#'
#' Use to filter reads that are most likely PRIMER DIMERS.
#' @param aln (data.frame) Should contain events from alignments in GRanges
#' style with columns eg. seqnames, width, start, end.
#' @param cfgT (data.frame) Needs columns Forward_Primer, ReversePrimer and
#' Amplicon.
#' @param PRIMER_DIMER (numeric) Value specifying buffer for PRIMER DIMER
#' detection. For a given read it will be recognized as PRIMER DIMER when
#' alignment will introduce gap of size bigger than: \cr
#' \code{length of amplicon - (lengths of PRIMERS + PRIMER_DIMER value)}
#' @return (logical) Where TRUE indicates event classified as PRIMER DIMER
#' @export
#' @family filters
#' @seealso \code{\link{findEOP}} \code{\link{findLQR}}
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "raw_events.csv", package = "amplican")
#' aln <- data.table::fread(file_path)
#' cfgT <- data.table::fread(
#'   system.file("extdata", "results", "config_summary.csv",
#'               package = "amplican"))
#' findPD(aln, cfgT)
#'
findPD <- function(aln, cfgT, PRIMER_DIMER = 30) {

  PD_cutoff <- nchar(cfgT$Amplicon) -
    (nchar(cfgT$Forward_Primer) + nchar(cfgT$Reverse_Primer) + PRIMER_DIMER)
  PD_cutoff <- PD_cutoff[match(aln$seqnames, cfgT$ID)]

  aln$width > PD_cutoff
}


#' Filters out sequences which have bad base quality readings.
#'
#' @keywords internal
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param min (numeric) This is the minimum quality that we accept for
#' every nucleotide. For example, if we have a sequence with nucleotides which
#' have quality 50-50-50-50-10, and we set the minimum to 30, the whole sequence
#' will be a bad sequence. The minimum is set to 0 by default.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#'
goodBaseQuality <- function(reads, min = 20) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(matrixStats::rowMins(methods::as(quality(reads), "matrix"),
                              na.rm = TRUE) >= min)
}


#' This filters out sequences which have bad average quality readings.
#'
#' @keywords internal
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @param avg (numeric) This is what the average score of the quality of
#' sequence should be. For example, if we have a sequence with nucleotides which
#' have quality 70-70-70, the average would be 70. If set the average to 70 or
#' less the sequence will pass. If we set the average to 71 the sequence will
#' not pass. The average is set to 0 by default.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#'
goodAvgQuality <- function(reads, avg = 30) {
  if (is.logical(reads)) {
    return(reads)
  }
  return(Matrix::rowMeans(methods::as(quality(reads), "matrix"),
                          na.rm = TRUE) >= avg)
}


#' This filters out sequences which have nonstandard nucleotides.
#'
#' @keywords internal
#' @param reads (ShortRead object) Loaded reads from fastq.
#' @return (boolean) Logical vector with the valid rows as TRUE.
#'
alphabetQuality <- function(reads) {
  if (is.logical(reads)) {
    return(reads)
  }

  nucq <- ShortRead::nFilter()
  return(as.logical(nucq(reads)))
}
