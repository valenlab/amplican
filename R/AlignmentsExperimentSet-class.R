methods::setClassUnion("listOrNULL", members = c("list", "NULL"))
methods::setClassUnion("data.frameOrNULL", members = c("data.frame", "NULL"))

#' An S4 class to represent alignments from multiple experiments
#'
#' Class \code{AlignmentsExperimentSet} holds data from multiple alignments for
#' many experiments. Allows to examine alignments in great detail.
#'
#'
#' @param x,object (AlignmentsExperimentSet)
#' @param value Represents assignment values for setter methods.
#' @param i,j,name,drop (numeric, missing, character, missing)
#' AlignmentsExperimentSet object can be subsetted using
#' names of the experiments eg. \code{x$name} or \code{x[i]}, resulting in
#' AlignmentsExperimentSet object that has only one experiment. During this
#' subsetting, values of unassignedData and barcodeData are dropped.
#' @param use_parallel (boolean) When using \code{extractEvents} you can
#' use multicore back-end through \code{\link[BiocParallel]{register}} as this
#' is very slow function (despite vectorization).
#' @param ... pass any number of AlignmentsExperimentSet objects, make sure
#' experiment IDs can be unique after merging
#' @inheritParams writeAlignments
#' @inheritParams lookupAlignment
#' @return depending on the function used
#'
#' @slot fwdReads,rveReads (list) Named list where each element is of class
#' \code{\link{PairwiseAlignmentsSingleSubject}}. Names correspond
#' to the experiment ID. Contains alignments of reads against amplicons.
#' @slot fwdReadsType,rveReadsType (list) Named list where each element is of
#' logical vector, so far TRUE corresponds to HDR events. Names correspond
#' to the experiment ID. Contains type of read - HDR/NHEJ.
#' @slot readCounts (list) Named list where each element is numeric vector that
#' describes how many reads are compressed into unique representation
#' before alignment in \code{fwdReads} and/or \code{rveReads}.
#' @slot unassignedData (data.frame) Contains reads that failed to be assigned
#' to any of the experiments. Alignment of forward against reverse reads may
#' give hint whether these reads are compromised in any way.
#' @slot experimentData (data.frame) Expands on configuration file and provides
#' information about cut rates, frameshifts, PRIMER DIMER detection etc. Each
#' row corresponds to experiment ID.
#' @slot barcodeData (data.frame) Information that is gathered on the barcode
#' level is gathered in this data.frame, mainly quality filtering statistics.
#'
#' @include helpers_general.R
#'
#' @name AlignmentsExperimentSet-class
#' @rdname AlignmentsExperimentSet-class
#' @exportClass AlignmentsExperimentSet
#' @examples
#' exampleAlignments <- Biostrings::pairwiseAlignment(
#'   Biostrings::DNAStringSet(c("ACTGACTG", "CGACGACG")), "ACGTACGTACGT")
#' new("AlignmentsExperimentSet",
#'      fwdReads = list(ID_1 = exampleAlignments, ID_2 = exampleAlignments),
#'      rveReads = list(ID_1 = exampleAlignments, ID_2 = exampleAlignments),
#'      fwdReadsType = list(ID_1 = c(FALSE, FALSE), ID_2 = c(FALSE, FALSE)),
#'      rveReadsType = list(ID_1 = c(FALSE, FALSE), ID_2 = c(FALSE, FALSE)),
#'      readCounts = list(ID_1 = c(2, 20), ID_2 = c(30, 100)),
#'      unassignedData = NULL,
#'      experimentData = data.frame(ID = c("ID_1", "ID_2"),
#'                                  Barcode = c("B1", "B1"),
#'                                  whatever = c(50, 100)),
#'      barcodeData = data.frame(Barcode = "B1", statistic1 = 100))
setClass("AlignmentsExperimentSet",
         representation(fwdReads = "listOrNULL",
                        rveReads = "listOrNULL",
                        fwdReadsType = "listOrNULL",
                        rveReadsType = "listOrNULL",
                        readCounts = "listOrNULL",
                        unassignedData = "data.frameOrNULL",
                        experimentData = "data.frameOrNULL",
                        barcodeData = "data.frameOrNULL"),
         prototype(fwdReads = NULL,
                   rveReads = NULL,
                   fwdReadsType = NULL,
                   rveReadsType = NULL,
                   readCounts = NULL,
                   unassignedData = NULL,
                   experimentData = NULL,
                   barcodeData = NULL))


listClasses <- function(lst) {
  unique(sapply(lst, class))
}


setValidity("AlignmentsExperimentSet", function(object) {
  err <- character()

  if (!is.null(readCounts(object)) &&
      any(!listClasses(readCounts(object)) %in%
          c("integer", "numeric", "NULL"))) {
    err <- c(err, paste("'readCounts' has to have all elements",
                        "of class 'integer' or 'numeric'"))
  }

  if (!is.null(readCounts(object)) &&
      any(sapply(readCounts(object), function(x) any(x == 0)))) {
    err <- c(err, paste("'readCounts' has to have all elements",
                        "higher than 0."))
  }

  if (!is.null(fwdReads(object)) &&
      any(!listClasses(fwdReads(object)) %in%
          c("PairwiseAlignmentsSingleSubject", "NULL"))) {
    err <- c(err, paste("'fwdReads' has to have all elements",
                        "of class PairwiseAlignmentsSingleSubject"))
  }

  if (!is.null(rveReads(object)) &&
      any(!listClasses(rveReads(object)) %in%
          c("PairwiseAlignmentsSingleSubject", "NULL"))) {
    err <- c(err, paste("'rveReads' has to have all elements",
                        "of class PairwiseAlignmentsSingleSubject"))
  }

  if (!is.null(fwdReads(object)) && !is.null(rveReads(object)) &&
      length(fwdReads(object)) != length(rveReads(object))) {
    err <- c(err, paste("'fwdReads' and 'rveReads' should have",
                        "the same number of experiments (same length)."))
  }

  if (!is.null(readCounts(object)) && !is.null(experimentData(object)) &&
      length(readCounts(object)) != dim(experimentData(object))[1]) {
    err <- c(err, paste("Number of rows in 'experimentData' should be same",
                        "as length of 'readCounts'."))
  }

  if (!is.null(rveReads(object)) &&
      length(readCounts(object)) != length(rveReads(object))) {
    err <- c(err, paste("Length of 'readCounts' should be same",
                        "as length of 'rveReads'."))
  }

  if (!is.null(fwdReads(object)) &&
      length(readCounts(object)) != length(fwdReads(object))) {
    err <- c(err, paste("Length of 'readCounts' should be same",
                        "as length of 'fwdReads'."))
  }

  if (!is.null(experimentData(object)) && is.null(experimentData(object)$ID)) {
    err <- c(err, paste("There should be ID column in 'experimentData'."))
  }

  if (!is.null(experimentData(object)) &&
      is.null(experimentData(object)$Barcode)) {
    err <- c(err, paste("There should be Barcode column in 'experimentData'."))
  }

  if (!is.null(unassignedData(object)) &&
      is.null(unassignedData(object)$Barcode)) {
    err <- c(err, paste("There should be Barcode column in 'unassignedData'."))
  }

  if (!is.null(barcodeData(object)) &&
      is.null(barcodeData(object)$Barcode)) {
    err <- c(err, paste("There should be Barcode column in 'barcodeData'."))
  }

  if (!all(c(barcodeData(object)$Barcode, unassignedData(object)$Barcode) %in%
           c(experimentData(object)$Barcode, NULL))) {
    err <- c(err, paste("There are unknown Barcodes without any experiments.",
                        "Check your 'experimentData' in relation to",
                        "'unassigedData' and 'barcodeData'."))
  }

  if (length(unique(experimentData(object)$ID)) !=
      length(experimentData(object)$ID)) {
    err <- c(err, paste("All experiment names should be unique."))
  }

  if ((!is.null(fwdReads(object)) &&
       !all(names(fwdReads(object)) == experimentData(object)$ID)) ||
      (!is.null(rveReads(object)) &&
       !all(names(rveReads(object)) == experimentData(object)$ID))) {
    err <- c(err, paste("Names of 'fwdReads' should be same ",
                        "as names of 'rveReads', ''readCounts' ",
                        "and 'experimentData' ID column."))
  }

  if (length(err) == 0) TRUE else err
})


#' @name AlignmentsExperimentSet
#' @rdname AlignmentsExperimentSet-class
#' @export
AlignmentsExperimentSet <- function(...) {
  methods::new("AlignmentsExperimentSet", ...)
}


#' @aliases length,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("length", "AlignmentsExperimentSet", function(x) {
  if (!is.null(experimentData(x))) dim(experimentData(x))[1] else 0
})


#' Alignments for forward reads.
#'
#' Get alignments for forward reads.
#' @name fwdReads
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' fwdReads(aln)
setGeneric("fwdReads", function(x) standardGeneric("fwdReads"))
#' @aliases fwdReads,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReads", "AlignmentsExperimentSet", function(x) x@fwdReads)


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name fwdReads<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' fwdReads(aln) <- fwdReads(aln) # replace with the same values
setGeneric("fwdReads<-", function(x, value) standardGeneric("fwdReads<-"))
#' @aliases fwdReads<-,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReads<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, fwdReads = value)
})


#' Alignments for reverse reads.
#'
#' Get alignments for reverse reads.
#' @name rveReads
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' rveReads(aln)
setGeneric("rveReads", function(x) standardGeneric("rveReads"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases rveReads,AlignmentsExperimentSet-method
setMethod("rveReads", "AlignmentsExperimentSet", function(x) x@rveReads)


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name rveReads<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' rveReads(aln) <- rveReads(aln) # replace with the same values
setGeneric("rveReads<-", function(x, value) standardGeneric("rveReads<-"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases rveReads<-,AlignmentsExperimentSet-method
setMethod("rveReads<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, rveReads = value)
})


#' Type of forward reads.
#'
#' Get type of forward reads.
#' @name fwdReadsType
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' fwdReadsType(aln)
setGeneric("fwdReadsType", function(x) standardGeneric("fwdReadsType"))
#' @aliases fwdReadsType,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReadsType", "AlignmentsExperimentSet", function(x) x@fwdReadsType)


#' Read type for forward reads.
#'
#' Set read type for forward reads.
#' @name fwdReadsType<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' fwdReadsType(aln) <- fwdReadsType(aln) # replace with the same values
setGeneric("fwdReadsType<-", function(x, value) {
  standardGeneric("fwdReadsType<-")
})
#' @aliases fwdReadsType<-,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReadsType<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, fwdReadsType = value)
})


#' Type of reverse reads.
#'
#' Get type of reverse reads.
#' @name rveReadsType
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' rveReadsType(aln)
setGeneric("rveReadsType", function(x) standardGeneric("rveReadsType"))
#' @aliases rveReadsType,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("rveReadsType", "AlignmentsExperimentSet", function(x) x@rveReadsType)


#' Read type for reverse reads.
#'
#' Set read type for reverse reads.
#' @name rveReadsType<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' rveReadsType(aln) <- rveReadsType(aln) # replace with the same values
setGeneric("rveReadsType<-", function(x, value) {
  standardGeneric("rveReadsType<-")
})
#' @aliases rveReadsType<-,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("rveReadsType<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, rveReadsType = value)
})


#' Unassigned read information.
#'
#' Get unassigned reads and their characteristics.
#' @name unassignedData
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (data.frameOrNULL)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' unassignedData(aln)
setGeneric("unassignedData", function(x) standardGeneric("unassignedData"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases unassignedData,AlignmentsExperimentSet-method
setMethod("unassignedData", "AlignmentsExperimentSet", function(x) {
  x@unassignedData
})


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name unassignedData<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' unassignedData(aln) <- unassignedData(aln) #replace with the same values
setGeneric("unassignedData<-", function(x, value) {
  standardGeneric("unassignedData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases unassignedData<-,AlignmentsExperimentSet-method
setMethod("unassignedData<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, unassignedData = value)
})


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name readCounts
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (listOrNULL)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' readCounts(aln)
setGeneric("readCounts", function(x) standardGeneric("readCounts"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases readCounts,AlignmentsExperimentSet-method
setMethod("readCounts", "AlignmentsExperimentSet", function(x) {
  x@readCounts
})


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name readCounts<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' @return (AlignmentsExperimentSet)
#' \code{\link{PairwiseAlignmentsSingleSubject}} class.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' readCounts(aln) <- readCounts(aln) # replace with the same values
setGeneric("readCounts<-", function(x, value) {
  standardGeneric("readCounts<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases readCounts<-,AlignmentsExperimentSet-method
setMethod("readCounts<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, readCounts = value)
})


#' Experiment data.
#'
#' Get experiment data.frame with information on the experiment level.
#' @name experimentData
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (data.frameOrNULL)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' experimentData(aln)
setGeneric("experimentData", function(x) standardGeneric("experimentData"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases experimentData,AlignmentsExperimentSet-method
setMethod("experimentData", "AlignmentsExperimentSet", function(x) {
  x@experimentData
})


#' Experiment data.
#'
#' Set experiment data.frame with information on the experiment level.
#' @name experimentData<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (data.frame)
#' @return (AlignmentsExperimentSet)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' experimentData(aln) <- experimentData(aln) # replace with the same values
setGeneric("experimentData<-", function(x, value) {
  standardGeneric("experimentData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases experimentData<-,AlignmentsExperimentSet-method
setMethod("experimentData<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, experimentData = value)
})


#' Barcode data.
#'
#' Get barcode data.frame with information on the barcode level.
#' @name barcodeData
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (data.tableOrNULL)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' barcodeData(aln)
setGeneric("barcodeData", function(x) standardGeneric("barcodeData"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases barcodeData,AlignmentsExperimentSet-method
setMethod("barcodeData", "AlignmentsExperimentSet", function(x) {
  x@barcodeData
})


#' Barcode data.
#'
#' Set barcode data.frame with information on the barcode level.
#' @name barcodeData<-
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param value (data.frame)
#' @return (AlignmentsExperimentSet)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' barcodeData(aln) <- barcodeData(aln) #replace with the same values as before
setGeneric("barcodeData<-", function(x, value) {
  standardGeneric("barcodeData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases barcodeData<-,AlignmentsExperimentSet-method
setMethod("barcodeData<-", "AlignmentsExperimentSet", function(x, value) {
  initialize(x, barcodeData = value)
})


#' Get count of unassigned reads.
#'
#' @name unassignedCount
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (numeric)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' unassignedCount(aln)
setGeneric("unassignedCount", function(x) standardGeneric("unassignedCount"))
#' @rdname AlignmentsExperimentSet-class
#' @keywords internal
#' @aliases unassignedCount,AlignmentsExperimentSet-method
setMethod("unassignedCount", "AlignmentsExperimentSet", function(x) {
  if (!is.null(unassignedData(x))) dim(unassignedData(x))[1] else 0
})


#' Get count of assigned reads.
#'
#' @name assignedCount
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @return (numeric)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' writeAlignments(aln, file.path(tempdir(), "aln.txt"))
#'
setGeneric("assignedCount", function(x) standardGeneric("assignedCount"))
#' @rdname AlignmentsExperimentSet-class
#' @keywords internal
#' @aliases assignedCount,AlignmentsExperimentSet-method
setMethod("assignedCount", "AlignmentsExperimentSet", function(x) {
  if (length(readCounts(x)) == 0) return(0)
  sum(sapply(readCounts(x), sum))
})


#' @aliases names,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("names", "AlignmentsExperimentSet", function(x) {
  as.character(experimentData(x)$ID)
})


#' @rdname AlignmentsExperimentSet-class
#' @aliases c,AlignmentsExperimentSet-method
setMethod("c", "AlignmentsExperimentSet", function(x, ...) {
  args <- if (missing(x)) list(...) else (list(x, ...))
  aln <- methods::new("AlignmentsExperimentSet",
               fwdReads = do.call(c, lapply(args, fwdReads)),
               rveReads = do.call(c, lapply(args, rveReads)),
               fwdReadsType = do.call(c, lapply(args, fwdReadsType)),
               rveReadsType = do.call(c, lapply(args, rveReadsType)),
               readCounts = do.call(c, lapply(args, readCounts)),
               unassignedData = Reduce(rbind, lapply(args, unassignedData)),
               experimentData = Reduce(rbind, lapply(args, experimentData)),
               barcodeData = Reduce(rbind, lapply(args, barcodeData)))
})

init <- function(x, i) {
  # drop barcodeData and unassignedData as they are probable not needed when
  # user subsets and can lead to problems if user later uses "c" as we
  # can't know which barcodeData and unassignedData belongs to which experiment
  # which would lead to combining these data.frames with repetitions
  initialize(x,
             fwdReads = fwdReads(x)[i],
             rveReads = rveReads(x)[i],
             fwdReadsType = fwdReadsType(x)[i],
             rveReadsType = rveReadsType(x)[i],
             readCounts = readCounts(x)[i],
             experimentData = experimentData(x)[i, ],
             barcodeData = NULL,
             unassignedData = NULL)
}
#' @rdname AlignmentsExperimentSet-class
setMethod("[", c("AlignmentsExperimentSet", "numeric", "missing", "missing"),
          function(x, i, j, ..., drop=TRUE) {
            if (any(i > length(x))) stop("Subsetting out of bonds.")
            init(x, i)
          })


#' @rdname AlignmentsExperimentSet-class
#' @method as.list AlignmentsExperimentSet
#' @export
as.list.AlignmentsExperimentSet <- function(x, ...) {
  lapply(seq_along(x), function(i) x[i])
}


#' @aliases $,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("$", "AlignmentsExperimentSet", function(x, name) {
  i = which(names(x) == name)
  init(x, i)
})


#' Write alignments to file.
#'
#' Saves alignments into txt or fasta file.
#' @name writeAlignments
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param file (connection or string) Destination file. When empty, defaults to
#' standard output.
#' @param aln_format ("txt" or "fasta") Specifies format of the file.
#' @return (invisible)
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' writeAlignments(aln, file.path(tempdir(), "aln.txt"))
#'
setGeneric("writeAlignments", function(x, file = "", aln_format = "txt") {
  standardGeneric("writeAlignments")
})
#' @section View alignments:
#' Write out all alignments in "fasta" or "txt" format.:  \cr
#' \code{writeAlignments(x, file = "", aln_format = "txt")}
#' @rdname AlignmentsExperimentSet-class
#' @aliases writeAlignments,AlignmentsExperimentSet-method
#'
setMethod("writeAlignments", "AlignmentsExperimentSet", function(
  x, file = "", aln_format = "txt") {

  if (!((is.character(file) && length(file) == 1L && !is.na(file)) |
        methods::is(file, "connection"))) {
    stop("'file' must be a single string or a connection object")
  }

  if (file == "") {
    file <- stdout()
  } else if (substring(file, 1L, 1L) == "|") {
    file <- pipe(substring(file, 2L), "w")
    on.exit(close(file))
  } else {
    file <- file(file, "w")
    on.exit(close(file))
  }

  if (aln_format == "txt") {
    for (ID in x@experimentData$ID) {
      if (is.null(readCounts(x)[[ID]])) next()
      writeLines(as.vector(rbind(
        paste("ID:", ID,
              "read_id:",
              format(seq_len(length(readCounts(x)[[ID]]))),
              "Count:", format(readCounts(x)[[ID]])),
        if (length(fwdReads(x)[[ID]]) > 0) rbind(
          as.character(Biostrings::pattern(fwdReads(x)[[ID]])),
          as.character(Biostrings::subject(fwdReads(x)[[ID]]))),
        if (length(fwdReads(x)[[ID]]) == length(rveReads(x)[[ID]])) "",
        if (length(rveReads(x)[[ID]]) > 0) rbind(
          as.character(Biostrings::pattern(rveReads(x)[[ID]])),
          as.character(Biostrings::subject(rveReads(x)[[ID]]))), "")),
        file)
    }
  }

  if (aln_format == "fasta") {
    for (ID in x@experimentData$ID) {
      counts = readCounts(x)[[ID]]
      if (is.null(counts)) next()
      writeLines(as.vector(rbind(
        if (length(fwdReads(x)[[ID]]) > 0)
          rbind(paste(">Forward read ID:", ID,
                      "read_id:",
                      format(seq_len(length(fwdReads(x)[[ID]]))),
                      "Count:", format(counts),
                      "Type:", format(fwdReadsType(x)[[ID]])),
                as.character(Biostrings::pattern(fwdReads(x)[[ID]])),
                paste(">Forward amplicon ID:", ID,
                      "read_id:",
                      format(seq_len(length(fwdReads(x)[[ID]]))),
                      "Count:", format(counts)),
                as.character(Biostrings::subject(fwdReads(x)[[ID]]))),
        if (length(rveReads(x)[[ID]]) > 0)
          rbind(paste(">Reverse read ID:", ID,
                      "read_id:",
                      format(seq_len(length(rveReads(x)[[ID]]))),
                      "Count:", format(counts),
                      "Type:", format(rveReadsType(x)[[ID]])),
                as.character(Biostrings::pattern(rveReads(x)[[ID]])),
                paste(">Reverse amplicon ID:", ID,
                      "read_id:",
                      format(seq_len(length(rveReads(x)[[ID]]))),
                      "Count:", format(counts)),
                as.character(Biostrings::subject(rveReads(x)[[ID]]))))
      ), file)
    }
  }

  invisible(file)
})


#' Show alignment in human readable format.
#'
#' Prints alignments in blast-like style for human examination.
#' @name lookupAlignment
#' @keywords internal
#' @param x (AlignmentsExperimentSet)
#' @param ID (string) Experiment Identifier
#' @param read_id (numeric) Read Identifier. Reads are sorted by frequency.
#' Defaults to 1, most abundant read.
#' @return (print to view)
#' @export
#' @examples
#' # load example object
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' # look at most frequent reads aligned from experiment ID_1
#' lookupAlignment(aln, "ID_1")
#'
setGeneric("lookupAlignment", function(x, ID, read_id = 1){
  standardGeneric("lookupAlignment")
})
#' @section View alignments:
#' Write out human readable alignments for given experiment and read_id.:  \cr
#' \code{lookupAlignment(x, ID, read_id = 1)}
#' @rdname AlignmentsExperimentSet-class
#' @aliases lookupAlignment,AlignmentsExperimentSet-method
setMethod("lookupAlignment", "AlignmentsExperimentSet", function(
  x, ID, read_id = 1) {

  if (length(readCounts(x)[[ID]]) == 0) return(NULL)
  if (length(fwdReads(x)[[ID]]) >= read_id) {
    fwd <- utils::capture.output(Biostrings::writePairwiseAlignments(
      fwdReads(x)[[ID]][as.numeric(read_id)]))
    fwd[7] <- "# Aligned sequences:"
    fwd[8] <- "# Forward read: P1"
    fwd[9] <- "# Amplicon sequence: S1"
    fwd <- c(fwd[1], paste0("# Forward read for ID ", ID, " read_id ", read_id),
             fwd[2:length(fwd)])
  } else {fwd <- NULL}

  if (length(rveReads(x)[[ID]]) >= read_id) {
    rve <- utils::capture.output(Biostrings::writePairwiseAlignments(
      rveReads(x)[[ID]][as.numeric(read_id)]))
    rve[7] <- "# Aligned sequences:"
    rve[8] <- "# Reverse read: P1"
    rve[9] <- "# Amplicon sequence: S1"
    rve <- c(rve[1], paste0("# Reverse read for ID ", ID, " read_id ", read_id),
             rve[2:length(rve)])
  } else {rve <- NULL}

  if (is.null(rve) & is.null(fwd)) stop("Subsetting out of bounds.")
  writeLines(c(fwd, rve))
})


getEventInfoObj <- function(object) {
  if (length(object) != 1) stop("Only length 1 objects.")
  ID <- names(object)[1]
  cfg <- experimentData(object)
  cfg <- cfg[cfg$ID == ID]
  fwdPrPos <- if (is.null(cfg$fwdPrPos)) 1 else cfg$fwdPrPos
  tempGR <- c(getEventInfo(fwdReads(object)[[ID]], ID, fwdPrPos, "+"),
              getEventInfo(rveReads(object)[[ID]], ID, fwdPrPos, "-"))
  tempGR$counts <- readCounts(object)[[ID]][as.integer(tempGR$read_id)]
  if (length(tempGR) > 0) tempGR$readType <- FALSE

  plus_strand <- as.vector(strand(tempGR) == "+")
  fwd_ids <- as.integer(tempGR$read_id)[plus_strand]
  fRT <- fwdReadsType(object)[[ID]]
  if (!is.null(fRT) & sum(plus_strand) > 0 & length(fwd_ids) > 0) {
    tempGR$readType[plus_strand] <- fRT[fwd_ids]
  }

  minus_strand <- as.vector(strand(tempGR) == "-")
  rve_ids <- as.integer(tempGR$read_id)[minus_strand]
  rRT <- rveReadsType(object)[[ID]]
  if (!is.null(rRT) & sum(minus_strand) > 0 & length(rve_ids) > 0) {
    tempGR$readType[minus_strand] <- rRT[rve_ids]
  }
  tempGR
}
#' Extract AlignmentsExperimentSet events into data.frame.
#'
#' Extracts events (insertions, deletions, mismatches) from alignments into
#' data.frame. Can use multiple cores as process is quite slow. All events are
#' relative towards forward strand. "-" in strand column indicates which events
#' were from reverse reads.
#' @name extractEvents
#' @keywords internal
#' @param object (AlignmentsExperimentSet)
#' @param use_parallel (boolean) Set to TRUE, if you have registered
#' multicore back-end with \code{\link[BiocParallel]{register}}.
#' @return (data.frame) Compatible with \code{\link{GRanges}}
#' style.
#' @export
#' @examples
#' file_path <- system.file("extdata", "results", "alignments",
#'                          "AlignmentsExperimentSet.rds", package = "amplican")
#' aln <- readRDS(file_path)
#' extractEvents(aln)
#'
setGeneric("extractEvents", function(object, use_parallel = FALSE){
  standardGeneric("extractEvents")
})
#' @section Coercion based on events:
#' Coerce to \code{data.frame} compatible with
#' \code{\link{GRanges}} .: \cr
#' as.data.frame(x)
#' @aliases extractEvents,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
#' @examples
#' # Coercion
#' extractEvents(AlignmentsExperimentSet())
#' GenomicRanges::GRanges(extractEvents(AlignmentsExperimentSet()))
#'
setMethod("extractEvents", "AlignmentsExperimentSet", function(
  object, use_parallel = FALSE) {

  if (length(readCounts(object)) == 0) {
    return(BiocGenerics::as.data.frame(GenomicRanges::GRanges()))
  }

  p <- if (!use_parallel) {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::bpparam()
  }
  finalGR <- BiocParallel::bplapply(object, FUN = getEventInfoObj, BPPARAM = p)
  finalGR <- unlist(GenomicRanges::GRangesList(finalGR), use.names = FALSE)
  flipRanges(GenomicRanges::as.data.frame(finalGR, row.names = NULL),
             experimentData(object))
})


setMethod("show", "AlignmentsExperimentSet", function(object){
  cat("\nAlignmentsExperimentSet object with", length(object), "experiments.\n")
  if (length(object) > 0) {
    cat("Showing only first experiment with ID: ",
        names(object)[1], "\n\n", sep = "")
    cat("Slot fwdReads:\n")
    print(fwdReads(object)[[1]])
    cat("\nSlot rveReads:\n")
    print(rveReads(object)[[1]])
    cat("\nSlot fwdReadsType:\n")
    print(table(fwdReadsType(object)[[1]]))
    cat("\nSlot rveReadsType:\n")
    print(table(rveReadsType(object)[[1]]))
    cat("\nSlot readCounts:\n")
    utils::str(readCounts(object)[[1]])
    cat("\nSlot experimentData:\n")
    utils::str(experimentData(object)[1, ])
  }
})

