#' An S4 class to represent alignments from multiple experiments
#'
#' Class \code{AlignmentsExperimentSet} holds data from multiple alignments for
#' many experiments. Allows to examine alignments in great detail.
#'
#'
#' @param x,object (AlignmentsExperimentSet)
#' @param value Represents assignment values for setter methods.
#' @param ... pass any number of AlignmentsExperimentSet objects, make sure
#' experiment IDs can be unique after merging
#' @inheritParams writeAlignments
#' @inheritParams lookupAlignment
#'
#' @slot fwdReads,rveReads (list) Named list where each element is of class
#' \code{\link[Biostrings]{PairwiseAlignmentsSingleSubject}}. Names correspond
#' to the experiment ID. Contains alignments of reads against amplicons.
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
#' a <- new("AlignmentsExperimentSet")
#' summary(a)
#'
setClass("AlignmentsExperimentSet",
         representation(fwdReads = "list",
                        rveReads = "list",
                        readCounts = "list",
                        unassignedData = "data.frame",
                        experimentData = "data.frame",
                        barcodeData = "data.frame"),
         prototype(fwdReads = list(),
                   rveReads = list(),
                   readCounts = list(),
                   unassignedData = data.frame(),
                   experimentData = data.frame(),
                   barcodeData = data.frame()))

listClasses <- function(lst) {
  unique(sapply(lst, class))
}

setValidity("AlignmentsExperimentSet", function(object) {
  err <- character()

  if (length(object@readCounts) != 0 &
      any(!listClasses(object@readCounts) %in% c("integer", "NULL"))) {
    err <- c(err, paste("'readCounts' has to have all elements ",
                        "of class 'integer'"))
  }

  if (length(object@fwdReads) != 0 &
      any(!listClasses(object@fwdReads) %in%
          c("PairwiseAlignmentsSingleSubject", "NULL"))) {
    err <- c(err, paste("'fwdReads' has to have all elements ",
                        "of class PairwiseAlignmentsSingleSubject"))
  }

  if (length(object@rveReads) != 0 &
      any(!listClasses(object@rveReads) %in%
          c("PairwiseAlignmentsSingleSubject", "NULL"))) {
    err <- c(err, paste("'rveReads' has to have all elements ",
                        "of class PairwiseAlignmentsSingleSubject"))
  }

  if (length(object@fwdReads) != 0 & length(object@rveReads) != 0 &
      length(object@fwdReads) != length(object@fwdReads)) {
    err <- c(err, paste("'fwdReads and 'rveReads' should have ",
                        "the same number of experiments (same length)."))
  }

  if (length(object@readCounts) != dim(object@experimentData)[1]) {
    err <- c(err, paste("Number of rows in 'experimentData' should be same ",
                        "as length of 'readCounts'."))
  }

  if (length(object@rveReads) != 0 &
      length(object@readCounts) != length(object@rveReads)) {
    err <- c(err, paste("Length of 'readCounts' should be same ",
                        "as length of 'rveReads'."))
  }

  if (length(object@fwdReads) != 0 &
      length(object@readCounts) != length(object@fwdReads)) {
    err <- c(err, paste("Length of 'readCounts' should be same ",
                        "as length of 'fwdReads'."))
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
  dim(x@experimentData)[1]
})


#' @aliases names,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("names", "AlignmentsExperimentSet", function(x) {
  x@experimentData$ID
})


#' @aliases names<-,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setReplaceMethod("names", "AlignmentsExperimentSet", function(x, value) {
  if (length(x) != length(value)) {
    stop("New names must be of the same length as previous.")
  }
  x@experimentData$ID <- value
  #names(fwdReads(x)) <- value todo
  #names(rveReads(x)) <- value
})


#' Alignments for forward reads.
#'
#' Get alignments for forward reads.
#' @name fwdReads
#' @param x (AlignmentsExperimentSet)
#' @export
setGeneric("fwdReads", function(x) standardGeneric("fwdReads"))
#' @aliases fwdReads,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReads", "AlignmentsExperimentSet", function(x) x@fwdReads)


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name fwdReads<-
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' \code{\link[Biostrings]{PairwiseAlignmentsSingleSubject}} class.
#' @export
setGeneric("fwdReads<-", function(x, value) standardGeneric("fwdReads<-"))
#' @aliases fwdReads<-,AlignmentsExperimentSet-method
#' @rdname AlignmentsExperimentSet-class
setMethod("fwdReads<-", "AlignmentsExperimentSet", function(x, value) {
  x@fwdReads <- value # change to initialize for validity check and speed
  x
})


#' Alignments for reverse reads.
#'
#' Get alignments for reverse reads.
#' @name rveReads
#' @param x (AlignmentsExperimentSet)
#' @export
setGeneric("rveReads", function(x) standardGeneric("rveReads"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases rveReads,AlignmentsExperimentSet-method
setMethod("rveReads", "AlignmentsExperimentSet", function(x) x@rveReads)


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name rveReads<-
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' \code{\link[Biostrings]{PairwiseAlignmentsSingleSubject}} class.
#' @export
setGeneric("rveReads<-", function(x, value) standardGeneric("rveReads<-"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases rveReads<-,AlignmentsExperimentSet-method
setMethod("rveReads<-", "AlignmentsExperimentSet", function(x, value) {
  x@rveReads <- value
  x
})


#' Unassigned read information.
#'
#' Get unassigned reads and their characteristics.
#' @name unassignedData
#' @param x (AlignmentsExperimentSet)
#' @export
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
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' \code{\link[Biostrings]{PairwiseAlignmentsSingleSubject}} class.
#' @export
setGeneric("unassignedData<-", function(x, value) {
  standardGeneric("unassignedData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases unassignedData<-,AlignmentsExperimentSet-method
setMethod("unassignedData<-", "AlignmentsExperimentSet", function(x, value) {
  x@unassignedData <- value
  x
})


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name readCounts
#' @param x (AlignmentsExperimentSet)
#' @export
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
#' @param x (AlignmentsExperimentSet)
#' @param value (list) Named (experiment IDs) list with elements of
#' \code{\link[Biostrings]{PairwiseAlignmentsSingleSubject}} class.
#' @export
setGeneric("readCounts<-", function(x, value) {
  standardGeneric("readCounts<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases readCounts<-,AlignmentsExperimentSet-method
setMethod("readCounts<-", "AlignmentsExperimentSet", function(x, value) {
  x@readCounts <- value
  x
})


#' Experiment data.
#'
#' Get experiment data.frame with information on the experiment level.
#' @name experimentData
#' @param x (AlignmentsExperimentSet)
#' @export
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
#' @param x (AlignmentsExperimentSet)
#' @param value (data.frame)
#' @export
setGeneric("experimentData<-", function(x, value) {
  standardGeneric("experimentData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases experimentData<-,AlignmentsExperimentSet-method
setMethod("experimentData<-", "AlignmentsExperimentSet", function(x, value) {
  x@experimentData <- value
  x
})


#' Barcode data.
#'
#' Get barcode data.frame with information on the barcode level.
#' @name barcodeData
#' @param x (AlignmentsExperimentSet)
#' @export
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
#' @param x (AlignmentsExperimentSet)
#' @param value (data.frame)
#' @export
setGeneric("barcodeData<-", function(x, value) {
  standardGeneric("barcodeData<-")
})
#' @rdname AlignmentsExperimentSet-class
#' @aliases barcodeData<-,AlignmentsExperimentSet-method
setMethod("barcodeData<-", "AlignmentsExperimentSet", function(x, value) {
  x@barcodeData <- value
  x
})


#' Get count of unassigned reads.
#'
#' @name unassignedCount
#' @param x (AlignmentsExperimentSet)
#' @export
setGeneric("unassignedCount", function(x) standardGeneric("unassignedCount"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases unassignedCount,AlignmentsExperimentSet-method
setMethod("unassignedCount", "AlignmentsExperimentSet", function(x) {
  dim(unassignedData(x))[1]
})


#' Get count of assigned reads.
#'
#' @name assignedCount
#' @param x (AlignmentsExperimentSet)
#' @export
setGeneric("assignedCount", function(x) standardGeneric("assignedCount"))
#' @rdname AlignmentsExperimentSet-class
#' @aliases assignedCount,AlignmentsExperimentSet-method
setMethod("assignedCount", "AlignmentsExperimentSet", function(x) {
  if (length(readCounts(x)) == 0) return(0)
  sum(sapply(readCounts(x), sum))
})


#' @rdname AlignmentsExperimentSet-class
#' @aliases c,AlignmentsExperimentSet-method
setMethod("c", "AlignmentsExperimentSet", function(x, ...) {
  args <- if (missing(x)) list(...) else (list(x, ...))
  methods::new("AlignmentsExperimentSet",
               fwdReads = Reduce(c, lapply(args, fwdReads)),
               rveReads = Reduce(c, lapply(args, rveReads)),
               readCounts = Reduce(c, lapply(args, readCounts)),
               unassignedData = Reduce(rbind, lapply(args, unassignedData)),
               experimentData = Reduce(rbind, lapply(args, experimentData)),
               barcodeData = Reduce(rbind, lapply(args, barcodeData)))
})


#' @rdname AlignmentsExperimentSet-class
#' @aliases summary,AlignmentsExperimentSet-method
setMethod("summary", "AlignmentsExperimentSet", function(object){
  writeLines(
    c(paste("AlignmentsExperimentSet object with", length(object),
            "experiments."),
      paste("Unassgined reads %:",
            (unassignedCount(object)/assignedCount(object)) * 100),
      paste("Number of unique guideRNA's:",
            length(unique(object@experimentData$guideRNA)))))
})


#' Alignments for forward reads.
#'
#' Set alignments for forward reads.
#' @name writeAlignments
#' @param x (AlignmentsExperimentSet)
#' @param file (connection or string) Destination file. When empty, defaults to
#' standard output.
#' @param aln_format ("txt" or "fasta") Specifies format of the file.
#' @export
setGeneric("writeAlignments", function(x, file = "", aln_format = "txt") {
  standardGeneric("writeAlignments")
})
#' @section View alignments:
#' Write out all alignments in "fasta" or "txt" format.:  \cr
#' \code{writeAlignments(x, file = "", aln_format = "txt")}
#' @rdname AlignmentsExperimentSet-class
#' @aliases writeAlignments,AlignmentsExperimentSet-method
setMethod("writeAlignments", "AlignmentsExperimentSet",
          function(x, file = "", aln_format = "txt") {

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
                  if (length(fwdReads(x)) > 0) rbind(
                    as.character(Biostrings::pattern(fwdReads(x)[[ID]])),
                    as.character(Biostrings::subject(fwdReads(x)[[ID]]))),
                  if (length(fwdReads(x)) == length(rveReads(x))) "",
                  if (length(rveReads(x)) > 0) rbind(
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
                  if (length(fwdReads(x)) > 0)
                    rbind(paste(">Forward read ID:", ID,
                                "read_id:",
                                format(seq_len(length(fwdReads(x)[[ID]]))),
                                "Count:", format(counts)),
                          as.character(Biostrings::pattern(fwdReads(x)[[ID]])),
                          paste(">Forward amplicon ID:", ID,
                                "read_id:",
                                format(seq_len(length(fwdReads(x)[[ID]]))),
                                "Count:", format(counts)),
                          as.character(Biostrings::subject(fwdReads(x)[[ID]]))),
                  if (length(rveReads(x)) > 0)
                    rbind(paste(">Reverse read ID:", ID,
                                "read_id:",
                                format(seq_len(length(rveReads(x)[[ID]]))),
                                "Count:", format(counts)),
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
#' @param x (AlignmentsExperimentSet)
#' @param ID (string) Experiment Identifier
#' @param read_id (numeric) Read Identifier. Reads are sorted by frequency.
#' Defaults to 1, most abundant read.
#' @export
setGeneric("lookupAlignment", function(x, ID, read_id = 1){
  standardGeneric("lookupAlignment")
})
#' @section View alignments:
#' Write out human readable alignments for given experiment and read_id.:  \cr
#' \code{lookupAlignment(x, ID, read_id = 1)}
#' @rdname AlignmentsExperimentSet-class
#' @aliases lookupAlignment,AlignmentsExperimentSet-method
setMethod("lookupAlignment", "AlignmentsExperimentSet",
          function(x, ID, read_id = 1) {

            if (length(readCounts(x)[[ID]]) == 0) return(NULL)
            if (length(fwdReads(x)) > 0) {
              fwd <- capture.output(Biostrings::writePairwiseAlignments(
                fwdReads(x)[[ID]][as.numeric(read_id)]))
              fwd[7] <- "# Aligned sequences:"
              fwd[8] <- "# Forward read: P1"
              fwd[9] <- "# Amplicon sequence: S1"
              fwd <- c(fwd[1], paste0("# Forward read for ID ", ID, " read_id ", read_id),
                       fwd[2:length(fwd)])
            } else {fwd <- NULL}

            if (length(rveReads(x)) > 0) {
              rve <- capture.output(Biostrings::writePairwiseAlignments(
                fwdReads(x)[[ID]][as.numeric(read_id)]))
              rve[7] <- "# Aligned sequences:"
              rve[8] <- "# Reverse read: P1"
              rve[9] <- "# Amplicon sequence: S1"
              rve <- c(rve[1], paste0("# Reverse read for ID ", ID, " read_id ", read_id),
                       rve[2:length(rve)])
            } else {rve <- NULL}

            writeLines(c(fwd, rve))
          })


#' @section Coercion:
#' Coerce to \code{data.frame} compatible with
#' \code{\link[GenomicRanges]{GRanges}} .: \cr
#' as(x, "data.frame")
#' @name as
#' @rdname AlignmentsExperimentSet-class
#' @examples
#' # Coercion
#' as(AlignmentsExperimentSet(), "data.frame")
#' GenomicRanges::GRanges(as(AlignmentsExperimentSet(), "data.frame"))
#'
setAs("AlignmentsExperimentSet", "data.frame", function(from) {
  if (length(readCounts(from)) == 0) {
    return(BiocGenerics::as.data.frame(GenomicRanges::GRanges()))
  }
  cfg <- experimentData(from)
  finalGR <- unlist(GenomicRanges::GRangesList(
    lapply(experimentData(from)$ID, function(ID) {
      tempGR <- c(if (length(fwdReads(from)[[ID]]) > 0) {
        getEventInfo(fwdReads(from)[[ID]], ID, cfg$fwdPrPos[cfg$ID == ID], "+")
      } else {
        GenomicRanges::GRanges()
      }, if (length(rveReads(from)[[ID]]) > 0) {
        getEventInfo(rveReads(from)[[ID]], ID,
                     cfg$rvePrPosEnd[cfg$ID == ID], "-")
      } else {
        GenomicRanges::GRanges()
      })

      tempGR$counts <- readCounts(from)[[ID]][as.integer(tempGR$read_id)]
      tempGR
    })), use.names = FALSE)

  flipRanges(GenomicRanges::as.data.frame(finalGR, row.names = NULL),
             experimentData(from))
})



