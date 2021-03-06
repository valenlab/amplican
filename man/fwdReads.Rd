% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AlignmentsExperimentSet-class.R
\name{fwdReads}
\alias{fwdReads}
\title{Alignments for forward reads.}
\usage{
fwdReads(x)
}
\arguments{
\item{x}{(AlignmentsExperimentSet)}
}
\value{
(listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
}
\description{
Get alignments for forward reads.
}
\examples{
file_path <- system.file("extdata", "results", "alignments",
                         "AlignmentsExperimentSet.rds", package = "amplican")
aln <- readRDS(file_path)
fwdReads(aln)
}
\keyword{internal}
