% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AlignmentsExperimentSet-class.R
\name{rveReads}
\alias{rveReads}
\title{Alignments for reverse reads.}
\usage{
rveReads(x)
}
\arguments{
\item{x}{(AlignmentsExperimentSet)}
}
\value{
(listOrNULL) list with objects of PairwiseAlignmentsSingleSubject
}
\description{
Get alignments for reverse reads.
}
\examples{
file_path <- system.file("extdata", "results", "alignments",
                         "AlignmentsExperimentSet.rds", package = "amplican")
aln <- readRDS(file_path)
rveReads(aln)
}
\keyword{internal}
