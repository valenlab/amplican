% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AlignmentsExperimentSet-class.R
\name{writeAlignments}
\alias{writeAlignments}
\title{Write alignments to file.}
\usage{
writeAlignments(x, file = "", aln_format = "txt")
}
\arguments{
\item{x}{(AlignmentsExperimentSet)}

\item{file}{(connection or string) Destination file. When empty, defaults to
standard output.}

\item{aln_format}{("txt" or "fasta") Specifies format of the file.}
}
\value{
(invisible)
}
\description{
Saves alignments into txt or fasta file.
}
\examples{
file_path <- system.file("extdata", "results", "alignments",
                         "AlignmentsExperimentSet.rds", package = "amplican")
aln <- readRDS(file_path)
writeAlignments(aln, file.path(tempdir(), "aln.txt"))

}
\keyword{internal}
