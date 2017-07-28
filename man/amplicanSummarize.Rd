% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amplicanSummarize.R
\name{amplicanSummarize}
\alias{amplicanSummarize}
\title{Summarize how many reads has frameshift, frameshift overlapping cut site and
how many reads with deletions overlaping cut site.}
\usage{
amplicanSummarize(aln, cfgT, overlaps = "overlaps")
}
\arguments{
\item{aln}{(data.frame) Contains events from the alignments.}

\item{cfgT}{(data.frame) Config file with the experiments details.}

\item{overlaps}{(character) Specifies which column of aln contains
output from \code{\link{amplicanOverlap}}. Defaults to "overlaps".}
}
\value{
(data.frame) As cfgT, but with extra columns.
}
\description{
Adds columns to cfgT:
\itemize{
\item{ReadsCut}{ Count of reads with deletions overlapping expected
cut site.}
\item{Reads_Frameshifted}{ Count of reads with frameshift overlapping
expected cut site.}
\item{Reads_Frameshifted_Overlapping}{ Count of reads with frameshift
overlapping expected cut site.}
}
}
\seealso{
Other analysis steps: \code{\link{amplicanAlign}},
  \code{\link{amplicanFilter}},
  \code{\link{amplicanNormalize}},
  \code{\link{amplicanOverlap}},
  \code{\link{amplicanPipeline}},
  \code{\link{amplicanReport}},
  \code{\link{map_to_relative}}
}