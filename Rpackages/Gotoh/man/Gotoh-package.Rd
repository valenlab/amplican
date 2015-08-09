\name{Gotoh-package}
\alias{Gotoh-package}
\alias{Gotoh}
\docType{package}
\title{
C++ wrapper for R which does alignments using Gotoh algorithm.
}
\description{
Implement a Gotoh approach to sequence alignment. The Gotoh algorithm was first
published by Osamu Gotoh in 1982; for more information please refer to:
<doi:10.1016/0022-2836(82)90398-9>
We are using 'Rcpp', which decrease greatly the CPU time with respect other
alignment implementations. This was done in aneffort to speed up the
calculations for the project AmpliCan. The improvement we got compared with
'DECIPHER' is in between 300x and 500x many fold for sequences between 100 to
200 nucleotides, and between 100x and 200x for sequences of between 200 to 300
nucleotides. And up to 1.5x for sequences of size 10.000 nucleotides and
beyond. This package contains the minimum files require for it to work. Please
download the repository at:
https://rafanozal@bitbucket.org/rafanozal/amplican.git
to get the a more verbose documentation and script benchmarks.
}
\details{
\tabular{ll}{
Package: \tab Gotoh\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2015-07-18\cr
License: \tab GPL (>= 2)\cr
}
In order to use the package you need to install the Rcpp library.
There are only one function: "gRCPP". For gRCPP you just need to give two
strings representing a sequence of nucleotides each; any other option have a
default value. Please read the help manual for a verbose description.
}
\author{
Rafael Nozal Canyadas <rca025@uib.no>
Eivind Valen <eivind.valen@gmail.com>

Maintainer: Rafael Nozal Canyadas <rca025@uib.no>
}
\references{
http://goo.gl/pIDucs

(If the link is dead, just search "Gotoh 1982".)
}
\keyword{ package }
\examples{
gRCPP("aaaaa", "aaddttaaaa")
gRCPP("aaaaa", "aaddttaaaa", "NUC44", 30, 5)
gRCPP("aaaaa", "aaddttaaaa", "NUC44", 30, 5, TRUE, TRUE)
}
