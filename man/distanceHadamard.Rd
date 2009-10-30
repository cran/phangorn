\name{distanceHadamard}

\alias{distanceHadamard}
\title{Distance Hadamard}
\description{
Distance Hadamard produces spectra of splits from a distance matrix.
}
\usage{
distanceHadamard(dm)
}
\arguments{
  \item{dm}{A distance matrix.}
}

\value{
\code{distanceHadamard} returns a matrix. The first column 
contains the distance spectra, the second one the edge-spectra.
}
\references{Hendy, M. D. and Penny, D. (1993). Spectral Analysis of Phylogenetic
Data. \emph{Journal of Classification}, \bold{10}, 5-24.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}, Tim White}
\seealso{\code{\link[phangorn]{hadamard}}}

\examples{
data(yeast)
dm = dist.hamming(yeast)
dm = as.matrix(dm)
fit = distanceHadamard(dm)
}
\keyword{ cluster }% at least one, from doc/KEYWORDS


