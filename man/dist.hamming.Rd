\name{dist.hamming}
\alias{dist.hamming}
\alias{dist.logDet}
\alias{dist.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Distances from Sequences}
\description{
\code{dist.hamming} and \code{dist.logDet} compute pairwise distances 
for an object of class \code{phyDat}. \code{dist.ml} fits distances
for amino acid models.
}
\usage{
dist.hamming(x, ratio = TRUE)
dist.logDet(x)
dist.ml(x, model="JC69",...)
}
\arguments{
  \item{x}{an object of class \code{dist.logDet}}
  \item{ratio}{compute uncorrected ('p') distance or character difference.}
  \item{model}{One of "JC69", "WAG", "JTT", "LG" or "Dayhoff"}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
  an object of class \code{dist}
}
\references{ 
Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994)
Recovering evolutionary trees under a more realistic model of sequence
evolution. \emph{Molecular Biology and Evolution}, \bold{11}, 605--602.
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{For more distance methods for nucleotide data see \code{\link[ape]{dist.dna}} }
\examples{
data(Laurasiatherian)
dm1 <- dist.hamming(Laurasiatherian)
tree1 <- NJ(dm1)
dm2 <- dist.logDet(Laurasiatherian)
tree2 <- NJ(dm2)
treedist(tree1,tree2)
}

\keyword{ cluster }
