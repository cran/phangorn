\name{upgma}
\alias{upgma}
\title{ UPGMA }
\description{
  UPGMA clustering. Just a wrapper function around \code{\link[stats]{hclust}}.
}
\usage{
upgma(D, method = "average", ...)
}
\arguments{
  \item{D}{A distance matrix.}
  \item{method}{The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
A phylogenetic tree of class \code{phylo}. 
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{ \code{\link[stats]{hclust}}, \code{\link[ape]{as.phylo}} }
\examples{
data(Laurasiatherian)
dm = dist.logDet(Laurasiatherian)
tree = upgma(dm)
plot(tree)

}
\keyword{cluster}
