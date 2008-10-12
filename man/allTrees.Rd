\name{allTrees}
\alias{allTrees}
\alias{nni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute all trees topologies and nearest neighbor interchange.}
\description{
\code{allTrees} computes all tree topologies for rooted or unrooted trees with up to 10 tips.
\code{nni} returns a list of all trees which are one nearest neighbor interchange away. 
Both methods assume bifurcating trees.
}
\usage{
allTrees(n, rooted = FALSE, tip.label = NULL)
nni(tree,  edge.length = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of tips (<=10).}
  \item{rooted}{Rooted or unrooted trees (default: rooted).  }
  \item{tip.label}{Tip labels.}
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}.}
  \item{edge.length}{Keep edge lengths?}
}
\value{
   an object of class phylo. 
}
\author{ Klaus Schliep \email{K.P.Schliep@massey.ac.nz} }

\seealso{ \code{\link[ape]{fastme}}, \code{\link[phangorn]{distanceHadamard}} }
\examples{
trees <- allTrees(5)
par(mfrow = c(3,5))
for(i in 1:15)plot(trees[[i]])

neighborhood <- nni(trees[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
