% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeManipulation.R
\name{allTrees}
\alias{allTrees}
\title{Compute all trees topologies.}
\usage{
allTrees(n, rooted = FALSE, tip.label = NULL)
}
\arguments{
\item{n}{Number of tips (<=10).}

\item{rooted}{Rooted or unrooted trees (default: rooted).}

\item{tip.label}{Tip labels.}
}
\value{
an object of class \code{multiPhylo}.
}
\description{
\code{allTrees} computes all bifurcating tree topologies for rooted or
unrooted trees with up to 10 tips. The number of trees grows fast.
}
\examples{

trees <- allTrees(5)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(3,5))
for(i in 1:15)plot(trees[[i]])
par(old.par)

}
\seealso{
\code{\link[ape]{rtree}}, \code{\link{nni}},
\code{\link[ape]{howmanytrees}}, \code{\link{dfactorial}}
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\keyword{cluster}
