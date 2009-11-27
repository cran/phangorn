\name{designTree}
\alias{designTree}
\alias{designSplits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute a design matrix }
\description{
\code{designTree} and \code{designSplits} compute design matrices for the 
estimation of edge length of (phylogenetic) trees using linear models. 
\code{designTree} also computes a contrast matrix if the method is "rooted".
}
\usage{
designTree(tree, method = "unrooted", ...)
designSplits(x, splits = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{an object of class \code{phylo} }
  \item{method}{design matrix for an "unrooted" or "rooted" ultrametric tree}
  \item{x}{number of taxa. }
  \item{splits}{one of "all", "star"}
  \item{\dots}{further arguments, passed to other methods.  }
}
\value{
 a matrix. 
}
\author{ Klaus Schliep \email{klaus.schliep@gmail.com} }

\seealso{ \code{\link[ape]{fastme}}, \code{\link[phangorn]{distanceHadamard}} }
\examples{
example(NJ)
dm <-  as.matrix(dm)
y <- dm[lower.tri(dm)]
X <- designTree(tree)
summary(lm(y~X-1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
