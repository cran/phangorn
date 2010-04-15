\name{as.splits}
\alias{as.splits}
\alias{as.splits.phylo}
\alias{as.splits.multiPhylo}
\alias{compatible}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splits representation of graphs and trees.
}
\description{
\code{as.splits} produces a list of splits or bipartitions. 
}
\usage{
as.splits(x, ...)
\method{as.splits}{phylo}(x, ...)
\method{as.splits}{multiPhylo}(x, ...)
compatible(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class phylo or multiPhylo.}
  \item{\dots}{Further arguments passed to or from other methods.}
  \item{obj}{an object of class splits.}  
}
\value{
\code{as.splits} returns an object of class splits, which is mainly a list of splits.
\code{compatible} return a lower triangular matrix where an 1 indicates that two splits are incompatible.   
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\note{
The internal representation is likely to change. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prop.part}}, \code{\link{lento}}, \code{\link{distanceHadamard}}
}
%% as.phylo.splits
\examples{
as.splits(rtree(5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

