\name{as.splits}
\alias{as.splits}
\alias{as.splits.phylo}
\alias{as.splits.multiPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splits representation of graphs and trees
}
\description{
\code{as.splits} 
}
\usage{
as.splits(x, ...)
\method{as.splits}{phylo}(x, ...)
\method{as.splits}{multiPhylo}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class phylo or multiPhylo.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
An object of class splits, which is mainly a list of splits.  
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\note{
The internal representation is likely to change. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prop.part}}, \code{\link{lento}}, \code{\link{distanceHadamard}}
}
\examples{
as.splits(rtree(5))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

