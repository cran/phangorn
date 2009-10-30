\name{bootstrap.pml}
%\Rdversion{1.1}
\alias{bootstrap.pml}
\alias{bootstrap.phyDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap }
\description{
\code{bootstrap.pml} performs (non-parametric) bootstrap analysis and \code{bootstrap.phyDat} produces a list of bootstrapped datasets. }
\usage{
bootstrap.pml(x, bs = 3, trees = TRUE, ...)
bootstrap.phyDat(x, FUN, bs = 100, ...) 
}
\arguments{
  \item{x}{
an object of class pml or phyDat.
}
  \item{bs}{
number of bootstrap samples.
}
  \item{trees}{
return trees only (default) or whole pml objects.  
}
  \item{\dots}{
further parameters used by optim.pml.
}
  \item{FUN}{ 
the function to estimate the trees.
}
}
\value{
returns an object of class multi.phylo or a list where each 
element is an object of class pml.
}
\author{
Klaus Schliep \email{klaus.schliep@gmail.com}
}
\seealso{
\code{\link{optim.pml}}, \code{\link{pml}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

