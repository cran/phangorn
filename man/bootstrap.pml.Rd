\name{bootstrap.pml}
%\Rdversion{1.1}
\alias{bootstrap.pml}
\alias{bootstrap.phyDat}
\alias{plotBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap }
\description{
\code{bootstrap.pml} performs (non-parametric) bootstrap analysis and \code{bootstrap.phyDat} produces a list of bootstrapped datasets. Under Linux the bootstrap is performed in parallel, if the multicore package is loaded. \code{plotBS} plots a phylogenetic tree with the with the bootstrap values assigned to the (internal) edges.}
\usage{
bootstrap.pml(x, bs = 100, trees = TRUE, ...)
bootstrap.phyDat(x, FUN, bs = 100, ...)
plotBS(tree, BStrees, type="unrooted", bs.col="black", bs.adj=c(0.5, 0.5), ...)
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
further parameters used by optim.pml or plot.phylo.
}
  \item{FUN}{ 
the function to estimate the trees.
}
  \item{tree}{ 
The tree on which edges the bootstrap values are plotted.
}
  \item{BStrees}{ 
a list of trees (object of class "multiPhylo"). 
}
  \item{type}{ 
the type of tree.
}
  \item{bs.col}{ 
color of bootstrap support labels.
}
  \item{bs.adj}{ 
one or two numeric values specifying the horizontal and vertical justification of the bootstrap labels.
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
\code{\link{optim.pml}}, \code{\link{pml}}, \code{\link{plot.phylo}}
}
\examples{
\dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit=pml(tree,Laurasiatherian)
fit = optim.pml(fit,TRUE)

set.seed(1)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
plotBS(fit$tree,bs)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

