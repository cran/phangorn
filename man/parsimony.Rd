\name{parsimony}
\alias{parsimony}
\alias{optim.parsimony}
\alias{sankoff}
\alias{fitch}
\alias{PNJ}
\alias{pace}
\title{Parsimony tree.}
\alias{CI}
\alias{RI}
\description{

\code{parsimony} returns the parsimony score of a tree.
\code{optim.parsimony} tries to find the maximum parsimony tree using Nearest Neighbor Interchange (NNI) 
rearrangements. \code{pace} returns a (logical) matrix of the ancestral states
of the root node. \code{CI} and \code{RI} computes the consistency and retention index.  
}
\usage{
parsimony(tree, data, method="fitch", ...)
optim.parsimony(tree, data, method="fitch", cost=NULL, trace=1, ...)
fitch(tree, data, site = "pscore")
sankoff(tree, data, cost = NULL, site = "pscore")
pace(tree, data, ...)
CI(tree, data)
RI(tree, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A object of class phyDat containing (dna) sequences.}
  \item{tree}{ tree to start the nni search from.}
  \item{method}{one of 'fitch' or 'sankoff'.}
  \item{cost}{A cost matrix for the transitions between two states.} 
  \item{site}{return either 'pscore' or 'site' wise parsimony scores. } 
  \item{trace}{defines how much information is printed during optimisation. } 
  \item{...}{Further arguments passed to or from other methods.} 
}
\value{
  \code{parsimony} returns the maximum parsimony score (pscore). 
  \code{optim.parsimony} returns a tree after NNI rearrangements.   
}
\references{Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{ \code{\link{ancestral.pml}}, \code{\link{nni}}, \code{\link{NJ}}, \code{\link{pml}}, \code{\link{getClans}} }
\examples{
data(Laurasiatherian)
dm = dist.logDet(Laurasiatherian)
tree = NJ(dm)
parsimony(tree, Laurasiatherian)
ptree <- optim.parsimony(tree, Laurasiatherian)
}
\keyword{cluster}
