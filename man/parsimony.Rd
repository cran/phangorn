\name{parsimony}
\alias{parsimony}
\alias{optim.parsimony}
\alias{sankoff}
\alias{fitch}
\title{Parsimony tree.}
\description{

\code{parsimony} returns the parsimony score of a tree.
\code{optim.parsimony} tries to find the mamximum parsimony tree using NNI 
rearrangements.  
}
\usage{
parsimony(tree, data, method="sankoff", ...)
optim.parsimony(tree, data, cost=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A (dna) sequences.}
  \item{tree}{ tree to start the nni search from.}
  \item{method}{one of 'fitch' or 'sankoff'.}
  \item{cost}{A cost matrix for the transitions between two states.}  
  \item{...}{Further arguments passed to or from other methods.} 
}
\value{
  \code{parsimony} returns the maximum parsimony score (pscore). 
  \code{optim.parsimony} returns a tree after NNI rearrangements.   
}
\references{Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}
\seealso{ \code{\link{sankoff}}, \code{\link{NJ}}}
\examples{
data(Laurasiatherian)
dm = dist.dna(as.DNAbin(as.character(Laurasiatherian)))
tree = NJ(dm)
parsimony(tree, Laurasiatherian)
ptree <- optim.parsimony(tree, Laurasiatherian)
}
\keyword{cluster}
