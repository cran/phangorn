\name{ancestral.pml}
\alias{ancestral.pml}
\alias{ancestral.pars}
\alias{pace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancestral character reconstruction. 
}
\description{
Marginal reconstruction of the ancestral character states (of the root).
}
\usage{
ancestral.pml(object, type = "ml", ...)
ancestral.pars(tree, data, eps= 1e-08, ...)
pace(tree, data, eps= 1e-08, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class pml}
  \item{type}{either "ml" or "bayes"}
  \item{tree}{a tree, i.e. an object of class pml}
  \item{data}{an object of class phyDat} 
  \item{eps}{a small value to prevent rounding errors}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
A matrix containing the 
}
\references{
Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\seealso{
\code{pml}, \code{parsimony}, \code{ace}, \code{root}
}
\examples{
example(NJ)
fit = pml(tree, Laurasiatherian)
anc.ml = ancestral.pml(fit, type = "ml")
anc.p = ancestral.pars(tree, Laurasiatherian)
\dontrun{
require(seqLogo)
seqLogo(t(anc.ml[1:10,]), ic.scale=FALSE)
seqLogo(t(anc.p[1:10,]), ic.scale=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
