\name{pmlMix}
\alias{pmlMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Phylogenetic mixture model }
\description{
  Phylogenetic mixture model. 
}
\usage{
pmlMix(formula, fit, m = 2, omega = rep(1/m, m), ...)
}
\arguments{
  \item{formula}{ a formula object (see details).}
  \item{fit}{ an object of class \code{pml}. }
  \item{m}{ number of mixtures. }
  \item{omega}{ mixing weights.}
  \item{\dots}{ Further arguments passed to or from other methods. }
}
\details{
    The \code{formula} object allows to specify which parameter get     
    optimised. The formula is generally of the form \code{edge + bf + Q 
    ~ rate + shape + \dots}, on the left side are the parameters which 
    get optimised over all mixtures, on the right the parameter which
    are optimised specific to each mixture. The parameters available 
    are \code{"nni", "bf", "Q", "inv", "shape", "edge", "rate"}.
    Each parameters can be used only once in the formula. 
    \code{"rate"} and \code{"nni"} are only available for the right side
    of the formula. On the other hand parameters for invariable sites 
    are only allowed on the lefthand side. 
    The convergence of the algorithm is very slow and is likely that the
    algorithm can get stuck in local optima. 
}
\value{
  \code{pmlMix} returns a list with elements
  \item{logLik}{log-likelihood of the fit}
  \item{omega}{mixing weights.} 
  \item{fits}{fits for the final mixtures.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{  \code{\link{pml}},\code{\link{pmlPart}},\code{\link{pmlCluster}} }
\examples{
X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4)
weights <- 1000*exp(fit$site) 
attr(X, "weight") <- weights
fit1 <- update(fit, data=X, k=1)
fit2 <- update(fit, data=X)

(fitMixture <- pmlMix(edge~rate, fit1 , m=4))
(fit2 <- optim.pml(fit2, optGamma=TRUE))


\dontrun{
data(Laurasiatherian)
dm <- dist.logDet(Laurasiatherian)
tree <- NJ(dm)
fit=pml(tree, Laurasiatherian)
fit = optim.pml(fit)

fit2 <- update(fit, k=4)
fit2 <- optim.pml(fit2, optGamma=TRUE)

fitMix = pmlMix(edge ~ rate, fit, m=4)
fitMix
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line