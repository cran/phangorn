\name{pml}
\alias{pml}
\alias{print.pml}
\alias{summary.pml}
\alias{print.summary.pml}
\alias{optim.pml}
\alias{logLik.pml}
\alias{predict.pml}
\alias{vcov.pml}
\title{Likelihood of a tree.}
\description{
\code{optim.pml} computes the likelihood of a phylogenetic tree 
given a sequence alignment and a model. \code{optim.pml} optimizes the 
different model parameters.
}
\usage{
pml(tree, data, bf=rep(1/length(levels),length(levels)), shape=1, 
    k=1, Q=rep(1,length(levels)*(length(levels)-1)/2), levels=attr(data,"levels"), 
	inv=0, g=NULL, w=NULL, eig=NULL, INV=NULL, ll.0=NULL, ...)     
optim.pml(object, optNni=FALSE, optBf=FALSE, optQ=FALSE,
    optInv=FALSE, optGamma=FALSE, optEdge=TRUE, control = list(maxit=10, eps=0.001))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}. }
  \item{data}{The (DNA) alignment.}
  \item{bf}{Base frequencies.}
  \item{shape}{Shape parameter of the gamma distribution.}
  \item{k}{Number of intervalls of the discrete gamma distribution.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.}
  \item{levels}{Factor levels of the data.}
  \item{inv}{Proportion of invariable sites.}
  \item{g}{Rate.}
  \item{w}{Weight of the mixture.}
  \item{eig}{The eigenvalues and eigenvectors of the transition matrix. .}
  \item{INV}{For internal use.}
  \item{object}{An object of class \code{pml}.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{optBf}{Logical value indicating whether base frequencies gets optimized.}
  \item{optQ}{Logical value indicating whether rate matrix gets optimized.
  For dna this means: FALSE=Jukes-Cantor, TRUE=GTR}
  \item{optInv}{Logical value indicating whether proportion of variable size gets optimized.}
  \item{optGamma}{Logical value indicating whether gamma rate parameter gets optimized.}
  \item{optEdge}{Logical value indicating the edge lengths gets optimized.}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{ll.0}{For internal use}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
The input variables w, g, eig, INV are used to speed up computation, 
especially during the optimising used in \code{optim.pml},
but are not intended to be used by the enduser. The toppology search
uses a nearest neighbour interchange (NNI) and is similar to phyML. 
}
\value{
Returns a list of class \code{ll.phylo}
  \item{logLik}{Log likelihood of the tree.}
  \item{siteLik}{Site log likelihoods.}
  \item{root}{likelihood in the root node.}
  \item{weight}{Weight of the site patterns.}
}
\references{ Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maxumum
likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17}, 368--376. 

Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.

Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University Press, Oxford.
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}
\note{For small trees the likelihood seems to be similar to Paup* or PhyML.}
\seealso{ For a different implementation see \code{\link[ape]{mlphylo}}.} 
\examples{
  example(NJ)
# Jukes-Cantor + Gamma + I   
  fitJC <- pml(tree, Laurasiatherian, k=4, inv=.25) 
  summary(fitJC)
# just optimise edge length parameter     
  fitJC <- optim.pml(fitJC) 
  summary(fitJC)
  plot(fitJC$tree)
  
\dontrun{
# optimise parameter     
  fitJC <- optim.pml(fitJC, optNni=TRUE, optGamma=TRUE, optInv=TRUE) 
  summary(fitJC)
# F81 + Gamma + I - model
  fitF81 <- optim.pml(fitJC, optNni=TRUE, optGamma=TRUE, optInv=TRUE, optBf=TRUE)
  summary(fitF81) 
# GTR + Gamma + I - model
  fitGTR <- optim.pml(fitF81, optNni=TRUE, optGamma=TRUE, optInv=TRUE, optBf=TRUE, optQ=TRUE) 
  summary(fitGTR)
}

# 2-state data (RY-coded)  
    
  dat <- as.character(Laurasiatherian)
  # RY-coding
  dat[dat=="a"] <- "r"
  dat[dat=="g"] <- "r"
  dat[dat=="c"] <- "y"
  dat[dat=="t"] <- "y"
  dat <- phyDat(dat, levels=c("r","y"))
  fit2ST <- pml(tree, dat, k=4, inv=.25) 
  fit2ST <- optim.pml(fit2ST,optNni=TRUE, optGamma=TRUE, optInv=TRUE) 
  
  methods(class="pml")  
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
