\name{pml}
\alias{pml}
\alias{summary.pml}
\alias{optim.pml}
\alias{update.pml}
\title{Likelihood of a tree.}
\description{
\code{pml} computes the likelihood of a phylogenetic tree 
given a sequence alignment and a model. \code{optim.pml} optimizes the 
different model parameters.
}
\usage{
pml(tree, data, bf=NULL, Q=NULL, inv=0, k=1, shape=1, rate=1, ...)     
optim.pml(object, optNni=FALSE, optBf=FALSE, optQ=FALSE,
    optInv=FALSE, optGamma=FALSE, optEdge=TRUE, optRate=FALSE, 
    control = pml.control(maxit=10, eps=0.001, trace=TRUE))
}
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}. }
  \item{data}{The (DNA) alignment.}
  \item{bf}{Base frequencies.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.}
  \item{inv}{Proportion of invariable sites.}
  \item{k}{Number of intervalls of the discrete gamma distribution.}
  \item{shape}{Shape parameter of the gamma distribution.}
  \item{rate}{Rate.}
  \item{object}{An object of class \code{pml}.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{optBf}{Logical value indicating whether base frequencies gets optimized.}
  \item{optQ}{Logical value indicating whether rate matrix gets optimized.}
  \item{optInv}{Logical value indicating whether proportion of variable size gets optimized.}
  \item{optGamma}{Logical value indicating whether gamma rate parameter gets optimized.}
  \item{optEdge}{Logical value indicating the edge lengths gets optimized.}
  \item{optRate}{Logical value indicating the overall rate gets optimized.}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
The input variables \code{w0, pl0, rate} are mainly for the use
in mixture or partition models and are only used internally. 
The toppology search uses a nearest neighbour interchange (NNI) 
and the implenatation is similar to phyML. 
Here is a overview how to estimate different phylogenetic models 
with \code{pml}:  
\tabular{lll}{
model \tab optBf \tab optQ \cr
Jukes-Cantor \tab FALSE \tab FALSE \cr
F81 \tab TRUE \tab FALSE \cr
GTR \tab TRUE \tab TRUE
}

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
%\note{For small trees the likelihood seems to be similar to Paup* or PhyML.}
\seealso{ For a different implementation see \code{\link[ape]{mlphylo}}.} 
\examples{
  example(NJ)
# Jukes-Cantor (starting tree from NJ)  
  fitJC <- pml(tree, Laurasiatherian) 
# optimise edge length parameter     
  fitJC <- optim.pml(fitJC)
  summary(fitJC)

\dontrun{    
# search for a better tree using NNI rearrangements     
  fitJC <- optim.pml(fitJC, optNni=TRUE)   
  summary(fitJC)
  plot(fitJC$tree)

# JC + Gamma + I - model
  fitJC_GI <- update(fitJC, k=4, inv=.2)
# optimise shape parameter + proportion of invariant sites     
  fitJC_GI <- optim.pml(fitJC_GI, optGamma=TRUE, optInv=TRUE)
  summary(fitJC_GI) 
# GTR + Gamma + I - model
  fitGTR <- optim.pml(fitJC_GI, optNni=TRUE, optGamma=TRUE, optInv=TRUE, optBf=TRUE, optQ=TRUE) 
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
  
  # show some of the methods available for class pml
  methods(class="pml")  
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
