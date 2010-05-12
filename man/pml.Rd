\name{pml}
\alias{pml}
\alias{optim.pml}
\alias{pml.control}
\title{Likelihood of a tree.}
\description{
\code{pml} computes the likelihood of a phylogenetic tree 
given a sequence alignment and a model. \code{optim.pml} optimizes the 
different model parameters.
}
\usage{
pml(tree, data, bf=NULL, Q=NULL, inv=0, k=1, shape=1, rate=1, model="", ...)     
optim.pml(object, optNni=FALSE, optBf=FALSE, optQ=FALSE,
    optInv=FALSE, optGamma=FALSE, optEdge=TRUE, optRate=FALSE, 
    control = pml.control(maxit=10, eps=0.001, trace=TRUE), model = NULL, subs = NULL, ...)  
}
\arguments{
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}. }
  \item{data}{The (DNA) alignment.}
  \item{bf}{Base frequencies.}
  \item{Q}{A vector containing the lower triangular part of the rate matrix.}
  \item{inv}{Proportion of invariable sites.}
  \item{k}{Number of intervals of the discrete gamma distribution.}
  \item{shape}{Shape parameter of the gamma distribution.}
  \item{rate}{Rate.}
  \item{model}{Amino acid models: one of "WAG", "JTT", "Dayhoff" or "LG" or nucleotide model}
  \item{object}{An object of class \code{pml}.}
  \item{optNni}{Logical value indicating whether toplogy gets optimized (NNI).}
  \item{optBf}{Logical value indicating whether base frequencies gets optimized.}
  \item{optQ}{Logical value indicating whether rate matrix gets optimized.}
  \item{optInv}{Logical value indicating whether proportion of variable size gets optimized.}
  \item{optGamma}{Logical value indicating whether gamma rate parameter gets optimized.}
  \item{optEdge}{Logical value indicating the edge lengths gets optimized.}
  \item{optRate}{Logical value indicating the overall rate gets optimized.}
  \item{control}{A list of parameters for controlling the fitting process.}
  \item{subs}{a (integer) vector same length as Q to specify the optimization of Q}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
The topology search uses a nearest neighbor interchange (NNI) 
and the implementation is similar to phyML. 
The option model in pml is only used for amino acid models. 
The option model defines the nucleotide model which is getting optmised, 
all models which are included in modeltest can be chosen. Setting this option 
(e.g. "K81" or "GTR") overrules options optBf and optQ.  
Here is a overview how to estimate different phylogenetic models 
with \code{pml}:  
\tabular{lll}{
model \tab optBf \tab optQ \cr
Jukes-Cantor \tab FALSE \tab FALSE \cr
F81 \tab TRUE \tab FALSE \cr
symmetric \tab FALSE \tab TRUE \cr
GTR \tab TRUE \tab TRUE
}
Via model in optim.pml the following nucleotide models can be specified:  
JC, F81, K80, HKY, TrNe, TrN, TPM1, K81, TPM1u, TPM2, TPM2u, TPM3, TPM3u, 
TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3, TVMe, TVM, SYM and GTR. 
For how these s models are specified see Posada (2008).

So far 4 amino acid models are supported ("WAG", "JTT", "Dayhoff" and "LG"). 
}
\value{
Returns a list of class \code{ll.phylo}
  \item{logLik}{Log likelihood of the tree.}
  \item{siteLik}{Site log likelihoods.}
  \item{root}{likelihood in the root node.}
  \item{weight}{Weight of the site patterns.}
}
\references{ 
Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maxumum
likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17}, 368--376. 

Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer Associates, Sunderland.

Yang, Z. (2006). \emph{Computational Molecular evolution}. Oxford University Press, Oxford.

Whelan, S. and Goldman, N. (2001) A general empirical model of 
protein evolution derived from multiple protein families using 
a maximum-likelihood approach. \emph{Molecular Biology and Evolution},  \bold{18}, 691--699

Le, S.Q. and Gascuel, O. (2008) LG: An Improved, General 
Amino-Acid Replacement Matrix \emph{Molecular Biology and Evolution} \bold{25(7)}, 1307--1320

Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. \emph{Molecular Biology and Evolution} \bold{25}: 1253-1256 
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
% \note{For small trees the likelihood seems to be very similar to Paup* or PhyML.}
% \seealso{ For a different implementation see \code{\link[ape]{mlphylo}}.} 
\examples{
  example(NJ)
# Jukes-Cantor (starting tree from NJ)  
  fitJC <- pml(tree, Laurasiatherian)  
# optimize edge length parameter     
  fitJC <- optim.pml(fitJC)
  fitJC 
  
\dontrun{    
# search for a better tree using NNI rearrangements     
  fitJC <- optim.pml(fitJC, optNni=TRUE)
  fitJC   
  plot(fitJC$tree)

# JC + Gamma + I - model
  fitJC_GI <- update(fitJC, k=4, inv=.2)
# optimize shape parameter + proportion of invariant sites     
  fitJC_GI <- optim.pml(fitJC_GI, optGamma=TRUE, optInv=TRUE)
# GTR + Gamma + I - model
  fitGTR <- optim.pml(fitJC_GI, optNni=TRUE, optGamma=TRUE, optInv=TRUE, optBf=TRUE, optQ=TRUE) 
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
  fit2ST
  # show some of the methods available for class pml
  methods(class="pml")  
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
