\name{treedist}
\alias{treedist}
\alias{RF.dist}
%\alias{print.treedist}
\title{ Distances between trees }
\description{
  \code{treedist} computes different tree distance methods and \code{RF.dist} the Robinson-Foulds distance.
}
\usage{
treedist(tree1, tree2)
RF.dist(tree1, tree2, check.labels = FALSE)
}
\arguments{
  \item{tree1}{ A phylogenetic tree. }
  \item{tree2}{ A phylogenetic tree. }
  \item{check.labels}{compares labels of the trees.}
}
\value{
  \code{treedist} returns a vector containing the following tree distance methods
  \item{symmetric.difference}{symmetric.difference or  Robinson-Foulds distance}
  \item{branch.score.difference}{branch.score.difference}
  \item{path.difference}{path.difference}
  \item{weighted.path.difference}{weighted.path.difference}
}
\references{Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparision metrics - some new results}, Syst. Biol.,42(2), 126-141}
\author{ Klaus P. Schliep \email{K.P.Schliep@massey.ac.nz}} 

\keyword{ classif }% at least one, from doc/KEYWORDS
