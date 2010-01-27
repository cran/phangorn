\name{Ancestors}
\alias{Ancestors}
\alias{Children}
\alias{Descendants}
\alias{Siblings}
\title{tree utility function}
\description{
  Functions for describing relationships among
  phylogenetic nodes.
}
\usage{
Ancestors(x, node, type=c("all","parent"))
Children(x, node)
Siblings(x, node, include.self=FALSE)
Descendants(x, node, type=c("tips","children","all"))
}
\arguments{
  \item{x}{a tree (a phylo object).}
  \item{node}{an integer corresponding to a node ID}
  \item{type}{specify whether to return just direct
    ancestor / parents or all }
  \item{include.self}{whether to include self in list of siblings}
}

\value{
   a vector containing the indices of the nodes. 
}
\seealso{these functions are inspired by \code{treewalk} in phylobase package}
\examples{
tree = rtree(10)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()
Ancestors(tree, 1, "all")
Children(tree, 11)
Descendants(tree, 11, "tips")
Siblings(tree, 3)
}
\keyword{misc}
