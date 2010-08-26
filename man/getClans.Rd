\name{getClans}
\alias{getClans}
\alias{getClips}
\alias{getSlices}
\alias{getDiversity}
\title{
Clans and Clips
}
\description{
Functions for clanistics, which compute clan, slices and clips for unrooted trees. 
}
\usage{
getClans(tree)
getClips(tree, all=TRUE)
getSlices(tree)
getDiversity(tree, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{An object of class phylo or multiPhylo (getDiversity).}
\item{all}{A logical, return all or just the largest clip. }
\item{x}{An object of class phyDat. }
}
\details{
Every edge in an unrooted tree defines two complemantary clans.

Slices are defined by a pair of splits or tripartitions, which are not 
clans. The number of distinguishable slices for a binary tree with 
n tips is \eqn{2n^2 - 10n + 12}. 
   
A clip is a group of tips for which the maximum of the pairwise 
distance within these group is smaller than the distince between 
any member of the clip and any other tip.

getDiversity computes how well a vector of characters (natives,  intruders) 
fits on a tree. getDiversity computes a normalized version of the shannon diversity 
index (E-value): \eqn{E = -\sum_{i}(p_i * log(p_i)) /log(\sum_{i}(p_i))}{E = -sum(p * log(p)) /log(sum(p))},
where \eqn{p_i} are the sizes of the homogeneous clans of intruders. 
If the data supplied by x can be separated by an edge of the tree than E-value will be zero.
}
\value{
getClans, getSlices and getClips return a matrix of bipartitions.
getDiversity a matrix of the divergence index for all tree.
}

\references{ 
Lapointe, F.-J., Lopez, P., Boucher, Y., Koenig, J., Bapteste, E. (2010) Clanistics: a multi-level perspective for harvesting unrooted gene trees. \emph{Trends in Microbiology} (in press) 

Wilkinson, M., McInerney, J.O., Hirt, R.P., Foster, P.G., Embley, T.M. (2007) Of clades and clans: terms for phylogenetic relationships in unrooted trees. \emph{Trends in Ecology and Evolution} 22: 114-115 
}
\author{
Klaus Schliep \email{klaus.schliep@snv.jussieu.fr} 

Francois-Joseph Lapointe \email{francois-joseph.lapointe@umontreal.ca}
}

\seealso{
\code{\link{parsimony}}, \code{\link{phyDat}}
}


\examples{
set.seed(1)
tree = rtree(10)
getClans(tree)
getClips(tree, all=TRUE)
getSlices(tree)
X = matrix(c(1,1,1,1,1,2,2,2,2,2), ncol=1, dimnames=list(paste('t',1:10, sep=""),NULL))
x = phyDat(X, type = "USER", levels = c(1,2))
plot(tree, "u", tip.color = X[tree$tip,])  # red intruder 
getDiversity(tree, x)
}
\keyword{ cluster }

