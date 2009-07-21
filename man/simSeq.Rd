\name{simSeq}
\alias{simSeq}
\title{ Simulate sequences. }
\description{
Simulate sequences for a given evolutionary tree.
}
\usage{
simSeq(tree, l=1000, Q=NULL, bf=NULL, rootseq=NULL, type="DNA",
    model="", levels=NULL, rate=1, ancestral=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{ a phylogenetic tree \code{tree}, an object of class \code{phylo}. }
\item{l}{  length of the sequence to simulate. }
\item{Q}{ the rate matrix.  }
\item{bf}{ base frequences. }
%  \item{ancestral}{ write out the ancestral \code{ancestral} sequnces? }
\item{rootseq}{a vector of length l containing the root sequence, 
other root sequence is randomly generated.}
\item{type}{Type of sequences ("DNA", "AA" or "USER").}
\item{model}{Amino acid models: one of "WAG", "JTT", "Dayhoff" or "LG"}
\item{levels}{ \code{levels} takes a character vector of the different bases,
default is for nucleotide sequences, only used when type = "USER".}
\item{rate}{rate. }
\item{ancestral}{Return ancestral sequences?}
}
\details{
\code{simSeq} simulates sequence alingments. So far rate variation 
is not yet implemented, but one can combine different alignments 
having their own rate.  
In fact it is possible to generate DNA, RNA, amino acid, or 0/1. 
% It models using a nonhomogeneous markov process.   
% It is possible to generate sequences with changin base composition 
% Barry and Hartigan to model complicated models (e.g. covarion). 
}
\value{
\code{simSeq} returns an object of class phyDat.
}
\author{ Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{  \code{\link[phangorn]{phyDat}}}
\examples{
tree = rtree(5)
plot(tree)
nodelabels()

# Example for simple DNA alignment
data = simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6)
as.character(data)

# Example to simulate discrete Gamma rate variation
rates = phangorn:::discrete.gamma(1,4)
data1 = simSeq(tree, l = 100, type="AA", model="WAG", rates[1])
data2 = simSeq(tree, l = 100, type="AA", model="WAG", rates[2])
data3 = simSeq(tree, l = 100, type="AA", model="WAG", rates[3])
data4 = simSeq(tree, l = 100, type="AA", model="WAG", rates[4])
data <- c(data1,data2, data3, data4)

write.phyDat(data, file="temp.dat", format="sequential",nbcol = -1, colsep = "")
unlink("temp.dat") 
}
\keyword{ cluster }% at least one, from doc/KEYWORDS

