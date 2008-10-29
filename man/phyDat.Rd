\name{phyDat}
\alias{phyDat}
\alias{as.phyDat.DNAbin}
\alias{as.data.frame.phyDat}
\alias{as.character.phyDat}
\alias{allSitePattern}
\title{Conversion amoung Sequence Formats}
\description{
These functions transform several DNA formats into the \code{phyDat} format. 
\code{allSitePattern} generates an alignment of all possible site patterns.
}
\usage{
phyDat(data, type = "DNA", levels = NULL, return.index=FALSE, ...) 
\method{as.phyDat}{DNAbin}(x, ...)
\method{as.character}{phyDat}(x, ...)
\method{as.data.frame}{phyDat}(x, ...)
allSitePattern(n, levels=c("a","c","g","t"), names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object containing sequences.}
  \item{x}{An object containing sequences.}
  \item{type}{Type of sequences ("DNA", "AA" or "USER").}
  \item{levels}{Level atttributes.}
  \item{return.index}{If TRUE returns a index of the site patterns.}
  \item{n}{Number of sequences.}
  \item{names}{Names of sequences.}
  \item{...}{further arguments passed to or from other methods.}
}
\details{
If \code{type} "USER" a vector has to be give to \code{levels}.
For example c("a", "c", "g", "t", "-") would create a data object that  
can be used in phylogenetic analysis with gaps as fifth state.  
\code{allSitePattern} returns all possible site patternsis and can
be usefull in simulation studies.  
}
\value{
The functions return an object of class \code{phyDat}. 
You may import data with \code{\link[ape]{read.dna}} or \code{\link[ape]{read.nexus.data}}. 
}

\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{ \code{\link[ape]{as.DNAbin}}, \code{\link[ape]{read.dna}} and \code{\link[ape]{read.nexus.data}} and 
the example of \code{\link[phangorn]{pmlMix}} for the use of \code{allSitePattern}}
\examples{
data(Laurasiatherian)
class(Laurasiatherian)
Laurasiatherian
# transform into old ape format
LauraChar <- as.character(Laurasiatherian)
# and back 
Laura <- phyDat(LauraChar, return.index=TRUE)
all.equal(Laurasiatherian, Laura)
allSitePattern(5)
}
\keyword{cluster}
