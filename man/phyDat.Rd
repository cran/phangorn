\name{phyDat}
\alias{phyDat}
\alias{as.phyDat.DNAbin}
\alias{phyDat.DNA}
\alias{print.phyDat}
\alias{as.data.frame.phyDat}
\alias{as.character.phyDat}
\title{Conversion amoung Sequence Formats}
\description{
These functions transform several DNA formats into the \code{phyDat} format. 
}
\usage{
phyDat(data, levels, return.index=FALSE)
phyDat.DNA(data, return.index=FALSE) 
as.phyDat.DNAbin(data, ...)
\method{as.character}{phyDat}(x, ...)
\method{as.data.frame}{phyDat}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object containing sequences.}
  \item{levels}{Level atttributes.}
  \item{return.index}{If TRUE returns a index of the site patterns.}
  \item{x}{An object of class \code{phyDat}.}
  \item{...}{further arguments passed to or from other methods.}
}
\value{
The functions return an object of class \code{phyDat}. 
You may import data with \code{\link[ape]{read.dna}} or \code{\link[ape]{read.nexus.data}}. 
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}

\seealso{ \code{\link[ape]{as.DNAbin}}, \code{\link[ape]{read.dna}} and \code{\link[ape]{read.nexus.data}} }
\examples{
data(Laurasiatherian)
class(Laurasiatherian)
Laurasiatherian
# transform into old ape format
LauraChar <- as.character(Laurasiatherian)
# and back 
Laura <- phyDat.DNA(LauraChar, return.index=TRUE)
all.equal(Laurasiatherian, Laura)

}
\keyword{cluster}
