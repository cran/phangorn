\name{phyDat}
\alias{phyDat}
\alias{as.phyDat.DNAbin}
%\alias{as.data.frame.phyDat}
\alias{as.character.phyDat}
\alias{as.DNAbin.phyDat}
\alias{read.phyDat}
\alias{write.phyDat}
\alias{allSitePattern}
\alias{as.phyDat}
\alias{[.phyDat}
\title{Conversion amoung Sequence Formats}
\description{
These functions transform several DNA formats into the \code{phyDat} format. 
\code{allSitePattern} generates an alignment of all possible site patterns.
}
\usage{
phyDat(data, type = "DNA", levels = NULL, return.index=TRUE, ...) 
read.phyDat(file, format="phylip", type="DNA", ...)
write.phyDat(x, file, format="phylip",...)
\method{as.phyDat}{DNAbin}(x, ...)
\method{as.character}{phyDat}(x, ...)
%\method{as.data.frame}{phyDat}(x, ...)
\method{as.DNAbin}{phyDat}(x, ...)
allSitePattern(n, levels=c("a","c","g","t"), names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object containing sequences.}
  \item{x}{An object containing sequences.}
  \item{type}{Type of sequences ("DNA", "AA" or "USER").}
  \item{levels}{Level atttributes.}
  \item{return.index}{If TRUE returns a index of the site patterns.}
  \item{file}{A file name.}
  \item{format}{File format of the sequence alignment (see details).}
  \item{n}{Number of sequences.}
  \item{names}{Names of sequences.}
  \item{...}{further arguments passed to or from other methods.}
}
\details{
If \code{type} "USER" a vector has to be give to \code{levels}.
For example c("a", "c", "g", "t", "-") would create a data object that  
can be used in phylogenetic analysis with gaps as fifth state.  
\code{allSitePattern} returns all possible site patterns and can
be usefull in simulation studies.  
\code{write.phyDat} calls the function write.dna or write.nexus.data and  
\code{read.phyDat} calls the function read.dna, read.aa or read.nexus.data
see for more details over there.
  
You may import data directly with \code{\link[ape]{read.dna}} or \code{\link[ape]{read.nexus.data}}
and convert the data to class phyDat. 

The generic function 'c' can be used to to combine sequences and unique to get 
all unique sequences or haplotides. 
}
\value{
The functions return an object of class \code{phyDat}. 
}

\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

\seealso{ \code{\link{DNAbin}}, \code{\link{as.DNAbin}}, \code{\link{read.dna}},  \code{\link{read.aa}} and \code{\link{read.nexus.data}} and 
the example of \code{\link{pmlMix}} for the use of \code{allSitePattern}}
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
