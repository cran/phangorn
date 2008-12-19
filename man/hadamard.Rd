\name{hadamard}
\alias{hadamard}
\alias{fhm}
\alias{h4st}
\alias{h2st}
\alias{write.nexus.splits}
\title{Hadamard Matrices and Fast Hadamard Multiplication}
\description{
A collection of functions to perform Hadamard conjugation. 
%Hv of a Hadamard matrix H with a vector v using fast Hadamard multiplication.
}
\usage{
hadamard(x)
fhm(v)
h2st(obj, levels=c("r","y"))
h4st(obj, levels = c("a","c","g","t"))
write.nexus.splits(obj, file="")
}
\arguments{
  \item{x}{ a vector of length \eqn{2^n}, where n is an integer. }
  \item{v}{ a vector of length \eqn{2^n}, where n is an integer. }
  \item{obj}{ a data.frame or character matrix, typical a sequence alignment.}
  \item{levels}{ levels of the sequences.}
  \item{file}{ a file name.}
}

\details{
\code{h2st} and \code{h4st} perform Hadamard conjugation 
for 2-state (binary, RY-coded) or 4-state (DNA/RNA) data.
\code{write.nexus.splits} writes splits returned from  \code{h2st} or  
\code{\link[phangorn]{distanceHadamard}} to a nexus file, which can be processed by 
Spectronet or Splitstree.
}

\value{
\code{hadamard} returns a Hadamard matrix. \code{fhm} returns the 
fast Hadamard multiplication. 
}
\references{Hendy, M.D. (1989). The relationship between simple evolutionary tree
models and observable sequence data. \emph{Systematic Zoology}, \bold{38} 310--321.

Hendy, M. D. and Penny, D. (1993). Spectral Analysis of Phylogenetic
Data. \emph{Journal of Classification}, \bold{10}, 5--24.

Hendy, M. D. (2005). Hadamard conjugation: an analytical tool for phylogenetics.
In O. Gascuel, editor, \emph{Mathematics of evolution and phylogeny}, Oxford University Press, Oxford

Waddell P. J. (1995). Statistical methods of phylogenetic analysis: Including hadamard conjugation, 
LogDet transforms, and maximum likelihood. \emph{PhD thesis}.
}
\author{Klaus Schliep \email{K.P.Schliep@massey.ac.nz}}
\seealso{\code{\link[phangorn]{distanceHadamard}}}

\examples{
H = hadamard(3)
v = 1:8
H%*%v
fhm(v)
\dontrun{
data(Laurasiatherian)
dat = as.character(Laurasiatherian)[1:10,]
# RY-coding
dat[dat=="a"] = "r"
dat[dat=="g"] = "r"
dat[dat=="c"] = "y"
dat[dat=="t"] = "y"

fit = h2st(dat)
write.nexus.splits(fit[fit$edges>1e-4, ][-1,], file = "test.nxs")
# read this file into Spectronet or Splitstree to show the network
unlink("test.nxs") 
}

}
\keyword{ cluster }% at least one, from doc/KEYWORDS

