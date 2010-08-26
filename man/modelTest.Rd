\name{modelTest}
\alias{modelTest}
\title{
ModelTest
}
\description{
Comparison of different substition models 
}
\usage{
modelTest(tree, data, model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, maxit = 3, trace = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{a phylogenetic tree.}
  \item{data}{an object of class phyDat}
  \item{model}{substitution models to compare with each other}
  \item{G}{(discrete) Gamma model}
  \item{I}{invariant sites}
  \item{k}{rate classes}
  \item{control}{A list of parameters for controlling the fitting process.}
}
\value{
a data.frame containing the AIC and BIC for the tested models. 
}
\references{
Posada, D. (2008) jModelTest: Phylogenetic Model Averaging. \emph{Molecular Biology and Evolution} \bold{25}: 1253-1256
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pml}}, \code{\link{anova}}
}
\examples{
\dontrun{    
example(NJ)
modelTest(tree, Laurasiatherian)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

