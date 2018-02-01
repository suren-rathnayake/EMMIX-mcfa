\name{EMMIXmcfa-package}
\alias{EMMIXmcfa-package}
\alias{emmixmcfa-package}
\alias{EMMIXmcfa}
\alias{emmixmcfa}
\docType{package}
\title{
Mixture Models with Common Factor Loadings
}
\description{
Fits a mixture of factor analyzers with a common component matrix for the 
factor loadings before the transformation of the latent factors to be 
white noise. It is designed specifically for the task of displaying the 
observed data points in a
lower (q-dimensional) space, where q is the number of factors adopted in the 
factor-analytic representation of the observed vector.

It also provides a greater reduction in the number of parameters in the model.
Component distributions can either be from the family of multivariate normals 
or from the family of multivariate \emph{t}-distributions. 
Maximum likelihood estimators of
model parameters are obtained using the Expectation-Maximization algorithm.
}
\details{
\tabular{ll}{
Package: \tab EMMIXmcfa\cr
Type: \tab Package\cr
Version: \tab 2.0.8\cr
Date: \tab 2017-10-18\cr
License: \tab GPL\cr
}
%%~~ An overview of how to use the package, including the most important ~~
%%~~ functions ~~
}
\author{
%Who wrote it
Suren Rathnayake, Jangsun Baek, Geoffrey McLachlan 
}
\references{
Baek J, McLachlan GJ, and Flack LK (2010). Mixtures of factor analyzers 
with common factor loadings: applications to the clustering and visualisation 
of high-dimensional data. \emph{IEEE Transactions on Pattern Analysis and 
Machine Intelligence} \strong{32}, 2089--2097.

Baek J, and McLachlan GJ (2011). Mixtures of common \emph{t}-factor analyzers 
for clustering highdimensional microarray data. 
\emph{Bioinformatics} \strong{27}, 1269--1276.

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.
}
\keyword{ package }
\keyword{ clustering }
\keyword{ multivariate }
\keyword{ model }
\seealso{
\code{mcfa}
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
\dontrun{
set.seed(1)
Y <- scale(iris[, -5])
mcfa_model <- mcfa(Y, g = 3, q = 3)
mctfa_model <- mctfa(Y, g = 3, q = 3)
}
}
