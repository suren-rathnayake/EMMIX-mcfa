\name{rmcfa}
\alias{rmcfa}
\alias{rmctfa}
\title{
 Random Deviates from mcfa and mctfa Models
}
\description{
Random number generators for \code{mcfa} and \code{mctfa}  models.
}
\usage{
rmcfa(n, model, ...)
rmctfa(n, model, ...)
}
\arguments{
  \item{model}{
  Object of class \code{mcfa}, or \code{mctfa}.
}
  \item{n}{
  Number of sample to generate.
}
  \item{...}{
  Can be used to pass optional arguments to \code{rmvnorm} function of
  the \pkg{mvtnorm} package.
  }
}
\details{
This function uses the \code{rmvnorm} function
in the \pkg{mvtnorm} package to generate samples 
from the mixture components.

Algorithm works by first randomly selecting a component in the
mixture based on the mixing probabilities, and then
drawing a sample from the selection component.

}
\value{
  \item{dat}{
  Matrix with samples drawn in rows.
}

}

\references{
McLachlan GJ, and Peel D (2000).
\emph{Finite Mixture Models}. New York: Wiley.

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.
}
\author{
Geoffrey McLachlan, Suren Rathnayake
}
%%\note{
%%}
\examples{
\dontrun{
 set.seed(1)
 model <- mcfa(iris[, -5], g=3, q=2, nkmeans=1, nrandom=1)
 dat <- rmcfa(n = 10, model = model)
}
}
\keyword{ cluster}
\keyword{ multivariate}
\keyword{ models}
