\name{rmix}
\alias{rmix}
\title{
 Random Deviates from EMMIXmcfa Models
}
\description{
Random number generator for \code{EMMIXmcfa} models.
}
\usage{
rmix(n, model, ...)
}
\arguments{
  \item{model}{
  Model of class \code{mcfa, mctfa, mfa}, or \code{mtfa}.
}
  \item{n}{
  Number of sample to generate.
}
  \item{...}{
  Not used.
  }
}
\details{
This function uses the \code{rdemmix2} function
in the \pkg{EMMIX} package to generate samples 
from the mixture components.

Algorithm works by first drawing a component based on
the mixture proprotion in the model, and then drawing
a sample from the component distribution.

}
\value{
  \item{dat}{
  Matrix with samples drawn in rows.
}

}

\references{
\url{https://people.smp.uq.edu.au/GeoffMcLachlan/mix_soft/EMMIX_R/}

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
 dat <- rmix(n = 10, model = model)
}
}
\keyword{ cluster}
\keyword{ multivariate}
\keyword{ models}
