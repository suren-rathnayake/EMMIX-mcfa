\name{factor_scores}
\alias{factor_scores}
\alias{factor_scores.mcfa}
\alias{factor_scores.mctfa}
\title{
 Computes Factor Scores.
}
\description{
This function computes factor scores given a data set and
a \code{EMMIXmcfa} model.
}
\usage{
factor_scores(model, Y, ...)
\method{factor_scores}{mcfa}(model, Y, tau = NULL, clust= NULL, ...)
\method{factor_scores}{mctfa}(model, Y, tau = NULL, clust= NULL, ...)
}
\arguments{
  \item{model}{
  Model of class \code{"mcfa"} or \code{"mctfa"}.
}
  \item{Y}{
  Data matrix with variables in columns in the same order as used in
  model estimation.
}
  \item{tau}{
  Optional. Posterior probabilities of belonging to the components
  in the mixture model. If not provided, they will be estimated. 
}
 \item{clust}{
 Optional. Indicators of belonging to the components. 
 If not provided, will be estimated using \code{tau}.
}
 \item{\dots}{
   Not used.
}     
}
\details{
Factor scores can be used in visualization of the data 
in the factor space.
}
\value{
  \item{U}{
    Estimated conditional expected component scores of the 
    unobservable factors given the data and the component membership.
    Size is \eqn{n \times q \times g}, where \code{n} is the number of sample,
    \code{q} is the number of factors and \code{g} is the number components.
}
  \item{Fmat}{
   Means of the estimated conditional expected factors scores over
   estimated posterior distributions. Size \eqn{n \times q}.
}
 \item{UC}{
 Alternative estimate of \code{Fmat} where the posterior probabilities
 for each sample are replaced by component indicator vectors 
 which contain one in the element corresponding to the highest posterior 
 probability while others zero.  Size \eqn{n \times q}.
}

}

\references{
McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.

McLachlan GJ, and Peel D (2000). \emph{Finite Mixture Models}. 
New York: Wiley.
}
\author{
Geoffrey McLachlan, Suren Rathnayake, Jungsun Baek
}
%%\note{
%%}
\examples{

# Fit a MCFA model to a subset 
set.seed(1)
samp_size <- dim(iris)[1]
sel_subset <- sample(1 : samp_size, 75)
model <- mcfa(iris[sel_subset, -5], g=3, q=2, nkmeans=1, nrandom=0)

# plot the data points in the factor space
plot(model)

# Allocating new samples to the clusters
Y <- iris[-c(sel_subset), -5]
Y <- as.matrix(Y)
clust <- predict(model, Y)

fa_scores <- factor_scores(model, Y)
# Visualizing new data in factor space
plot_factors(fa_scores, type="Fmat", clust=clust)
}
\keyword{ cluster}
\keyword{ multivariate}
\keyword{ models}
