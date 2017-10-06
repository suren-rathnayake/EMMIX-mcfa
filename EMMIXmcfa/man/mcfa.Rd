\name{mcfa}
\alias{mcfa}
\alias{mcfa.default}
\alias{print.emmixmcfa}
\alias{summary.emmixmcfa}
\alias{predict.emmixmcfa}
\alias{plot.emmixmcfa}
\alias{mctfa}
\alias{mctfa.default}
\title{
Mixture of Common Factor Analyzers
}
\description{
Functions for fitting of Mixtures Common Factor Analyzers (MCFA) and
Mixtures of Common \emph{t}-Factor Analyzers (MC\emph{t}FA).
Maximum Likelihood estimates of the model parameters are obtained
using the Expectation--Maximization algorithm.

MCFA adds the following 
restrictions to,
\deqn{
\Sigma_i= A~ \Omega_i ~A^T + D \quad (i=1,\,\ldots,\,g),
}
and
\deqn{
\mu_i = A \xi_i \quad (i=1,\,\ldots,\,g)
}

where \eqn{A} is a \eqn{p \times q} matrix, \eqn{\xi_i} is a
\eqn{q}-dimensional vector, \eqn{\Omega_i} is a 
\eqn{q \times q} positive
definite symmetric matrix, and \eqn{D} is a diagonal 
\eqn{p \times p} matrix.

With this representation,
the component distribution of \eqn{Y_j} is modeled as
\deqn{
Y_j = A ~U_{ij} + e_{ij}
}
with prob. \eqn{\pi_i (i = 1,\,\ldots,\,g)}
for \eqn{j = 1,\,\ldots,\,n}, where the (unobservable) factors
\eqn{U_{i1},\,\ldots,\,U_{in}} are distributed independently
\eqn{N({\xi_i},\,\Omega_i)}, independently of
the \eqn{e_{ij}}, which are distributed independently
\eqn{N(0,\,D)}, where \eqn{D} is a diagonal matrix,
\eqn{(i = 1,\,\ldots,\,g)}.
}
\usage{
mcfa(Y, g, q, ...)
mctfa(Y, g, q, ...)
\method{mcfa}{default}(Y, g, q, itmax = 500, nkmeans = 20, nrandom = 20, 
  tol = 1.e-5, init_clust = NULL, init_para = NULL, 
  init_method = 'randA', conv_measure = 'diff', 
  warn_messages = TRUE, ...) 
\method{mctfa}{default}(Y, g, q, itmax = 500, nkmeans = 20, nrandom = 20, 
  tol = 1.e-5, df_init = rep(30, g), df_update = TRUE, 
  init_clust = NULL, init_para = NULL, init_method = 'randA', 
  conv_measure = 'diff', warn_messages = TRUE, ...)
\method{print}{emmixmcfa}(x, ...)
\method{summary}{emmixmcfa}(object, ...)
\method{plot}{emmixmcfa}(x, ...)
\method{predict}{emmixmcfa}(object, Y, ...)
}
\arguments{
  \item{Y}{
  A matrix or a data frame of which rows correspond to
  observations and columns to variables.
}
\item{x, object}{
   An object of class \code{mcfa} or \code{mctfa}.
}
  \item{g}{
Number of components.
}
  \item{q}{
Number of factors.
}
  \item{itmax}{
Maximum number of EM iterations.
}
  \item{nkmeans}{
  The number of times the k-means algorithm to be used in partition 
  the data into \code{g} groups. These groupings are then used in 
  initializing the parameters for the EM algorithm.
}
  \item{nrandom}{
  The number of random \code{g}-group partitions for the data to be used 
  initializing the EM algorithm.  
}
  \item{tol}{
  The EM algorithm terminates if the measure of convergence falls below 
  this value. 
}
  \item{init_clust}{
  A vector or matrix consisting of partition of samples to be used
  in the EM algorithm. For matrix of partitions, columns must corresponds 
  individual partitions of the data. Optional.
}
  \item{init_para}{
  A list containing model parameters to be used as initial 
  parameter estimates for the EM algorithm. Optional.
}
  \item{init_method}{
  To determine how the initial parameter values are computed. See Details.
}  
  \item{conv_measure}{
  The default \code{'diff'} stops the EM iterations if 
  |\eqn{l^{(k+1)}} - \eqn{l^{(k)}}| < \code{tol} where
  \eqn{l^{(k)}} is the log-likelihood at the \eqn{k}th EM iteration. 
  If \code{'ratio'}, then the convergence of the EM steps is measured 
  using the |(\eqn{l^{(k+1)}} - \eqn{l^{(k)}})/\eqn{l^{(k+1)}}|. 
}
  \item{df_init}{
  Initial values of the degree of freedom parameters for \code{mctfa}. 
}
  \item{df_update}{
  If \code{df_update = TRUE} (default), then the degree of freedom parameters
  values will be updated during the EM iterations. 
  Otherwise, if \code{df_update = FALSE}, they will be fixed at the initial 
  values specified in \code{df_init}.
}
  \item{warn_messages}{
  If \code{warn_messages = TRUE} (default), the output would
  include error messages for instances, if any, where the
  model fitting function failed to provide estimates of parameters.
  Otherwise the messages will not be stored.
}
 \item{\dots}{
Not used.
}
}
\details{
With the default \code{init_method = "randA"}, initialization of
the parameters is done by using the procedure in
Baek et al. (2010) where initial values for elements of 
\eqn{A} are drawn from the \eqn{N(0, 1)} distribution. 
With this method, the \strong{columns of the 
data need to be normalized}.
}
\value{
Object of class \code{c("emmixmcfa", "mcfa")} or \code{c("emmixmcfa",
"mctfa")} containing the fitted model parameters is returned.
Details of the components are as fellows:
 \item{g}{
 Number of mixture components.
}
 \item{q}{
  Number of factors.
}
 \item{pivec}{
  Mixing proportions of the components.
}
  \item{A}{
  Loading matrix. Size \eqn{p \times q}.
}
  \item{xi}{
  Matrix containing factor means for components in columns. 
  Size \eqn{q \times g}.
}
  \item{omega}{
  Array containing factor covariance matrices for components.
  Size \eqn{q \times q \times g}.
}
  \item{D}{
  Error covariance matrix. Size \eqn{p \times p.}
}
  \item{U}{
    Estimated conditional expected component scores of the 
    unobservable factors given the data and the component membership.
    Size is Size \eqn{n \times q \times g}.
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
  \item{clust}{Cluster labels.
} 
  \item{tau}{Posterior probabilities.
}
  \item{logL}{Log-likelihood of the model.
}
  \item{BIC}{Bayesian Information Criteria.
}
  \item{warn_msg}{Description of error messages, if any.
}
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
\author{
Suren Rathnayake, Jangsun Baek, Geoffrey McLachlan
}
%%\note{
%%  ~~further notes~~
%%}
%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
mcfa_fit <- mcfa(iris[, -5], g = 3, q = 3, 
              itmax = 250, nkmeans = 5, nrandom = 5, tol = 1.e-5)

plot(mcfa_fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{multivariate}
\keyword{models}% __ONLY ONE__ keyword per line
