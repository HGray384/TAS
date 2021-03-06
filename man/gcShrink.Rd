% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcShrink.R
\name{gcShrink}
\alias{gcShrink}
\title{Bayesian Gaussian conjugate (GC) single target linear shrinkage covariance estimator}
\usage{
gcShrink(
  X,
  target = "none",
  var = 2,
  cor = 1,
  alpha = seq(0.01, 0.99, 0.01),
  plots = TRUE,
  weighted = FALSE,
  ext.data = FALSE
)
}
\arguments{
\item{X}{\code{matrix} -- data matrix with variables in rows and 
observations in columns. This method performs best when there are more 
variables than observations.}

\item{target}{\code{character} or \code{matrix} -- if "none" then a default
target specified by \code{var} and \code{cor} will be used for shrinkage. If
\code{matrix} then this will be used as the target for shrinkage. The
target must be a real symmetric positive definite matrix.}

\item{var}{\code{numeric} -- c(1, 2, 3) variance structure for the target matrix.
1 sets all variances equal to 1. 2 sets all variances equal to their sample mean.
3. sets all variances to their sample values.}

\item{cor}{\code{numeric} -- c(1, 2, 3) correlation structure for the 
target matrix. 1 sets the correlations to 0. 2 sets the correlations equal 
to their sample mean. 3 sets the correlations equals to an autocorrelation
structure with parameter equal to the sample mean.}

\item{alpha}{\code{list} -- the grid of shrinkage intensities in (0, 1) to be used. 
Recommended to be an equidistant grid that covers the whole interval. A short
comparison of estimation accuracy versus granularity is provided in ...}

\item{plots}{\code{logical} -- if TRUE then plots the log-marginal likelihood
for each value of alpha with the value of alpha that maximises this highlighted.}

\item{weighted}{\code{logical} -- if TRUE then average over all values of 
alpha and their respective marginal likelihood value as in Gray et al (submitted). 
If FALSE then only use the value of alpha that maximises the log-marginal likelihood
as in Hannart and Naveau (2014).}

\item{ext.data}{\code{matrix} -- an external data matrix used a surrogate to estimate 
the parameters in the default target set for X. Never recommended 
unless there is a belief that ext.data is informative of the covariances of
 X.}
}
\value{
{\code{list} --
\describe{
\item{sigmahat}{\code{matrix} -- the estimated covariance matrix.}
\item{optimalpha}{\code{numeric} -- the value of alpha that maximises the 
log-marginal likelihood.}
\item{target}{\code{matrix} -- the target matrix used for shrinkage.}
\item{logmarg}{\code{numeric} -- the values of the log marginal 
likelihood for each (target, alpha) pair. }
}}
}
\description{
Implements a Bayesian Gaussian conjugate (GC) single target linear shrinkage 
covariance estimator as in Gray et al. (2018)  and Hannart and Naveau (2014). 
It is most useful when the observed data is high-dimensional 
(more variables than observations) and allows a user-specified target matrix.
}
\examples{
set.seed(102)
X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
t1 <- gcShrink(X, var=1, cor=1) # apply shrinkage and view likelihood for T1
t2 <- gcShrink(X, var=2, cor=2) # apply shrinkage and view likelihood for T2
norm(t1$sigmahat-diag(10), type="F") # calculate loss
norm(t2$sigmahat-diag(10), type="F") # calculate loss
# one target clearly better but how to choose this a priori?
}
\references{
Gray, H., Leday, G.G., Vallejos, C.A. and Richardson, S.,
 2018. Shrinkage estimation of large covariance matrices using 
 multiple shrinkage targets. \href{https://arxiv.org/abs/1809.08024}{arXiv preprint}.

Hannart, A. and Naveau, P., 2014. Estimating high dimensional 
covariance matrices: A new look at the Gaussian conjugate framework. 
Journal of Multivariate Analysis, 131, pp.149-162. 
\href{http://dx.doi.org/10.1016/j.jmva.2014.06.001}{doi}.
}
