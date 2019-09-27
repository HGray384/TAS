\name{TAS-package}
\alias{TAS-package}
\alias{TAS}
\docType{package}
\title{
  Target-Averaged Linear Shrinkage estimation
}
\description{
  Conjugate Bayesian covariance matrix estimation using linear shrinkage with multiple target matrices (Gray et al., 2018). Most useful in high-dimensional data settings, where the number of variables is greater than the number of samples.
}
\details{
  This package contains functions for covariance estimation using a conjugate
  Bayesian model. Whilst the main functionality of the package is for multiple
  target linear shrinkage estimation, we also provide functionality for the 
  single target analogue (Hannart and Naveau, 2014; Gray et al., 2018).
  
  These shrinkage methods perform best when an external dataset is used to
  create a target matrix/target matrices that is informative of the actual
  dataset under examination. An example of this utility is provided in Gray et
  al. (2018), in which high-dimensional protein covariance matrices for
  various cancer types are greatly informed by large sample covariance matrices
  from 'similar' cancer types.
}
\author{
Harry Gray

Maintainer: Harry Gray <h.w.gray@dundee.ac.uk>
}
\references{
Gray, H., Leday, G.G.R., Vallejos, C.A. and Richardson, S.,
2018. Shrinkage estimation of large covariance matrices using 
multiple shrinkage targets. \href{https://arxiv.org/abs/1809.08024}{arXiv preprint}.

Hannart, A. and Naveau, P., 2014. Estimating high dimensional 
covariance matrices: A new look at the Gaussian conjugate framework. 
Journal of Multivariate Analysis, 131, pp.149-162. 
\href{http://dx.doi.org/10.1016/j.jmva.2014.06.001}{doi}.
}
\keyword{ package }

