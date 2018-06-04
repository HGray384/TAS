\name{TAS-package}
\alias{TAS-package}
\alias{TAS}
\docType{package}
\title{
  Target-Averaged Linear Shrinkage estimation
}
\description{
  Conjugate Bayesian covariance matrix estimation using linear shrinkage with multiple target matrices (Gray et al, submitted). Most useful in high-dimensional data settings, where the number of variables is greater than the number of samples.
}
\details{
  This package contains functions for covariance estimation using a conjugate
  Bayesian model. Whilst the main functionality of the package is for multiple
  target linear shrinkage estimation, we also provide functionality for the 
  single target analogue (Hannart and Naveau, 2014; Gray et al, submitted).
  
  These shrinkage methods perform best when an external dataset is used to
  create a target matrix/target matrices that is informative of the actual
  dataset under examination. An example of this utility is provided in Gray et
  al (submitted), in which high-dimensional protein covariance matrices for
  various cancer types are greatly informed by large sample covariance matrices
  from 'similar' cancer types.
}
\author{
Harry Gray

Maintainer: Harry Gray <harry.gray@mrc-bsu.cam.ac.uk>
}
\references{
  Harry Gray, Gwenael G.R. Leday, Catalina A.
Vallejos, Sylvia Richardson (submitted; available upon request).
Target-averaged linear Shrinkage: high dimensional covariance matrix
estimation in functional genomics. 

Alexis Hannart and Philippe Naveau (2014). 
Estimating high dimensional covariance matrices: 
A new look at the Gaussian conjugate framework. 
Journal of Multivariate Analysis.
}
\keyword{ package }

