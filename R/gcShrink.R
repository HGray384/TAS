#' Bayesian Gaussian conjugate (GC) single target linear shrinkage covariance estimator 
#' 
#' Implements a Bayesian Gaussian conjugate (GC) single target linear shrinkage 
#' covariance estimator as in Gray et al (submitted) (pre-print available upon
#' request) and Hannart and Naveau (2014). 
#' It is most useful when the observed data is high-dimensional 
#' (more variables than observations) and allows a user-specified target matrix.
#'
#' @param X \code{matrix} -- data matrix with variables in rows and 
#' observations in columns. This method performs best when there are more 
#' variables than observations.
#' @param target \code{character} or \code{matrix} -- if "none" then a default
#' target specified by \code{var} and \code{cor} will be used for shrinkage. If
#' \code{matrix} then this will be used as the target for shrinkage. The
#' target must be a real symmetric positive definite matrix.
#' @param var \code{numeric} -- c(1, 2, 3) variance structure for the target matrix.
#' 1 sets all variances equal to 1. 2 sets all variances equal to their sample mean.
#' 3. sets all variances to their sample values.
#' @param cor \code{numeric} -- c(1, 2, 3) correlation structure for the 
#' target matrix. 1 sets the correlations to 0. 2 sets the correlations equal 
#' to their sample mean. 3 sets the correlations equals to an autocorrelation
#' structure with parameter equal to the sample mean.
#' @param alpha \code{list} -- the grid of shrinkage intensities in (0, 1) to be used. 
#' Recommended to be an equidistant grid that covers the whole interval. A short
#' comparison of estimation accuracy versus granularity is provided in ...
#' @param plots \code{logical} -- if TRUE then plots the log-marginal likelihood
#' for each value of alpha with the value of alpha that maximises this highlighted.
#' @param weighted \code{logical} -- if TRUE then average over all values of 
#' alpha and their respective marginal likelihood value as in Gray et al (submitted). 
#' If FALSE then only use the value of alpha that maximises the log-marginal likelihood
#' as in Hannart and Naveau (2014).
#' @param ext.data \code{matrix} -- an external data matrix used a surrogate to estimate 
#' the parameters in the default target set for X. Never recommended 
#' unless there is a belief that ext.data is informative of the covariances of
#'  X.
#'
#' @return {\code{list} --
#' \describe{
#' \item{sigmahat}{\code{matrix} -- the estimated covariance matrix.}
#' \item{optimalpha}{\code{numeric} -- the value of alpha that maximises the 
#' log-marginal likelihood.}
#' \item{target}{\code{matrix} -- the target matrix used for shrinkage.}
#' \item{logmarg}{\code{numeric} -- the values of the log marginal 
#'likelihood for each (target, alpha) pair. }
#' }}
#' @references Harry Gray, Gwenael G.R. Leday, Catalina A.
#' Vallejos, Sylvia Richardson (submitted). Target-averaged
#' linear Shrinkage: high dimensional covariance matrix
#' estimation in functional genomics.
#' 
#' Alexis Hannart and Philippe Naveau (2014). 
#' Estimating high dimensional covariance matrices: 
#' A new look at the Gaussian conjugate framework. 
#' Journal of Multivariate Analysis.
#' @export
#'
#' @examples
#' set.seed(102)
#' X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
#' X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
#' t1 <- gcShrink(X, var=1, cor=1) # apply shrinkage and view likelihood for T1
#' t2 <- gcShrink(X, var=2, cor=2) # apply shrinkage and view likelihood for T2
#' norm(t1$sigmahat-diag(10), type="F") # calculate loss
#' norm(t2$sigmahat-diag(10), type="F") # calculate loss
#' # one target clearly better but how to choose this a priori?
gcShrink <- function(X, target="none", var=2, cor=1, alpha = seq(0.01, 0.99, 0.01),
                     plots = TRUE, weighted=FALSE, ext.data=FALSE)
{
  # input handling
  if(!is.numeric(X)){
    message("The data matrix must be numeric!")
    stop()
  }
  if(target!="none" && !is.matrix(target)){
    message("The target must be either 'none' or a matrix!")
    stop()
  }
  if(!is.numeric(var)){
    message("'var' must be numeric!")
    stop()
  } else if (!any(var==c(1, 2, 3))){
    message("'var' must be in c(1, 2, 3)!")
  }
  if(!is.numeric(cor)){
    message("'cor' must be numeric!")
    stop()
  } else if (!any(var==c(1, 2, 3))){
    message("'cor' must be in c(1, 2, 3)!")
  }
  if(!is.numeric(alpha)){
    message("The shrinkage parameters 'alpha' must be numeric!")
    stop()
  }
  if(any(alpha<=0) || any(alpha>=1)){
    message("The shrinkage parameters must be within, and not inclusive of, (0, 1)!")
    stop()
  }
  if(!is.logical(plots)){
    message("'plots' must TRUE or FALSE!")
    stop()
  }
  if(!is.logical(weighted)){
    message("'weighted' must TRUE or FALSE!")
    stop()
  }
  if(is.logical(ext.data)){
    if(ext.data){
      message("Instead of entering ext.data=TRUE,
              set ext.data to be your external data matrix")
      stop()
    }
  }else if(!is.numeric(ext.data)){
    message("ext.data should either be your external 
            data matrix, or FALSE")
    stop()
  }
  # data dimensions
  n <- ncol(X)
  # p <- nrow(X)
  
  ## center the data
  X <- as.matrix(X)
  X <- t(scale(t(X), scale=F, center=T))
  
  
  # if target not specified
  if (is.character(target) && target=="none")
  {
    # get the target
    if(!ext.data){
      target <- getTarget(X, var, cor)
    }
    else
    {
      ext.data <- as.matrix(ext.data)
      target <- getTarget(ext.data, var, cor)
    }
  }
  
  # compute log marginals
  logmarg <- logML(X, target, alpha)
  
  # max of these gives optimal alpha
  optimalpha <- alpha[which(logmarg==max(logmarg))]
  if (plots){
    plot(alpha, logmarg, col = 'blue', pch = 16,
         ylab = "log marginal likelihoods", xlab = expression(alpha))
    
    lines(x = rep(optimalpha, 2), y = c(min(logmarg), max(logmarg)), col='red')
  }
  
  # compute the estimate
  if (weighted==FALSE){
    # put all weight on alpha with highest log marginal
    sigmahat <- ((1 - optimalpha) * tcrossprod(X)/n) + (optimalpha * target)
  } else {
    # weight each value of alpha
    marg <- exp(logmarg - matrixStats::logSumExp(logmarg))
    weights <- alpha*marg
    sigmahat <- ((1-sum(weights))*tcrossprod(X)/n) + (sum(weights)*target)
  }
  
  # # compute the variance for each individual element of sigmahat
  # # parameters
  # nu <- (optimalpha/(1-optimalpha)) + p + 1
  # Psi <- sigmahat * (optimalpha/(1-optimalpha))
  #
  # # compute the numerator of the expression for the variances
  # variances <- matrix(0, nrow(sigmahat), ncol(sigmahat))
  # for (i in 1:nrow(sigmahat)){
  #   for (j in 1:ncol(sigmahat)){
  #     variances[i, j] <- (nu - p + 1)*Psi[i, j]^2 + (nu - p - 1)*Psi[i, i]*Psi[j, j]
  #   }
  # }
  #
  # # Divide all of these by the denominator of the variances
  # denvar <- (nu - p)*((nu - p - 1)^2)*(nu - p - 3)
  # variances <- variances/denvar
  
  # return the estimated matrix, the optimal shrinkage intensity,
  # the sample covariance matrix and the log-marginal likelihood values
  list("sigmahat" = sigmahat, "optimalpha" = optimalpha,
       "target"= target,
       #"variances"=variances,
       "logmarg" = logmarg)
}