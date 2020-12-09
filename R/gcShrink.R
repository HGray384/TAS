#' Bayesian Gaussian conjugate (GC) single target linear shrinkage covariance estimator 
#' 
#' Implements a Bayesian Gaussian conjugate (GC) single target linear shrinkage 
#' covariance estimator as in Gray et al. (2018)  and Hannart and Naveau (2014). 
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
#' @references Gray, H., Leday, G.G., Vallejos, C.A. and Richardson, S.,
#'  2018. Shrinkage estimation of large covariance matrices using 
#'  multiple shrinkage targets. \href{https://arxiv.org/abs/1809.08024}{arXiv preprint}.
#' 
#' Hannart, A. and Naveau, P., 2014. Estimating high dimensional 
#' covariance matrices: A new look at the Gaussian conjugate framework. 
#' Journal of Multivariate Analysis, 131, pp.149-162. 
#' \href{http://dx.doi.org/10.1016/j.jmva.2014.06.001}{doi}.
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
  # the data matrix X
  if(!is.numeric(X)){
    stop("X must be numeric!")
  }
  if(is.null(dim(X)) || length(dim(X))!=2){
    stop("X must have 2 dimensions.")
  }
  X <- as.matrix(X)
  if(length(X)==1){
    stop("X is a single number.")
  }
  if(any(is.na(X)) || any(is.nan(X))){
    stop("X contains missing values, consider 
         imputing these before performing shrinkage")
  }
  if(any(is.infinite(X))){
    stop("X cannot contain infinite values!")
  }
  if(nrow(X)<ncol(X)){
    warning("Shrinkage was designed for high-dimensional data analysis, but the 
            number of variables (p) in your X is less than the number 
            of samples (n). \n 
            If you know that this warning is incorrect then likely you need to 
            transpose X and run this again. Otherwise, know that
            shrinkage might not be the most suitable method here. ")
  }
  # target input
  if(is.character(target)){
    if (target!="none"){
      stop("If 'target' is a character, it must be 'none'")
    }
  }
  if(!is.character(target) && !is.matrix(target)){
    stop("The target must be either 'none' or a matrix!")
  }
  if(is.matrix(target)){
    if(dim(target)[1]!=dim(target)[2] || dim(target)[1]!=nrow(X)){
      stop("Dimensions 1 and 2 of the target matrix must be equal to p.")
    }
    if(!is.numeric(target)){
      stop("The target must be a numeric matrix.")
    }
    if(any(is.na(target)) || any(is.nan(target))){
      stop("The supplied target matrix contains missing values.")
    }
    if(any(is.infinite(target))){
      stop("The target matrix cannot contain infinite values!")
    }
    if(!isSymmetric.matrix(target)){
      stop("The target matrix is not symmetric:
           tested using isSymmetric.matrix().")
    }
    if(kappa(target, exact=TRUE)==0){
      stop("T target matrix is not positive definite.")
    }
    if(kappa(target, exact=TRUE)>1e12){
      stop("The target matrix is not positive definite.")
    }
    if(kappa(target, exact=TRUE)>1e4){
      warning("The target matrix is ill-conditioned, all results
              may contain numerical error.")
    }
  } else {
    if(length(var)!=1){
      stop("'var' must have length 1")
    }
    if(!is.numeric(var)){
      stop("'var' must be numeric!")
    }
    if (!any(var==c(1, 2, 3))){
      stop("'var' must be in c(1, 2, 3)!")
    }
    if(length(cor)!=1){
      stop("'cor' must have length 1")
    }
    if(!is.numeric(cor)){
      stop("'cor' must be numeric!")
    }
    if (!any(var==c(1, 2, 3))){
      stop("'cor' must be in c(1, 2, 3)!")
    }
  }
  # the alpha input
  if(!is.numeric(alpha)){
    stop("The shrinkage values 'alpha' must be numeric!")
  }
  if(any(is.na(alpha)) || any(is.nan(alpha))){
    stop("alpha contains missing values, please re-specify")
  }
  if(any(is.infinite(alpha))){
    stop("alpha cannot contain infinite values!")
  }
  if(any(alpha<=0) || any(alpha>=1)){
    stop("The shrinkage parameters must be within, and not inclusive of, (0, 1)!")
  }
  alpha <- sort(alpha)
  # the plot input
  if(is.na(plots) || is.nan(plots)){
    stop("'plots' is missing")
  }
  if(!is.logical(plots)){
    stop("'plots' must TRUE or FALSE!")
  }
  # the weighted input
  if(is.na(weighted) || is.nan(weighted)){
    stop("'weighted' is missing")
  }
  if(!is.logical(weighted)){
    stop("'weighted' must TRUE or FALSE!")
  }
  # the ext.data input
  if(is.logical(ext.data)){
    if(ext.data){
      stop("Instead of entering ext.data=TRUE,
             set ext.data to be your external data matrix")
    }
  } else {
    if(!is.numeric(ext.data)){
      stop("ext.data should either be your external 
               data matrix, or FALSE")
    }
    if(nrow(ext.data)!=nrow(X)){
      stop("ext.data should have the same number of rows as X.")
    }
    if(any(is.na(ext.data)) || any(is.nan(ext.data))){
      stop("The external data matrix contains missing values, consider 
               imputing these before running TAS")
    }
    if(any(is.infinite(ext.data))){
      stop("The external data matrix cannot contain infinite values!")
    }
  }
  # data dimensions
  n <- ncol(X)
  # p <- nrow(X)
  
  ## center the data
  X <- t(scale(t(X), scale=F, center=T))
  
  
  # if target not specified
  if (!is.matrix(target))
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
  logmarg <- logML(t(X), target)
  
  # max of these gives optimal alpha
  optimalpha <- logmarg$alphaOpt
  if (plots){
    plot(logmarg$gridAlpha[,2], logmarg$gridAlpha[,3], col = 'blue', pch = 16,
         ylab = "log marginal likelihoods", xlab = expression(alpha))
    
    lines(x = rep(optimalpha, 2), y = c(min(logmarg$gridAlpha[,3]), 
                                        max(logmarg$gridAlpha[,3])), col='red')
  }
  
  # compute the estimate
  if (weighted==FALSE){
    # put all weight on alpha with highest log marginal
    sigmahat <- ((1 - optimalpha) * tcrossprod(X)/n) + (optimalpha * target)
  } else {
    # weight each value of alpha
    marg <- exp(logmarg$gridAlpha[,3] - matrixStats::logSumExp(logmarg$gridAlpha[,3]))
    weights <- logmarg$gridAlpha[,2]*marg
    
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