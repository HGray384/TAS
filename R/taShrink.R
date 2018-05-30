#' Bayesian Target-Averaged linear Shrinkage (TAS) covariance estimator 
#' 
#' Implements a Bayesian target-averaged linear shrinkage covariance estimator
#' as in (Gray et al 2018). It is most useful when the observed data is 
#' high-dimensional (more variables than observations) and there are other 
#' datasets that can be used to include as prior data-driven targets to shrink 
#' towards.
#'
#' @param X \code{matrix} -- data matrix with variables in rows and 
#' observations in columns. This method performs best when there are more 
#' variables than observations.
#' @param targets \code{character} or \code{array} -- "default" creates 
#' a target set of common literature targets, or the user may specify an array
#' of targets to use, e.g. ones that have been derived from external data.
#' @param without \code{list} -- if targets=="default" then this indicates which of the 
#' default targets should be excluded from shrinkage. This can be useful
#' when exploring the shrinkage behaviour with a subset of targets (e.g. 
#' through simulation).
#' @param alpha \code{list} -- the grid of shrinkage intensities in (0, 1) to be used. 
#' Recommended to be an equidistant grid that covers the whole interval. A short
#' comparison of estimation accuracy versus granularity is provided in ...
#' @param plots \code{logical} -- if TRUE then create a barplot of the target-specific 
#' shrinkage weights. Recommend option FALSE if using many iterations.
#' @param ext.data \code{matrix} -- an external data matrix used a surrogate to estimate 
#' the parameters in the default target set for X. Never recommended 
#' unless there is a belief that ext.data is informative of the covariances of
#'  X. 
#'
#' @return {\code{list} --
#' \describe{
#' \item{sigmahat}{\code{matrix} -- the estimated covariance matrix.}
#' \item{targets}{\code{array} -- the targets used for shrinkage.}
#' \item{weights}{\code{matrix} -- the weight of each (target, alpha)
#' pair such that \code{sum(weights)=1}. The weights are calculated by 
#' normalising the log-marginal likelihood values below. }
#' \item{logmarginals}{\code{matrix} -- the values of the log marginal 
#'likelihood for each (target, alpha) pair. }
#' }}
#' @export
#'
#' @examples
#' set.seed(101)
#' X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
#' X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
#' tas <- taShrink(X) # apply shrinkage and view target weight bar plot
#' abs(tas$sigmahat - diag(10)) # inspect absolute differences
#' norm(tas$sigmahat-diag(10), type="F") # calculate loss
#' # compare this to each single target
#' norm(gcShrink(X, var=1, cor=1)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=2, cor=1)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=3, cor=1)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=1, cor=2)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=2, cor=2)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=3, cor=2)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=1, cor=3)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=2, cor=3)$sigmahat-diag(10), type="F")
#' norm(gcShrink(X, var=3, cor=3)$sigmahat-diag(10), type="F")
#' 
#' 
#' 
taShrink <- function(X, targets="default",  without=0,
                     alpha = seq(0.01, 0.99, 0.01), plots = TRUE, ext.data=FALSE)
{
  # data dimensions
  n <- ncol(X)
  p <- nrow(X)
  
  # center the data
  X <- as.matrix(X)
  X <- t(scale(t(X), scale=F, center=T))
  
  # create target set
  if (targets == "default")
  {
    models <- array(0, dim = c(p, p, 9))
    if(!ext.data){
      # default target set derived from X
      models <- getTargetSet(X)
    }else{
      # default target set derived from external data
      ext.data <- as.matrix(ext.data)
      models <- getTargetSet(ext.data)
    }
    # remove unwanted targets
    if (length(without) > 1)
    {
      models <- models[,,-without]
    }
  }
  # user-specified targets
  else
  {
    models <- targets
  }
  
  
  # compute the log-marginal likelihood
  logmargs <- t(apply(models, 3, function(x){logML(X=X, target=x, alpha=alpha)}))
  # compute the posterior model weights
  weights <- exp(logmargs - matrixStats::logSumExp(logmargs))
  # compute the shrinkage intensity weights
  shrinkageweights <- weights%*%alpha
  # the weight allocated to the sample covariance
  sweight <- 1-sum(shrinkageweights)
  shrinkageweights <- simplify2array(lapply(shrinkageweights, function(x){matrix(x, p, p)}))
  # compute the estimate
  sigmahat <- apply(shrinkageweights*models, c(1,2), sum)
  sigmahat <- sigmahat + (sweight*tcrossprod(X)/n)
  
  # make target-weight plot
  if (plots){
    barplot(rowSums(weights)/sum(weights), names.arg = 1:dim(models)[3], main = "Distribution of target weights",
            col = rainbow(dim(models)[3]), space = 0, xlab = "Target", ylab = "Weight")
  }
  # return results
  list("sigmahat" = sigmahat, "targets" = models, "weights" = weights, "logmarginals" = logmargs)
}