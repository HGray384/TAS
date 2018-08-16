#' Bayesian Target-Averaged linear Shrinkage (TAS) covariance estimator 
#' 
#' Implements a Bayesian target-averaged linear shrinkage covariance estimator
#' as in Gray et al (submitted) (pre-print available upon request).
#' It is most useful when the observed data is high-dimensional 
#' (more variables than observations) and there are other 
#' datasets that can be used to include as prior data-driven targets to shrink 
#' towards.
#'
#' @param X \code{matrix} -- data matrix with variables in rows and 
#' observations in columns. This method performs best when there are more 
#' variables than observations.
#' @param targets \code{character} or \code{array} -- "default" creates 
#' a target set of common literature targets, or the user may specify an array
#' of targets to use, e.g. ones that have been derived from external data. All
#' targets must be real symmetric positive definite matrices.
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
#' \item{alpha}{\code{list} -- the values of shrinkage intensities used.}
#' }}
#' @references Harry Gray, Gwenael G.R. Leday, Catalina A.
#' Vallejos, Sylvia Richardson (submitted). Target-averaged
#' linear Shrinkage: high dimensional covariance matrix
#' estimation in functional genomics.
#' @export
#'
#' @examples
#' set.seed(101)
#' X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
#' X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
#' tas <- taShrink(X, plots = FALSE) # apply shrinkage and view target weight bar plot
#' barplot(targetWeights(tas), names.arg = c(1:9, "S"),
#' main = "Target-specific shrinkage weights",
#' col = rainbow(dim(tas$targets)[3]+1), space = 0,
#' xlab = "Target", ylab = "Weight")
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
taShrink <- function(X, targets="default", without=0,
                     alpha = seq(0.01, 0.99, 0.01), 
                     plots = TRUE, ext.data=FALSE)
{
  ## input handling
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
            imputing these before running TAS")
  }
  if(any(is.infinite(X))){
    stop("X cannot contain infinite values!")
  }
  if(nrow(X)<ncol(X)){
    warning("TAS was designed for high-dimensional data analysis, but the 
            number of variables (p) in your X is less than the number 
            of samples (n). \n 
            If you know that this warning is incorrect then likely you need to 
            transpose X and run TAS again. Otherwise, know that
            TAS might not be the most suitable method here. ")
  }
  # the target input
  if(targets!="default" && !is.array(targets)){
    stop("The targets must be either 'default' or an array!")
  }
  if(is.array(targets)){
    if(dim(targets)[1]!=dim(targets)[2] || dim(targets)[1]!=nrow(X)){
      stop("Dimensions 1 and 2 of the 'targets' array must be equal to p.")
    }
    if(length(dim(targets))==2){
      message("only 1 target matrix specified... redirecting to gcShrink")
      return(gcShrink(X = X, target = targets, alpha = alpha, plots = plots,
                      weighted = TRUE, ext.data = ext.data))
    }
    if(any(!apply(targets, 3, is.numeric))){
      stop("The targets must be numeric matrices.")
    }
    if(any(is.na(targets)) || any(is.nan(targets))){
      stop("The supplied target matrices contain missing values.")
    }
    if(any(is.infinite(targets))){
      stop("The target matrices cannot contain infinite values!")
    }
    if(any(!apply(targets, 3, isSymmetric.matrix))){
      stop("At least one target matrix is not symmetric:
           tested using isSymmetric.matrix().")
    }
    if(any(apply(targets, 3, function(x){kappa(x, exact = TRUE)==0}))){
      stop("At least one target matrix is not positive definite.")
    }
    if(any(apply(targets, 3, function(x){kappa(x, exact = TRUE)>1e12}))){
      stop("At least one target matrix is not positive definite.")
    }
    if(any(apply(targets, 3, function(x){kappa(x, exact = TRUE)>1e4}))){
      warning("At least one target matrix is ill-conditioned, all results
              may contain numerical error.")
    }
  } else {
    # the without input
    if(!is.numeric(without)){
      stop("'without' must be numeric!")
    }
    if(any(is.na(without)) || any(is.nan(without))){
      stop("'without' contains missing values, please re-specify")
    }
    if(any(is.infinite(without))){
      stop("'without' cannot contain infinite values!")
    }
    without <- unique(without)
    without <- sort(intersect(without, 0:9))
    if (length(without)>1 && any(without==0)){
      without <- without[-which(without==0)]
    }
    if(length(without)==0){
      stop("'without' must be at least one of 0:9.")
    }
    if(length(without)==9){
      stop("'without' cannot be all of 1:9 - you will have no targets!")
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
  p <- nrow(X)
  
  # center the data
  X <- t(scale(t(X), scale=F, center=T))
  
  # create target set
  if (!is.array(targets))
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
    if (without!=0)
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
  logmargs <- t(apply(models, 3,
                      function(x){logML(X=X, target=x, alpha=alpha)}))
  
  # compute the posterior model weights
  weights <- exp(logmargs - matrixStats::logSumExp(logmargs))
  
  # compute the shrinkage intensity weights
  shrinkageweights <- weights%*%alpha
  
  # the weight allocated to the sample covariance
  sweight <- 1-sum(shrinkageweights)
  shrinkageweightslist <- simplify2array(lapply(shrinkageweights,
                                            function(x){matrix(x, p, p)}))
  
  # compute the estimate
  sigmahat <- apply(shrinkageweightslist*models, c(1,2), sum)
  sigmahat <- sigmahat + (sweight*tcrossprod(X)/n)
  
  # make target-weight plot
  if (plots){
    if (!is.array(targets)){
      if (length(which(without==0))==0){
        nm <- c(1:9)[-without]
      } else {
        nm <- 1:dim(models)[3]
      }
    } else {
      nm <- 1:dim(models)[3]
    }
    barplot(c(shrinkageweights, sweight), names.arg = c(nm, "S"),
            main = "Target-specific shrinkage weights",
            col = rainbow((dim(models)[3]+1)), space = 0, 
            xlab = "Target", ylab = "Weight")
  }
  # return results
  list("sigmahat" = sigmahat, "targets" = models, "weights" = weights, 
       "logmarginals" = logmargs, "alpha"=alpha)
}