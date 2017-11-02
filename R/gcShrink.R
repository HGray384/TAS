gcShrink <- function(X, target="none", var="eq", cor="id", alpha = seq(0.1,  0.9, 0.1),
                     plots = TRUE, weighted=FALSE, ext.data=FALSE)
{
  ## initialise useful values
  n <- ncol(X)
  
  ## center the data.. crucial
  X <- as.matrix(X)
  X <- t(scale(t(X), scale=F, center=T))
  S <- tcrossprod(X)/n
  
  
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
    sigmahat <- ((1 - optimalpha) * S) + (optimalpha * target)
  } else {
    # weight each value of alpha
    marg <- exp(logmarg - matrixStats::logSumExp(logmarg))
    weights <- alpha*marg
    sigmahat <- ((1-sum(weights))*S) + (sum(weights)*target)
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
  
  # return the estimated matrix, the optimal shrinkage intrensity,
  # the sample covariance matrix and the matrix of element-wise variances
  list("sigmahat" = sigmahat, "optimalpha" = optimalpha, "samplecov" = S,
       "target"= target,
       #"variances"=variances,
       "logmarg" = logmarg)
}