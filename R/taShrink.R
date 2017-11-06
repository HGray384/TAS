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
  }
  # user-specified targets
  else
  {
    models <- targets
  }
  # remove unwanted targets
  if (length(without) > 1)
  {
    models <- models[,,-without]
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