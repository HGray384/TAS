taShrink <- function(X, targets="default",  without=0,
                     alpha = seq(0.1, 0.9, 0.1), plots = TRUE, ext.data=FALSE)
{
  # enter targets as an array with the third dimension
  # indicating the target number
  
  ## initialise useful values
  n <- ncol(X)
  p <- nrow(X)
  
  ## center the data.. crucial
  X <- as.matrix(X)
  X <- t(scale(t(X), scale=F, center=T))
  S <- tcrossprod(X)/n
  
  
  if (is.character(targets) && targets == "default")
  {
    models <- array(0, dim = c(p, p, 9))
    if(!ext.data){
      models[,,1] <- diag(p)
      models[,,2] <- getTarget(X, var="eq", cor="id")
      models[,,3] <- getTarget(X, var="un", cor="id")
      models[,,4] <- getTarget(X, var="id", cor="cc")
      models[,,5] <- getTarget(X, var="eq", cor="cc")
      models[,,6] <- getTarget(X, var="un", cor="cc")
      models[,,7] <- getTarget(X, var="id", cor="arc")
      models[,,8] <- getTarget(X, var="eq", cor="arc")
      models[,,9] <- getTarget(X, var="un", cor="arc")
    }else{
      ext.data <- as.matrix(ext.data)
      models[,,1] <- diag(p)
      models[,,2] <- getTarget(ext.data, var="eq", cor="id")
      models[,,3] <- getTarget(ext.data, var="un", cor="id")
      models[,,4] <- getTarget(ext.data, var="id", cor="cc")
      models[,,5] <- getTarget(ext.data, var="eq", cor="cc")
      models[,,6] <- getTarget(ext.data, var="un", cor="cc")
      models[,,7] <- getTarget(ext.data, var="id", cor="arc")
      models[,,8] <- getTarget(ext.data, var="eq", cor="arc")
      models[,,9] <- getTarget(ext.data, var="un", cor="arc")
    }
  }
  else
  {
    models <- targets
  }
  
  if (length(without) > 1)
  {
    models <- models[,,-without]
  }
  
  # ntargets <- dim(models)[3]
  
  ## compute the estimate
  logmargs <- t(apply(models, 3, function(x){logML(X, x, alpha)}))
  
  weights <- exp(logmargs - matrixStats::logSumExp(logmargs))
  
  shrinkageweights <- weights%*%alpha
  sweight <- 1-sum(shrinkageweights)
  shrinkageweights <- simplify2array(lapply(shrinkageweights, function(x){matrix(x, p, p)}))
  sigmahat <- shrinkageweights*models
  
  sigmahat <- apply(sigmahat, c(1,2), sum)
  sigmahat <- sigmahat + (sweight*S)
  
  
  
  #   ## output handling
  #   # control flow variables
  #   vars <- 1
  #   cors <- 1
  #   structs <- rep(0, ntargets)
  #   variances <- c("unit", "mean", "sample")
  #   correlations <- c("zero", "mean", "auto-regressive")
  #
  #   for (m in 1:ntargets)
  #   {
  #     if (targets == "default")
  #     {
  #       structs[m] <- paste("Target ", m, " Variances: ", variances[vars], " Correlations: ", correlations[cors])
  #     }
  #     else
  #     {
  #       structs <- paste("Target ", m)
  #     }
  #     print(structs[m])
  #     print(paste("Alpha: ", alpha, " Normalised marginal: ", weights[m,]))
  #     if (vars == 3)
  #     {
  #       vars <- 0
  #       cors <- cors + 1
  #     }
  #     vars <- vars + 1
  #   }
  if (plots){
    barplot(rowSums(weights)/sum(weights), names.arg = 1:dim(models)[3], main = "Distribution of target weights",
            col = rainbow(dim(models)[3]), space = 0, xlab = "Target", ylab = "Weight")
  }
  
  list("sigmahat" = sigmahat, "targets" = models, "weights" = weights, "logmarginals" = logmargs)
}