#' Add a new target to TAS without re-running taShrink()
#'
#' @param X \code{matrix} -- data matrix with variables in rows and 
#' observations in columns. This method performs best when there are more 
#' variables than observations.
#' @param TASoutput \code{list} -- output from the taShrink function.
#' @param NEWtarget \code{matrix} -- a new target to add to the shrinkage
#' estimation. Must have the same dimensions as the other targets.
#'
#' @return \code{list} -- the updated TAS output having added the new target matrix
#'  to the target set.
#' @export
#' @seealso \code{\link{taShrink}}
#'
#' @examples
#'  set.seed(102)
#'  X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
#'  X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
#'  targets <- getTargetSet(X)[,,c(1, 4, 7)] # use unit variance targets
#'  alpha <- seq(0.01, 0.99, 0.01)
#'  tas <- taShrink(X, targets = targets[,,c(1, 3)], plots = FALSE)
#'  tw1 <- targetWeights(tas)
#'  barplot(tw1, names.arg = c("target1", "target2", "S"),
#'  main = "Target-specific shrinkage weights",
#'  col = c("red", "green", "purple"), space = 0, 
#'  xlab = "Target", ylab = "Weight")
#'  tas2 <- addTarget(X, tas, targets[,,2])
#'  tw2 <- targetWeights(tas2)
#'  par(mfrow=c(1, 2))
#'  barplot(tw1, names.arg = c("target1", "target2", "S"),
#'  main = "Target-specific shrinkage weights",
#'  col = c("red", "green", "purple"), space = 0, 
#'  xlab = "Target", ylab = "Weight")
#'  barplot(tw2, names.arg = c("target1", "target2", "target3", "S"),
#'  main = "Target-specific shrinkage weights",
#'  col = c("red", "green", "blue", "purple"), space = 0, 
#'  xlab = "Target", ylab = "Weight")
#'  par(mfrow=c(1, 1))
#'  plot(alpha, tas2$logmarginals[1,], col = 'red', pch = 16,
#'  ylab = "log marginal likelihoods", xlab = expression(alpha))
#'  points(alpha, tas2$logmarginals[2,], col = 'green', pch = 16)
#'  points(alpha, tas2$logmarginals[3,], col = 'blue', pch = 16)
#'  legend('bottomright', c("target1", "target2", "target3"), pch = 16,
#'    col=c('red', 'green', 'blue'))
addTarget <- function(X, TASoutput, NEWtarget){
  # get new target log marginal weights
  addMarg <- gcShrink(X = X, target = NEWtarget,
                      alpha = TASoutput$alpha, plots = FALSE,
                      weighted = TRUE, ext.data = FALSE)$logmarg
  oldMarg <- TASoutput$logmarginals
  # add them to the old ones
  newMarg <- rbind(oldMarg, t(addMarg))
  
  # calculate new posterior model weights
  weights <- exp(newMarg - matrixStats::logSumExp(newMarg))
  
  # compute the shrinkage intensity weights
  shrinkageweights <- weights%*%TASoutput$alpha
  
  # the weight allocated to the sample covariance
  sweight <- 1-sum(shrinkageweights)
  shrinkageweightslist <- simplify2array(
    lapply(shrinkageweights,function(x){matrix(x, nrow(X), nrow(X))})
    )
  
  # compute the estimate
  models <- array(c(TASoutput$targets, NEWtarget), 
                  dim=c(nrow(X), nrow(X), dim(TASoutput$targets)[3]+1))
  sigmahat <- apply(shrinkageweightslist*models, c(1,2), sum)
  sigmahat <- sigmahat + (sweight*tcrossprod(X)/ncol(X))
  
  # return results
  list("sigmahat" = sigmahat, "targets" = models, "weights" = weights, 
       "logmarginals" = newMarg, "alpha"=TASoutput$alpha)
}
