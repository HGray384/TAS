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
#'
#' @examples
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
