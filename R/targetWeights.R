#' Extract the target-specific and sample covariance shrinkage weights from TAS output
#'
#' @param TASoutput \code{list} -- output from the \code{taShrink} function.
#'
#' @return \code{list} -- the weights from each target and sample
#' covariance matrix in TAS.
#' @export
#' @seealso \code{\link{taShrink}}
#'
#' @examples
#'  set.seed(102)
#'  X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
#'  X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
#'  targets <- getTargetSet(X)[,,c(1, 4, 7)] # use unit variance targets
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
#'  xlab = "Target", ylab = "Weight", ylim = c(0, 0.6))
#'  barplot(tw2, names.arg = c("target1", "target2", "target3", "S"),
#'  main = "Target-specific shrinkage weights",
#'  col = c("red", "green", "blue", "purple"), space = 0, 
#'  xlab = "Target", ylab = "Weight", ylim = c(0, 0.6))
#'  par(mfrow=c(1, 1))
targetWeights <- function(TASoutput){
  # compute the shrinkage intensity weights
  shrinkageweights <- TASoutput$weights%*%TASoutput$alpha
  
  # the weight allocated to the sample covariance
  sweight <- 1-sum(shrinkageweights)
  
  return(c(shrinkageweights, sweight))
}