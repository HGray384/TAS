#' Extract the target-specific shrinkage weights from TAS output
#'
#' @param TASoutput \code{list} -- output from the taShrink function.
#'
#' @return \code{list} -- the weights from each target and sample
#' covariance matrix in TAS.
#' @export
#'
#' @examples
targetWeights <- function(TASoutput){
  # compute the shrinkage intensity weights
  shrinkageweights <- TASoutput$weights%*%TASoutput$alpha
  
  # the weight allocated to the sample covariance
  sweight <- 1-sum(shrinkageweights)
  
  return(c(shrinkageweights, sweight))
}