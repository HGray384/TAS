# This file implements a function that generates a target matrix

getTarget <- function(X, var="eq", cor="id"){
  # Get some other useful parameters from the input
  
  n <- ncol(X)
  p <- nrow(X)
  
  # SETTING THE
  #             STANDARD DEVIATION MATRIX
  X <- as.matrix(X)
  S <- tcrossprod(X)/n
  if (var=="id"){
    SDs <- diag(p)
  } else if (var=="eq"){
    SDs <- diag(sqrt(sum(diag(S))/p), p, p)
  } else if(var=="un"){
    SDs <- diag(sqrt(diag(S)), p, p)
  }
  
  # SETTING THE
  #             CORRELATION MATRIX
  # The diagonal elements of the correlation matrix are always 1
  cormat <- diag(p)
  if (cor=="cc"){
    rho <- mean(cov2cor(S)[upper.tri(S)])
    rhomat <- matrix(rho, p, p)
    cormat[which(upper.tri(cormat))] <- rhomat[which(upper.tri(rhomat))]
    cormat[which(lower.tri(cormat))] <- rhomat[which(lower.tri(rhomat))]
  } else if (cor=="arc"){
    rho <- mean(cov2cor(S)[upper.tri(S)])
    powers <- 1:p
    H <- abs(outer(powers, powers, "-"))
    cormat <- rho^H
  }
  
  # COMPUTE THE
  #             TARGET MATRIX
  # The following formula computes the target covariance matrix
  # from the variance matrix and the correlation matrix
  target <- SDs %*% cormat %*% SDs
  return(as.matrix(target))
}