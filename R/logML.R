#' @export
logML <- function(X, target, alpha)
{
  # useful variables
  X <- t(scale(t(X), center = T, scale=F))
  n <- ncol(X)
  p <- nrow(X)
  S <- tcrossprod(X)/n
  beta <- alpha/(1-alpha)
  nu <- beta*n + p + 1


  # first terms are a quotient of gamma functions
  r <- length(alpha)
  j <- 1:p
  y1 <- matrix(rep(nu, p), r, p)
  w <- t(matrix(rep(1-j, r), p, r))
  z1 <- lgamma((y1+w)/2)
  z1 <- rowSums(z1)

  y2 <- matrix(rep(nu+n, p), r, p)
  z2 <- lgamma((y2+w)/2)
  z2 <- rowSums(z2)

  gammaTerm <- z2-z1

  # terms of determinant quotient
  logdet1 <- nu/2*sapply(beta, function(x){determinant(x = x*target, logarithm=TRUE)$modulus})
  logdet2 <- (nu+n)/2*sapply(beta, function(x){determinant(x=x*target+S, logarithm=TRUE)$modulus})

  logdetTerm <- logdet1-logdet2

  # total
  gammaTerm+logdetTerm-((n*p)/2)*log(n*pi)
}
