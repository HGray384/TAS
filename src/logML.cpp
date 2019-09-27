#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//' Log-marginal likelihood of a Gaussian-inverse Wishart
//' conjugate model
//' 
//' Evaluate the log-marginal likelihood of a Gaussian-inverse Wishart
//' distribution parametrised in terms of its prior mean matrix and its
//' prior variance parameter. In the Bayesian linear shrinkage model,
//' these parameters correspond to the target matrix and the shrinkage
//' intensity (Hannart and Naveau, 2014). 
//' 
//' 
//' @param X \code{matrix} --data matrix with variables in rows and 
//' observations in columns.
//' @param target \code{matrix} -- prior mean matrix parameter of the
//' inverse-Wishart distribution. 
//' @param alpha \code{numeric} -- prior variance parameter of the
//' inverse-Wishart distribution. 
//' @return \code{numeric} -- log-marginal likelihood evaluated at
//'  (\code{target}, \code{alpha}). If \code{alpha} is a vector is a vector
//'  then the function returns a vector evaluated at each element of 
//'  \code{alpha}.
//' @seealso \code{\link{gcShrink}}, \code{\link{taShrink}}
//' @examples
//'   set.seed(102)
//'   X <- matrix(rnorm(50), 10, 5) # p=10, n=5, identity covariance
//'   X <- t(scale(t(X), center=TRUE, scale=FALSE)) # mean 0
//'   target <- getTarget(X)
//'   alpha <- seq(0.01, 0.99, 0.01)
//'   lml <- logML(X, target, alpha)
//'   plot(alpha, lml, col = 'blue', pch = 16,
//'   ylab = "log marginal likelihoods", xlab = expression(alpha))
//'   lines(x = rep(alpha[which(lml==max(lml))], 2), y = c(min(lml), max(lml)), col='red')
//' @references Alexis Hannart and Philippe Naveau (2014). 
//' Estimating high dimensional covariance matrices: 
//' A new look at the Gaussian conjugate framework. 
//' Journal of Multivariate Analysis. \href{http://dx.doi.org/10.1016/j.jmva.2014.06.001}{doi}.
// [[Rcpp::export]]
arma::vec logML(arma::mat X, arma::mat target, arma::vec alpha) {
  
  // useful quantities
  const float n = X.n_cols;
  const float p = X.n_rows;
  const float a = alpha.n_elem;
  const float constTerm = ((n * p) / 2) * log(n * arma::datum::pi);
  arma::mat S = X * X.t() / n;
  arma::vec beta = alpha / (1-alpha);
  arma::vec nu = beta * n + p + 1;
  arma::vec gammaTerm(a);
  arma::vec detTerm(a);
  double val1;
  double sign1;
  double val2;
  double sign2;
  
  // calculate the quotient of mvtgamma functions
  gammaTerm.zeros();
  detTerm.zeros();
  for (int r = 0; r < a; r++){
    for(int j = 0; j < p; j++){
      gammaTerm(r) += lgamma((nu(r) + n + (1- j))/2) - lgamma((nu(r) + (1- j))/2);
    }
    log_det(val1, sign1, beta(r) * target);
    detTerm(r) = nu(r) * val1;
    log_det(val2, sign2, beta(r) * target + S);
    detTerm(r) -= (nu(r) + n) * val2;
  }
  detTerm/=2;
  
  // return results
  //List ret;
  //ret["gammaTerm"] = gammaTerm;
  //ret["detTerm"] = detTerm;
  //ret["constTerm"] = constTerm;
  //ret["logmarg"] = gammaTerm + detTerm - constTerm;
  
  return gammaTerm + detTerm - constTerm;
}