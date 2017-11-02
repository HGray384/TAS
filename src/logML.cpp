#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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
  
  // calculate the quotient of mvtgamma functions
  gammaTerm.zeros();
  detTerm.zeros();
  for (int r = 0; r < a; r++){
    for(int j = 0; j < p; j++){
      gammaTerm(r) += lgamma((nu(r) + n + (1- j))/2) - lgamma((nu(r) + (1- j))/2);
    }
    detTerm(r) = nu(r) * log(arma::det(beta(r) * target));
    detTerm(r) -= (nu(r) + n) * log(arma::det(beta(r) * target + S));
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