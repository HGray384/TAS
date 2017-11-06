#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat getTarget(arma::mat X, int varNumber = 2, int corNumber = 1) {
  //std::cout << "Before quantities!" << "\n";
  // useful quantities
  const float n = X.n_cols;
  const float p = X.n_rows;
  arma::mat target(p, p);
  arma::vec sds(p);
  arma::mat cormat(p, p);
  
  if(varNumber==1){
    sds.ones();
  }
  else {
    arma::vec vars = arma::diagvec(X * X.t()) / n;
    if (varNumber==2){
    sds.fill(sqrt(sum(vars) / p));
    } else {
      sds = sqrt(vars);
    }
  }
  
  if (corNumber == 1){
    cormat = arma::eye(p, p);
  } else {
    const float rbar = (arma::accu(arma::cor(X.t(), 1))-p) / (p * (p + 1));
    if (corNumber == 2){
      cormat = rbar * arma::ones(p, p) - rbar * arma::eye(p, p) + arma::eye(p, p);
    } else {
      arma::vec rbarpow(p);
      for (int j = 0; j < p; j++){
        rbarpow(j) = pow(rbar, j);
      }
      cormat = arma::toeplitz(rbarpow, rbarpow);
    }
  }
  
  target = arma::diagmat(sds) * cormat * arma::diagmat(sds);
  
  return target;
}

