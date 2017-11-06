#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube getTargetSet(arma::mat X) {
  //std::cout << "Before quantities!" << "\n";
  // useful quantities
  const float n = X.n_cols;
  const float p = X.n_rows;
  arma::vec vars = arma::diagvec(X * X.t()) / n;
  const float vbar = sum(vars) / p;
  const float rbar = (arma::accu(arma::cor(X.t(), 1))-p) / (p * (p + 1));
  int ntargs = 9;
  arma::cube targets(p, p, ntargs);
  //std::cout << "After quantities! Before targets!" << "\n";

  // create targets
  targets.slice(0) = arma::eye(p, p);
  targets.slice(1) = vbar * arma::eye(p, p);
  targets.slice(2) = arma::diagmat(vars);
  targets.slice(3) = rbar * arma::ones(p, p) - rbar * arma::eye(p, p) + arma::eye(p, p);
  targets.slice(4) = targets.slice(3) - arma::eye(p, p) + targets.slice(1);
  targets.slice(5) = targets.slice(3) - arma::eye(p, p) + targets.slice(2);
  arma::vec rbarpow(p);
  
  //std::cout << "After targets! Before rbarpow!" << "\n";
  
  for (int j = 0; j < p; j++){
    rbarpow(j) = pow(rbar, j);
  }
  arma::mat cormat = arma::toeplitz(rbarpow, rbarpow);
  targets.slice(6) = cormat;
  targets.slice(7) = cormat - arma::eye(p, p) + targets.slice(1);
  targets.slice(8) = cormat - arma::eye(p, p) + targets.slice(2);
  
  return targets;
}

