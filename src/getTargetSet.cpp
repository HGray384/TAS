#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//' Construct a set of target matrices for Target-Averaged linear shrinkage
//'
//' Construct a set of popular target matrices from the linear shrinkage
//'  literature.
//' These nine targets consist of the combinations of variance and correlation
//'  structures; variance structures are unit, sample mean, and 
//'  sample; correlation structures are zero, sample mean, and
//' autocorrelation.
//' 
//' @param X \code{matrix} -- data matrix with variables in rows and 
//' observations in columns.
//' @return \code{array} -- a \code{p}x\code{p}x9 array of
//' target matrices, where \code{p} is the number of variables of \code{X}.
// [[Rcpp::export]]
arma::cube getTargetSet(arma::mat X) {
  //std::cout << "Before quantities!" << "\n";
  // useful quantities
  const float n = X.n_cols;
  const float p = X.n_rows;
  arma::vec sds1(p);
  sds1.ones(p);
  arma::vec vars = arma::diagvec(X * X.t()) / n;
  arma::vec sds2(p);
  sds2.fill(sqrt(sum(vars) / p));
  arma::vec sds3 = sqrt(vars);
  
  arma::mat cormat1 = arma::eye(p, p);
  const float rbar = (arma::accu(arma::cor(X.t(), 1))-p) / (p * (p + 1));
  arma::mat cormat2 = rbar * arma::ones(p, p) - rbar * arma::eye(p, p) + arma::eye(p, p);
  arma::vec rbarpow(p);
  for (int j = 0; j < p; j++){
    rbarpow(j) = pow(rbar, j);
  }
  arma::mat cormat3 = arma::toeplitz(rbarpow, rbarpow);

  int ntargs = 9;
  arma::cube targets(p, p, ntargs);
  //std::cout << "After quantities! Before targets!" << "\n";

  // create targets
  targets.slice(0) = arma::eye(p, p);
  targets.slice(1) = arma::diagmat(sds2) * cormat1 * arma::diagmat(sds2);
  targets.slice(2) = arma::diagmat(sds3) * cormat1 * arma::diagmat(sds3);
  targets.slice(3) = arma::diagmat(sds1) * cormat2 * arma::diagmat(sds1);
  targets.slice(4) = arma::diagmat(sds2) * cormat2 * arma::diagmat(sds2);
  targets.slice(5) = arma::diagmat(sds3) * cormat2 * arma::diagmat(sds3);
  targets.slice(6) = arma::diagmat(sds1) * cormat3 * arma::diagmat(sds1);
  targets.slice(7) = arma::diagmat(sds2) * cormat3 * arma::diagmat(sds2);
  targets.slice(8) = arma::diagmat(sds3) * cormat3 * arma::diagmat(sds3);
  
  return targets;
}

