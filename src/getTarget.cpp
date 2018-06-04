#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//' Construct a target matrix for single target linear shrinkage
//'
//' Construct a popular target matrix from the linear shrinkage literature.
//' These targets consist of a combination of variance and correlation
//'  structure. Possible variance structures are unit, sample mean, and 
//'  sample. Possible correlation structures are zero, sample mean, and
//'   autocorrelation.
//' 
//'
//' @param X \code{matrix} -- data matrix with variables in rows and 
//' observations in columns.
//' @param varNumber \code{numeric} -- c(1, 2, 3) variance structure for the
//'  target matrix.
//' 1 sets all variances equal to 1. 2 sets all variances equal to their sample
//'  mean (using \code{X}).
//' 3. sets all variances to their sample values (using \code{X}).
//' @param corNumber \code{numeric} -- c(1, 2, 3) correlation structure for
//'  the target matrix. 1 sets the correlations to 0. 2 sets the correlations
//'   equal 
//' to their sample mean (using \code{X}). 3 sets the correlations equals to an
//'  autocorrelation
//' structure with parameter equal to the sample mean (using \code{X}).
//' @return \code{matrix} -- target matrix for linear shrinkage
//'  estimation.
//'
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

