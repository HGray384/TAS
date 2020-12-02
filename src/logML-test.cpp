//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]
#include <tuple>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/beta.hpp>
#undef NDEBUG
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using boost::math::tools::brent_find_minima;
using boost::math::ibetac;
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
//'   lml <- logMLtest(X, target, alpha)
//'   plot(alpha, lml, col = 'blue', pch = 16,
//'   ylab = "log marginal likelihoods", xlab = expression(alpha))
//'   lines(x = rep(alpha[which(lml==max(lml))], 2), y = c(min(lml), max(lml)), col='red')
//' @references Alexis Hannart and Philippe Naveau (2014). 
//' Estimating high dimensional covariance matrices: 
//' A new look at the Gaussian conjugate framework. 
//' Journal of Multivariate Analysis. \href{http://dx.doi.org/10.1016/j.jmva.2014.06.001}{doi}.
double alphaToDelta(double alpha, int n, int p){
  return (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha);
}

double deltaToAlpha(double delta, int n, int p){
  return (delta-p-1)/(n+delta-p-1);
}

double lpvarGamma(const double x, const int p) {
  double ans = (p * (p - 1) * 0.25) * log(datum::pi);
  for(int j = 1; j < (p + 1); ++j){
    ans += std::lgamma(x - ((j - 1.0) * 0.5));
  }
  return ans;
}

double logMLt(const double delta, const int p, const int n, colvec eigs, const double logdettarget){
  double out = -0.5*n*p*std::log(datum::pi);
  out += lpvarGamma((delta+n)*0.5, p);
  out -= lpvarGamma(delta*0.5, p);
  out += 0.5*delta*p*std::log(delta-p-1);
  out -= 0.5*(delta+n)*sum(arma::log((delta-p-1)+eigs));
  if(logdettarget!=0){
    out -= 0.5*n*logdettarget;
  }
  return(out);
}

double getDeltaOpt(const int n, const int p, colvec eigs, const double logdettarget){
  const double lowerVal = alphaToDelta(0.001, n, p);
  const double upperVal = alphaToDelta(0.999, n, p);
  const auto obj = [p, n, eigs, logdettarget](double x) { return -logMLt(x, p, n, eigs, logdettarget); };
  boost::uintmax_t it = 1000;
  const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  //std::pair<double, double> result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  auto deltaOpt = 0.0, valOpt = 0.0;
  std::tie(deltaOpt, valOpt) = result;
  return(deltaOpt);
}


// [[Rcpp::export("logMLtest")]]
List logMLtest(arma::mat X, arma::mat target) {
  
  // Dimension data
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  
  // Sample variances
  rowvec s = sum(square(X), 0)/n;
  
  
  // Scatter matrix
  arma::mat XTX;
  XTX = X.t()*X;
  
  arma::mat targetinv, choltarget;
  double logdettarget = 0.0;
  bool isD = !all(target.diag()==0);
  if(isD){
    choltarget = chol(target);
    if(all(choltarget.diag()>0)){
      logdettarget = 2*sum(log(choltarget.diag()));
      choltarget = inv(choltarget);
      targetinv = choltarget*choltarget.t();
    }else{
      choltarget.clear();
      isD = false;
      Rcpp::Rcout << "Warning: target is not positive definite and has been set to the identity matrix" << std::endl;
    }
  }
  
  // Eigenvalue decomposition
  arma::colvec eigs = zeros(p);
  double deltaOpt;
  
  if(n >= p){
    
    // Compute eigenvalues
    if(XTX.is_empty()){
      XTX = X.t()*X;
    }
    
    if(isD){
      
      arma::mat targetinvXTX = targetinv*XTX;
      cx_vec eigsval = eig_gen(targetinvXTX);
      targetinvXTX.clear();
      eigs += sort(real(eigsval), "descend");
      eigsval.clear();
      
    }else{
      
      eigs += sort(eig_sym(XTX), "descend");
      
    }
    
    // Optimal shrinkage
    deltaOpt = getDeltaOpt(n, p, eigs, logdettarget);
    
    
  }else{
    
    if(isD){
      
      // Compute eigenvalues
      arma::mat XtargetinvXT = X*targetinv*X.t();
      eigs.subvec(0, n-1) += sort(eig_sym(XtargetinvXT), "descend");
      
      // Optimal shrinkage
      deltaOpt = getDeltaOpt(n, p, eigs, logdettarget);
      XtargetinvXT.clear();
      
    }else{
      
      // Compute eigenvalues
      arma::mat XXT = X*X.t();
      eigs.subvec(0, n-1) += sort(eig_sym(XXT), "descend");
      
      // Optimal shrinkage
      deltaOpt = getDeltaOpt(n, p, eigs, logdettarget);
      XXT.clear();
      
    }
    
  }
  
  X.clear();
  const double valOpt = logMLt(deltaOpt, p, n, eigs, logdettarget);
  const double alphaOpt = deltaToAlpha(deltaOpt, n, p);
  
  // Values of log-ML for a grid of alphas
  arma::mat gridAlpha(100, 3, fill::zeros);
  gridAlpha.col(1) += linspace(0.01, 0.99, 100);
  for(int j=0; j<100; j++){
    gridAlpha(j, 0) = alphaToDelta(gridAlpha(j, 1), n, p);
    gridAlpha(j, 2) = logMLt(gridAlpha(j, 0), p, n, eigs, logdettarget);
  }
  eigs.clear();
  
  // Output object
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("deltaOpt") = deltaOpt,
                                      Rcpp::Named("alphaOpt") = alphaOpt,
                                      Rcpp::Named("gridAlpha") = gridAlpha,
                                      Rcpp::Named("valOpt") = valOpt,
                                      Rcpp::Named("s") = s);
  
  return out;
}