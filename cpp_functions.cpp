#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat l_qr_covariance(arma::mat x, arma::colvec density, double alpha) {
  try {
    int n = x.n_rows;
    int k = x.n_cols;
    
    // Define some 0-matrices
    arma::mat lambda = arma::zeros<arma::mat>(k, k);
    arma::mat C = arma::zeros<arma::mat>(k, k);
    
    arma::mat xi, xx;
    
    // Compute the matrix elements
    for (int i = 0; i < n; i++) {
      xi = x.row(i);
      xx = xi.t() * xi;
      lambda += xx * density(i);
      C += (1-alpha) / alpha * xx * pow(alpha, 2);
    }
    
    // Compute the covariance
    arma::mat cov = inv(lambda/n) * (C/n) * inv(lambda/n) / n;
    
    return cov;
  } catch(...) {
    Rcpp::stop("Cannot compute the covariance!");
  }
}
