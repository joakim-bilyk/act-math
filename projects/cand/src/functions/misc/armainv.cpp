#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat armaInv(arma::mat x) {
  return arma::inv(x);
}