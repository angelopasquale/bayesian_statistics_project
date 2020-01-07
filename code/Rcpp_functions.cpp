#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double tnormRcpp(double lo, double hi, double mu, double sig){

  double q1, q2, z;

  q1 = Rf_pnorm5(lo,mu,sig,1,0);
  q2 = Rf_pnorm5(hi,mu,sig,1,0);
  z = q1 + unif_rand()*(q2-q1);
  z = Rf_qnorm5(z, mu, sig, 1, 0);

  if(z > hi){
    z = lo;
  }

  if(z < lo){
    z = hi;
  }
  return(z);
}
