//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace Rcpp;

using std::cout;
using std::endl;
using std::vector;
using std::setw;

typedef unsigned int uint;
typedef long double real;

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::LLT<MatrixXd> LLT;

VectorXd rmvnorm2 ( const VectorXd & mu, const LLT &SigmaInvLLT ) {
  VectorXd z  = as<VectorXd> ( rnorm ( mu.rows(), 0, 1 ) );
  return mu + SigmaInvLLT.matrixL().solve ( z );
}

// GJAM model using algorithm from Taylor-Rodríguez, 2017
//[[Rcpp::export]]

List GJAM_Gibbs_sampler ( MatrixXd v, uint ndraws, uint burnin, uint trim = 1,
                   int r, double alpha0 = 1, double nu0 = 1, int N) {

}

// v: a n-by-S matrix with, on the i-th row, the vector of the observations for
// all species at site i

// ndraws, burnin, trim: the iteration parameters; the Gibbs sampler will run for
// a total of burnin + ndraws iterations, the first burnin iterations will be
// discarded and the remaining ndraws will be considered for the output, subsampling
// one iteration every trim to remove autocorrelation

// r: the number of latent factors which gives the size of the S × r matrix A
// A = QZ where Z is a N × r matrix representing all the vector values that
// the lines of A may take. Z is sampled at step 1.

// alpha0: the prior hyperparameters for the beta distribution of the draws of p

// nu0: the prior hyperparameters for the distribution of sigma_eps^2

// N: truncation level for stick-breaking
