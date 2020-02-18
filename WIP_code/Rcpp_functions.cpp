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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat trMVNmatrixRcpp(arma::mat avec, arma::mat muvec,
                                arma::mat smat, arma::mat lo,
                                arma::mat hi, arma::uvec whichSample,
                                arma::uvec idxALL){
  int cindex;
  arma::rowvec av;
  arma::rowvec mv;
  arma::vec mAs(2);
  int nm = smat.n_rows;
  int nr = muvec.n_rows;
  arma::rowvec p1(nm-1);
  arma::mat sin(nm-1, nm-1);
  arma::uvec cid(1);
  arma::uvec idx;
  arma::mat m1(1,1);
  arma::mat s1(1,1);
  double tiny = min(smat.diag())*.0001;
  int nk = whichSample.n_elem;

  arma::mat A(nr, nm); A.fill(NA_REAL);
  arma::umat idxALLm(nm-1, nm);

  for(int j=0; j < nm; j++)

    idxALLm.col(j) = idxALL.elem( find(idxALL != j) );

  for(int i = 0; i < nr ; i++){

    for(int k = 0; k < nk; k++){

      cindex = whichSample[k] - 1;

      av = avec.row(i);
      mv = muvec.row(i);

      cid(0) = cindex;
      idx = idxALLm.col(cindex);
      sin = arma::inv_sympd(smat.submat(idx, idx));
      p1 = trans(smat.submat(idx, cid)) * sin;

      m1 = mv[cindex] + dot(p1, (av.elem(idx) - mv.elem(idx)));
      s1 = smat(cindex,cindex) - dot(p1, smat.submat(cid, idx)) ;

      mAs[0] = m1(0,0);
      mAs[1] = s1(0,0);
      if(mAs[1] < 0) mAs[1] = tiny;

      double sss = pow(mAs[1],.5);

      avec(i,cindex) = tnormRcpp(lo(i,cindex), hi(i,cindex), mAs[0], sss);
      A(i,cindex) = avec(i,cindex);

    }
  }
  return A;
}
