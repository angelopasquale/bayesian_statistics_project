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

// Fits the ANOVA-DDP model using algorithm 2 from Neal, 2000
//[[Rcpp::export]]
List fitAnovaDDP ( MatrixXd y, MatrixXd h, uint ndraws, uint burnin, uint trim = 1,
                   double alpha0 = 1, double a0 = 1, double b0 = 1, double nu0 = 1, double kappa0 = 1 ) {
  RNGScope scope;

  // Import functions from R environment
  Environment invgamma = Environment::namespace_env("invgamma");
  Function rinvgamma = invgamma["rinvgamma"];

  // Identify NAs
  // We store the (cell, time) coordinates of the missing values for when we resample them in the Gibbs sampler loop
  vector<uint> naCells, naTimes; // NA cells and times
  for ( uint i = 0; i < y.rows(); ++i ) {
    for ( uint j = 0; j < y.cols(); ++j ) {
      if ( NumericMatrix::is_na(y(i,j)) ) {
        naCells.push_back ( i );
        naTimes.push_back ( j );
        y(i,j) = 0.0;
      }
    }
  }

  // cout.precision(4); cout << std::fixed;
  cout << "Found " << naTimes.size() << " NA values\n";

  // Dimensions of the model
  const uint cells = y.rows(); // Number of cells
  const uint times = y.cols(); // Number of times
  const uint p = h.cols(); // Number of covariates

  // State of the Markov chain; for each variable of interest we keep both the current state and a record of the states
  // visited by the chain

  // Cluster labels for each cell
  vector<uint> labels ( cells, 0 );
  vector<vector<uint>> labelsRecord;

  // Regression coefficients beta for each cluster (i.e. unique regression coefficients beta)
  MatrixXd beta = MatrixXd::Zero ( 1, p );
  vector<MatrixXd> betaRecord;

  // Numerosity of the clusters
  vector<uint> clusterSize ( 1, cells );

  // Similarity matrix; similarityMatrix(i,j) == 1 iff labels(i) == labels(j)
  // Given the size of the dataset, it costs too much memory to save the matrix at each iteration. We return only the
  // posterior similarity matrix, i.e. the ergodic mean of the similarity matrices
  MatrixXd similarityMatrix = MatrixXd::Constant ( cells, cells, 1 );
  MatrixXd psm = MatrixXd::Zero ( cells, cells );

  // Variance of the regression residuals
  real sigma2 = 10;
  vector<real> sigmaRecord;

  // Precision of the regression coefficients' DP base measure
  real tauBeta = 1.0;
  vector<real> tauBetaRecord;

  // Here we pre-allocate and pre-compute some of the matrices and vectors we will need during the execution of the
  // algorithm, to save CPU time and memory

  VectorXd tmpBetaEigen ( p );
  NumericVector tmpBetaRcpp ( p );
  MatrixXd tmpSigmaInv ( p, p );
  LLT tmpSigmaInvLLT = tmpSigmaInv.llt();
  MatrixXd hth = h.transpose() * h;
  real sigma1 = sqrt(sigma2);
  MatrixXd sigmaInv = hth / sigma2; sigmaInv.diagonal() += VectorXd::Constant(p,1,tauBeta);
  LLT sigmaInvLLT = sigmaInv.llt();
  real detSigma = 1 / sigmaInv.determinant();
  MatrixXd hBetaT = h * beta.transpose();

  // MCMC iterations
  for ( uint iter = 0; iter < burnin + ndraws; ++iter ) {
    cout << endl;

    // Step 1: for each cell, resample the cluster label
    for ( uint i = 0; i < cells; ++i ) { // Loop on cells

      // Logging
      if ( i % 100 == 0 ) {
        uint c = 0;
        for ( uint j = 0; j < clusterSize.size(); ++j ) if ( clusterSize[j] > 0 ) c++;
        cout << "\riter: " << setw(4) << iter
             << ", cell: " << setw(4) << i
             << ", n. clusters: " << setw(4) << c
             << ", sigma2:" << setw(8) << sigma2;
      }

      tmpBetaEigen = sigmaInvLLT.solve ( h.transpose() * y.row(i).transpose() / sigma2 );
      real tmp = log(alpha0) + 0.5 * p * log(tauBeta) + 0.5 * tmpBetaEigen.transpose() * sigmaInv * tmpBetaEigen
               + 0.5 * log(detSigma);
      real probNewClusterSum = 0.0;

      // Compute the probabilities of assigning to a new cluster or to an existing one
      // We store the probability of assigning to cluster j into probNewCluster[j], and into the last element we store
      // the probability of assigning to a new cluster
      vector<real> probNewCluster ( clusterSize.size(), 0 );

      for ( uint j = 0; j < clusterSize.size(); ++j ) { // For each cluster
         if ( clusterSize[j] == 0 ) {
            probNewCluster[j] = 0;
            continue;
         }

         real r = (y.row(i).transpose() - hBetaT.col(j)).squaredNorm();
         uint n = ( j == labels[i] && clusterSize[j] > 1 ? clusterSize[j] - 1 : clusterSize[j] );

         probNewCluster[j] = n * exp ( - r / (2.0 * sigma2) );
         probNewClusterSum += probNewCluster[j];
      }

      // The last entry in probNewCluster is the probability of creating a new cluster for the i-th cell
      probNewCluster.push_back( exp(tmp - y.row(i).squaredNorm() / (2*sigma2)) );
      probNewClusterSum += probNewCluster.back();

      // We remove the i-th cell is removed cluster
      for ( uint j = 0; j < cells; ++j ) {
        similarityMatrix(i,j) = 0;
        similarityMatrix(j,i) = 0;
      }
      clusterSize[labels[i]] -= 1;

      // Draw the new cluster for the i-th cell
      real which = runif ( 1, 0, probNewClusterSum )[0];
      uint newCluster = 0;

      // Find the right draw
      for ( uint k = 0; k < probNewCluster.size(); ++k ) {
        if ( which < probNewCluster[k] ) {
          newCluster = k;
          break;
        }
      }

      // If we drew the last element, we have to assign to a new cluster
      if ( newCluster >= clusterSize.size() ) {
         // It is possible that some previously created cluster was emptied. If so, we reuse the memory we allocated
         // for it to store the information about the newly created cluster

         // Look for an empty cluster
         uint k = 0;
         for ( ; k < clusterSize.size(); ++k )
            if ( clusterSize[k] == 0 ) break;

         // If there were no free slots, allocate a new one
         if ( k == clusterSize.size() ) {
            labels[i] = k;
            clusterSize.push_back(1);
            beta.conservativeResize( beta.rows()+1, Eigen::NoChange );
         }
         // Else, reuse the free slot
         else {
            labels[i] = k;
            clusterSize[k] = 1;
         }

         similarityMatrix(i,i) = 1;

         // Since we just created a new cluster, we have to sample its beta coefficients
         beta.row(labels[i]) = rmvnorm2 ( tmpBetaEigen, sigmaInvLLT );
         hBetaT = h * beta.transpose();
      }

      // If we drew an existing cluster, we assign the cell to it
      else {
         labels[i] = newCluster;
         clusterSize[newCluster] += 1;

         // Update the similarity matrix
         // We have to find another cell that belongs to this cluster, update the corresponding row in the similarity
         // matrix, and then copy that row into the row and into the column corresponding to the i-th cell
         for ( uint k = 0; k < cells; ++k ) {
            if ( labels[k] == newCluster ) {
               similarityMatrix(k,i) = 1;
               for ( uint l = 0; l < cells; ++l ) {
                  similarityMatrix(i,l) = similarityMatrix(k,l);
                  similarityMatrix(l,i) = similarityMatrix(i,l);
               }

               break;
            }
         }
      }
   } // End loop on cells

    // Step 2: resample the parameter for each of the clusters

    // First, we compute for each cluster the sum of Y vectors that are in that cluster
    MatrixXd sumY = MatrixXd::Zero ( clusterSize.size(), times );
    for ( uint i = 0; i < cells; ++i )
      sumY.row(labels[i]) += y.row(i);

    for ( uint c = 0; c < clusterSize.size(); ++c ) { // Loop on clusters
      // If the cluster is empty we clean up its coefficients (just for output readability purposes)
      if ( clusterSize[c] == 0 ) {
        beta.row(c) = MatrixXd::Zero ( 1, p );
        continue;
      }

      tmpSigmaInv = clusterSize[c] / sigma2 * hth;
      tmpSigmaInv.diagonal() += VectorXd::Constant(p,1,tauBeta);
      tmpSigmaInvLLT = tmpSigmaInv.llt();
      tmpBetaEigen = tmpSigmaInvLLT.solve ( h.transpose() * sumY.row(c).transpose() / sigma2 );
      beta.row(c) = rmvnorm2 ( tmpBetaEigen, tmpSigmaInvLLT );
    } // End loop on clusters

    // Recompute H times beta
    hBetaT = h * beta.transpose();

    // Step 3: resample the variance of residuals
    real sumres2 = 0;
    for ( uint i = 0; i < cells; ++i )
      sumres2 += (y.row(i).transpose() - hBetaT.col(labels[i])).squaredNorm();

    sigma2 = NumericVector(rinvgamma ( 1, _["shape"] = a0 + 0.5 * cells * times,
                                          _["rate"] = b0 + 0.5 * sumres2 ))[0];

    // Step 4: resample the variance of the coefficients' DP base measure
    uint clusterCount = 0; // Number of non-empty clusters
    tmpBetaEigen = VectorXd::Zero(p,1);
    for ( uint c = 0; c < clusterSize.size(); ++c ) {
      if ( clusterSize[c] == 0 ) continue;
      clusterCount++;
      tmpBetaEigen += beta.row(c).cwiseProduct ( beta.row(c) );
    }

    // Resample tauBeta
    tauBeta = NumericVector ( rgamma ( 1, nu0 + 0.5 * p * clusterCount, kappa0 + 0.5 * tmpBetaEigen.sum() ) )[0];

    // Recompute the matrices and vectors associated to sigma and sigmaBeta now that they have been resampled
    sigma1 = sqrt ( sigma2 );
    sigmaInv = hth / sigma2; sigmaInv.diagonal() += VectorXd::Constant(p,1,tauBeta);
    sigmaInvLLT = sigmaInv.llt();
    detSigma = 1 / sigmaInv.determinant();

    // Step 5: resample NA values.
    // NA values are treated as parameters of the model, and therefore they are resampled
    // at each Gibbs iteration
    for ( uint k = 0; k < naTimes.size(); ++k )
      y(naCells[k], naTimes[k]) = rnorm ( 1, hBetaT(naTimes[k],labels[naCells[k]]), sigma1 )[0];

    // Store current state
    if ( iter >= burnin && (iter - burnin) % trim == 0 ) {
      psm += similarityMatrix;
      betaRecord.push_back ( beta );
      labelsRecord.push_back ( labels );
      sigmaRecord.push_back ( sigma2 );
      tauBetaRecord.push_back ( tauBeta );
    }
  } // End MCMC loop

  // Output the record of state variables
  return List::create ( Named("psm") = psm / betaRecord.size(),
                        Named("beta") = betaRecord,
                        Named("labels") = labelsRecord,
                        Named("sigma") = sigmaRecord,
                        Named("tauBeta") = tauBetaRecord );
}
