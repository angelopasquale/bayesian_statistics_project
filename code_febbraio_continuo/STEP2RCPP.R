STEP2RCPP<-function(sigmaeps2,A,r,V,B,x,n_sites,W){
  for ( i in 1:n_sites ) {  
    Sigma_W = solveRcpp(1/sigmaeps2 * t(A) %*% A + diag(r))
    mu_W = 1/sigmaeps2 * Sigma_W %*% t(A) %*% (t(t(V[i,])) - B %*% x[i,])
    #Sigma_W <- Sigma_W + diag(ncol(Sigma_W))*0.01
    Sigma_W <- make.positive.definite(Sigma_W, tol=1e-3)
    # make.positive.definite finds the closest positive definite matrix using the
    # algorithm of NJ Higham
    # we can find an implemented version in Rcpp here:
    # https://stackoverflow.com/questions/51490499/results-for-calculating-nearest-positive-definite-matrix-are-different-in-r-func
    Sigma_W[lower.tri(Sigma_W)] = t(Sigma_W)[lower.tri(Sigma_W)]
    W[i,] = rmvnormRcpp ( n = 1, mu = mu_W, sigma = Sigma_W )
    
  }
  W
}