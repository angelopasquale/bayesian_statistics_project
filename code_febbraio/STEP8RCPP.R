STEP8RCPP<-function(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star){
  # Step 8 : TO CHECK, and REALLY CAREFULLY !!!!!!!!! to do in Rcpp?
  for (j in seq(1,S,1)) {
    muBetaj = solveRcpp(1/(sigmaB)^2 * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x) %*% t(x) %*% (V_star[,j] - W %*% A[j,]) * 1/sigmaeps2
    sigmaBetaj = solveRcpp(1/(sigmaB^2) * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x)
    # sigmaBetaj <- sigmaBetaj + diag(ncol(sigmaBetaj))*0.01
    sigmaBetaj <- make.positive.definite(sigmaBetaj, tol=1e-3)
    sigmaBetaj[lower.tri(sigmaBetaj)] = t(sigmaBetaj)[lower.tri(sigmaBetaj)]
    B_star[j,] <- rmvnormRcpp ( n = 1, mu = muBetaj, sigma = sigmaBetaj )
  }
  
  B = solveRcpp(D)^(1/2) %*% B_star;
 B
}