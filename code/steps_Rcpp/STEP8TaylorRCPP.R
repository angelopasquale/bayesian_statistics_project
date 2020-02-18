STEP8TaylorRCPP<-function(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star){
  # Step 8 : TO CHECK, and REALLY CAREFULLY !!!!!!!!! to do in Rcpp?
  for (j in seq(1,S,1)) {
    muBetaStarj = t( V_star[j,] - W %*% A[j,] ) %*% x %*% solveRcpp(t(x) %*% x)
    sigmaBetaStar = sigmaeps2 * solveRcpp(t(x) %*% x)
    B_star[j,] <- rmvnormRcpp ( n = 1, mu = muBetaStarj, sigma = sigmaBetaStar )
  }
  
  B = solveRcpp(D)^(1/2) %*% B_star;
}