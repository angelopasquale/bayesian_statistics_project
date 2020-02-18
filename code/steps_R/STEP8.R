STEP8<-function(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star){
  tic("Step 8")
for (j in seq(1,S,1)) {
  muBetaj = solve(1/(sigmaB)^2 * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x) %*% t(x) %*% (V_star[,j] - W %*% A[j,]) * 1/sigmaeps2
  sigmaBetaj = solve(1/(sigmaB^2) * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x)
  # sigmaBetaj <- sigmaBetaj + diag(ncol(sigmaBetaj))*0.01
  sigmaBetaj <- make.positive.definite(sigmaBetaj, tol=1e-3)
  sigmaBetaj[lower.tri(sigmaBetaj)] = t(sigmaBetaj)[lower.tri(sigmaBetaj)]
  B_star[j,] <- rmvnorm ( n = 1, mean = muBetaj, sigma = sigmaBetaj )
}

B = solve(D)^(1/2) %*% B_star
toc(log=TRUE,quiet=TRUE)
log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
return(list("B"=B,"timer"=timings))

} 
