STEP8<-function(S,x,sigmaeps2,Y,W,A,sigmaB,n_cov,D){
  B<-matrix(0,nrow=S,ncol=2)
for (j in seq(1,S,1)) {
  muBetaj = solve(1/(sigmaB)^2 * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x) %*% t(x) %*% (Y[,j] - W %*% A[j,]) * 1/sigmaeps2
  sigmaBetaj = solve(1/(sigmaB^2) * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x)
  # sigmaBetaj <- sigmaBetaj + diag(ncol(sigmaBetaj))*0.01
  sigmaBetaj <- make.positive.definite(sigmaBetaj, tol=1e-3)
  sigmaBetaj[lower.tri(sigmaBetaj)] = t(sigmaBetaj)[lower.tri(sigmaBetaj)]
  B[j,] <- rmvnorm ( n = 1, mean = muBetaj, sigma = sigmaBetaj )
}
B
#B = solve(D)^(1/2) %*% B_star;

#list_B[[niter]] <- B
#list_A[[niter]] <- A
#list_Z[[niter]] <- Z
#list_k[[niter]] <- k
#list_R[[niter]] <- R
#list_sigmaeps2[[niter]] <- sigmaeps2
} 
