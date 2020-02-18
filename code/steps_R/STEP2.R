STEP2<-function(sigmaeps2,A,r,V,B,x,n_sites,W){
  tic("Step 2")
  for ( i in 1:n_sites ) {  
    Sigma_W = solve(1/sigmaeps2 * t(A) %*% A + diag(r))
    Sigma_W <- make.positive.definite(Sigma_W, tol=1e-3)
    mu_W = 1/sigmaeps2 * Sigma_W %*% t(A) %*% (t(t(V[i,])) - B %*% x[i,])
    Sigma_W[lower.tri(Sigma_W)] = t(Sigma_W)[lower.tri(Sigma_W)]
    W[i,] <- rmvnorm ( n = 1, mean = mu_W, sigma = Sigma_W )
  }
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("W"=W,"timer"=timings))
}
