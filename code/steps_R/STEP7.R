STEP7<-function(A,sigmaeps2,S,B,n_sites,W,Y,x){
  tic("Step 7")
  Sigma_star = A %*% t(A) + sigmaeps2 * diag(S)
  D = diag(diag(Sigma_star))
  R = solve(D)^(1/2) %*% Sigma_star %*% solve(D)^(1/2)
  B_star = (D)^(1/2) %*% B
  V_star=matrix(0,nrow=n_sites,ncol=S)
  
  for (i in seq(1,n_sites)) {
    for (j in seq(1,S)) {
      
      mean1 = as.vector(B_star[j,] %*% x[i,] + A[j,] %*% W[i,])
     
      if (Y[i,j] == 1){
        V_star[i,j] <- tnorm(1, 0, Inf, mean1, sigmaeps2)
      } else {
        V_star[i,j] <- tnorm(1, -Inf, 0, mean1, sigmaeps2)
      }
    }
  }
  
  V =  V_star%*%solve(D)^(1/2)
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("V" = V, "V_star" = V_star, "D"=D,"B_star"=B_star,"timer"=timings))
}

