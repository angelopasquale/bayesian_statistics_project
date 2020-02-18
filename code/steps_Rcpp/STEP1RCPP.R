STEP1RCPP<-function(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B,eta_h,Y){
  ##### Step 1 #####
  # This step computes matrix A
  # resample the cluster assignments
  tic("Step 1")
  nn   <- nrow(x)
  pp    <- ncol(x)
  SS    <- ncol(V)
  ntot <- nrow(V)
  
  covR <- solveRcpp( (1/sigmaeps2)*crossprod(Z[k,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[k,]/sigmaeps2,t(Y - x %*% t(B) ))
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  Dz<-Dz+0.001*diag(r)
  Z    <- fnZRcpp(kk=k, Yk=V, Xk=x, Dk=Dz, Bk=B, 
                  Wk=RR, sigmasqk=sigmaeps2, Nz=N_stick)
  
  A = Z[k,]
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("A"=A,"timer"=timings))
}