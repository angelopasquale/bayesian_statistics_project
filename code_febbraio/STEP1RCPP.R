STEP1RCPP<-function(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B,eta_h,Y){
  ##### Step 1 #####
  #This step computes matrix A
  #print("-----------------------------")
  #print("------------STEP1------------")
  # resample the cluster assignments
  # for each species j = 1,...,S do
  # Z <- matrix(data = 0, nrow = N_stick, ncol = r) 
  #Z <- rmvnormRcpp(N_stick,rep(0,r),1/S*diag(r))
  # Q <- matrix(data = 0, nrow = S, ncol = N_stick) 
  
  nn   <- nrow(x)
  pp    <- ncol(x)
  SS    <- ncol(V)
  ntot <- nrow(V)
  
  covR <- solveRcpp( (1/sigmaeps2)*crossprod(Z[k,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[k,]/sigmaeps2,t(Y - x %*% t(B) ))
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  Dz<-Dz+0.001*diag(r)
  #D    <- .riwish(df = (2 + r + N_stick), S = (crossprod(Z) + 2*2*diag(1/eta_h)))
  Z    <- fnZRcpp(kk=k, Yk=V, Xk=x, Dk=Dz, Bk=B, 
                  Wk=RR, sigmasqk=sigmaeps2, Nz=N_stick)
  
  A = Z[k,]
  A
}