STEP6<-function(r,Dz,Z,N_stick){
  eta_h <- 1/rgamma(r, shape = (2 + r )/2, 
                    rate = ((1/1000) + 2*diag(solveRcpp(Dz+0.001*diag(r))) ) )
  #eta_h <- 1/rgamma(r, shape = 1/2,  rate = 1/1e4  )
  
  S_Dz <- crossprod(Z) + 4 * diag(1/eta_h)
  # S_Dz <- S_Dz + diag(ncol(S_Dz))*0.01
  #S_Dz <- make.positive.definite(S_Dz, tol=1e-3)
  #S_Dz[lower.tri(S_Dz)] = t(S_Dz)[lower.tri(S_Dz)]
  Dz<-riwish( 2 + r + N_stick - 1, S_Dz)
  Dz
}