STEP5RCPP<-function(n_sites,V,B,x,A,W,S,nu,G,sigmaeps2,Z,k,r,Y){
  # Step 5 in R, todo in Rcpp
  
  nn   <- nrow(x)
   pp    <- ncol(x)
   SS    <- ncol(V)
   ntot <- nrow(V)
   
   covR <- solveRcpp( (1/sigmaeps2)*crossprod(Z[k,]) + diag(r) ) # Sigma_W
   z1   <- crossprod( Z[k,]/sigmaeps2,t(Y - x %*% t(B) ))        
   RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  # 
   rndEff <- RR%*%t(Z[k,])
   
   res        <- sum((V - x %*% t(B) - rndEff )^2)
   sigmaeps2 <- 1/rgamma(1,shape=(SS*nn + 1)/2, rate=res/2)  
   #print(sigmaeps2)
   sigmaeps2
}
