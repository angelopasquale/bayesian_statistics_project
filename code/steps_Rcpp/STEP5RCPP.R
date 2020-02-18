STEP5RCPP<-function(n_sites,V,B,x,A,W,S,nu,G,sigmaeps2,Z,k,r,Y){
  # Step 5 
  tic("Step 5")
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
   toc(log=TRUE,quiet=TRUE)
   log.txt <- tic.log(format = TRUE)
   log.lst <- tic.log(format = FALSE)
   tic.clearlog()
   timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
   return(list("sigmaeps2"=sigmaeps2,"timer"=timings))
}
