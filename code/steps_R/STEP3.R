STEP3<-function(S,N_stick,sigmaeps2,p,pl,logpl,V,B,W,Z,k){
  tic("Step 3")
  logpl<-matrix(0,nrow=S,ncol=N_stick)
  for ( l in 1:S ) { # Loop on species
    for (j in 1:N_stick) {
      logpl[l,j] <- log(p[j]) -1/(2*sigmaeps2) * norm(V[,l] - x %*% B[l,] - W %*% Z[j,])^2
    }
    logpl[l,]<- logpl[l,] - max(logpl[l,])   
    pl[l,]<- exp(logpl[l,])
    pl[l,]<- pl[l,] / sum(pl[l,])
    k[l]<- sample(N_stick, size = 1, replace=TRUE, prob = pl[l,])
  }
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("k"=k,"timer"=timings,"pl"=pl,"logpl"=logpl))
}