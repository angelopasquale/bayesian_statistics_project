STEP5<-function(n_sites,V,B,x,A,W,S,nu,G){
  tic("Step 5")
  dx = c(0)
  for (i in 1:n_sites) {
    dx = dx + norm(V[i,] - B %*% x[i,] - A %*% W[i,])^2
  }
  sigmaeps2<-1/rgamma(1,shape = (n_sites * S + nu)/2, rate = dx/2 + nu/G^2)
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("sigmaeps2"=sigmaeps2,"timer"=timings))
}