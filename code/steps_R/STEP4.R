STEP4<-function(p,N_stick,S,alpha0,k){
  tic("Step 4")
  xi <- compute_GD_prior(N_stick,alpha0,k)
  p[1] <- xi[1]
  p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
  p[N_stick] <- 1 - sum( p[1:N_stick-1] )
  pl <- matrix(0, nrow = S, ncol = N_stick) 
  logpl <- matrix(0, nrow = S, ncol = N_stick)
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("p"=p,"timer"=timings,"pl"=pl,"logpl"=logpl))
}