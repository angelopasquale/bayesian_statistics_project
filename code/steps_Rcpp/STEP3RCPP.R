STEP3RCPP<-function(p,V,Z,x,B,W,sigmaeps2,N_stick){
  ###STEP 3###
  #This function computes vector of labels a posteriori
  tic("Step 3")
  pl <- getPmatKRcpp(pveck = p,Yk = V, Zk = Z,
                     Xk = x, Bk = B, Wk = W,
                     sigmasqk = sigmaeps2)
  # and then the vector of labels k
  k <- unlist( apply(pl, 1, function(x)sample(1:N_stick, size=1, prob=x)) )
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("k"=k,"timer"=timings))
}