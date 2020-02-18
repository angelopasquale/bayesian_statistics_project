STEP4RCPP<-function(p,N_stick,S,alpha0,k){
  tic("Step 4")
  #This step computes probabilities p a posteriori
  p <- .sampleP(N = N_stick, avec = rep(alpha0/N_stick,(N_stick-1)),
                bvec = ((N_stick-1):1)*alpha0/N_stick, K = k) 
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("p"=p,"timer"=timings))
  
}