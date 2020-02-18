STEP6<-function(r,Dz,Z,N_stick){
  tic("Step 6")
  eta_h <- 1/rgamma(r, shape = (2 + r )/2, 
                    rate = ((1/1000) + 2*diag(solveRcpp(Dz+0.001*diag(r))) ) )
  S_Dz <- crossprod(Z) + 4 * diag(1/eta_h)
  Dz<-riwish( 2 + r + N_stick - 1, S_Dz)
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("Dz"=Dz,"timer"=timings))
}