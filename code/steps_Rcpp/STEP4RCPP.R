STEP4RCPP<-function(p,N_stick,S,alpha0,k){
  # Step 4 -> R
  # RESAMPLE p|k
  # OUR CODE
  # xi <- compute_GD_prior(N_stick,alpha0,k)
  # # c <- 2
  # # given the vector k = [1,..., 1] of length = S
  # p[1] <- xi[1]
  # p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
  # p[N_stick] <- 1 - sum( p[1:N_stick-1] )
  # BUT in Rcpp:
  p <- .sampleP(N = N_stick, avec = rep(alpha0/N_stick,(N_stick-1)),
                bvec = ((N_stick-1):1)*alpha0/N_stick, K = k) 
  
}