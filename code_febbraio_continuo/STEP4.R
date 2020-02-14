STEP4<-function(p,N_stick,S,alpha0){
  xi <- compute_GD_prior(N_stick,alpha0,k)
  # c <- 2
  # given the vector k = [1,..., 1] of length = S
  p[1] <- xi[1]
  p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
  p[N_stick] <- 1 - sum( p[1:N_stick-1] )
  # k = sample(1:N_stick, size = S, replace=TRUE, prob = p) #why we sample again k here?
 # pl <<- matrix(0, nrow = S, ncol = N_stick) # needed for step 4 in the Gibbs 
  #logpl <<- matrix(0, nrow = S, ncol = N_stick
  p
}