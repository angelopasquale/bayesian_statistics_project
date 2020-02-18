STEP1<-function(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B,Cardinality_S){
  tic("Step 1")
  for ( j in 1:N_stick ) {
    if (!(j %in% k)) {
      Z[j,] <- rmvnorm ( n = 1, mean = rep(0, times = r), sigma = Dz )
      Cardinality_S[j]<- 1
    }
    else { # if j in k
      Cardinality_S[j]<- sum(j==k)
      Sigma_Zj = solve(Cardinality_S[j]/sigmaeps2 * t(W) %*% W + solve(Dz) )
      mu_Zj = c(0)
      for (l in 1:S) {
        if (k[l] == j){
          mu_Zj = mu_Zj + 1/sigmaeps2 * Sigma_Zj %*% t(W) %*% ( t(t(V[,l])) - x %*% B[l,] )
        }
        Q[l,k[l]]<- 1
      }
      Sigma_Zj <- make.positive.definite(Sigma_Zj, tol=1e-3)
      # because of machine precision the matrix could seem to be not positive definite:
      # we overcome this by making sure your det(Sigma_Zj) returns a positive
      # value. One way is to add some variance in all directions:
      Sigma_Zj[lower.tri(Sigma_Zj)] = t(Sigma_Zj)[lower.tri(Sigma_Zj)]
      # because of machine precision the matrix doesn't seem symmetric
      Z[j,]<-rmvnorm ( n = 1, mean = mu_Zj, sigma = Sigma_Zj )
    }
  }
  A<- Q %*% Z
  toc(log=TRUE,quiet=TRUE)
  log.txt <- tic.log(format = TRUE)
  log.lst <- tic.log(format = FALSE)
  tic.clearlog()
  timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
  return(list("A"=A,"timer"=timings))
}