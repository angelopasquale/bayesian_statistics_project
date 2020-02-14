STEP1<-function(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B){
  for ( j in 1:N_stick ) {
    if (!(j %in% k)) {
      #print('entered if')
      Z[j,] <- rmvnorm ( n = 1, mean = rep(0, times = r), sigma = Dz )
      Cardinality_S[j]<- 1
      #problema: k non viene modificato, se ho un nuovo label aumento
      #la cardinalità del label rispettivo ma non inserisco nel k il nuovo label
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
      # Sigma_Zj <- Sigma_Zj + diag(ncol(Sigma_Zj))*0.01
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
  A
}