STEP3<-function(S,N_stick,sigmaeps2,p,pl,logpl){
  logpl<-matrix(0,nrow=S,ncol=N_stick)
  for ( l in 1:S ) { # Loop on species
    for (j in 1:N_stick) {
      pl[l,j] <- p[j] * exp( -1/(2*sigmaeps2) * norm(V[,l] - x %*% B[l,] - W %*% Z[j,])^2)
    }
    logpl[l,]<- log(pl[l,])   
    logpl[l,]<- logpl[l,] - max(logpl[l,])   
    pl[l,]<- exp(logpl[l,])
    pl[l,]<- pl[l,] / sum(pl[l,])
    k[l]<- sample(N_stick, size = 1, replace=TRUE, prob = pl[l,])
  }
}