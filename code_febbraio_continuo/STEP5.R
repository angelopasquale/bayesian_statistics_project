STEP5<-function(n_sites,V,B,x,A,W,S,nu,G){
  dx = c(0)
  for (i in 1:n_sites) {
    dx = dx + norm(V[i,] - B %*% x[i,] - A %*% W[i,])^2
  }
  sigmaeps2<-1/rgamma(1,shape = (n_sites * S + nu)/2, rate = dx/2 + nu/G^2)
  sigmaeps2
}