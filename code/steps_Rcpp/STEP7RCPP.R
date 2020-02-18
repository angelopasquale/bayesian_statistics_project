STEP7RCPP<-function(A,sigmaeps2,S,B,n_sites,W,Y,x){
  Sigma_star = A %*% t(A) + sigmaeps2 * diag(S)
  D = diag(diag(Sigma_star))
  R = solveRcpp(D)^(1/2) %*% Sigma_star %*% solveRcpp(D)^(1/2)
  B_star = (D)^(1/2) %*% B
  V_star=matrix(0,nrow=n_sites,ncol=S)
  
  for (i in seq(1,n_sites)) {
    for (j in seq(1,S)) {
      
      mean1 = as.vector(B_star[j,] %*% x[i,] + A[j,] %*% W[i,])
      
      if (Y[i,j] == 1){
        V_star[i,j] <- .tnorm(1, 0, Inf, mean1, sigmaeps2)
      } else {
        V_star[i,j] <- .tnorm(1, -Inf, 0, mean1, sigmaeps2)
      }
    }
  }
  #print(dim(solveRcpp(D)^(1/2)))
 # print(dim(V_star))
  V = V_star %*% solveRcpp(D)^(1/2)
  newList <- list("V" = V, "V_star" = V_star, "D"=D,"B_star"=B_star)
  return(newList)
  
}