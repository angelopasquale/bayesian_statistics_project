STEP7<-function(A,sigmaeps2,S,B,n_sites,W,Y,x){
Sigma_star = A %*% t(A) + sigmaeps2 * diag(S)
D = diag(diag(Sigma_star))
R = solve(D)^(1/2) %*% Sigma_star %*% solve(D)^(1/2)
B_star = D^(1/2) %*% B
V_star=matrix(0,nrow=n_sites,ncol=S)
#Y 400*15
for (i in (1:n_sites)) {
  #print("ghgh")
  for (j in (1:(S))) {
    #print("j")
    mean1 =  B_star[j,]%*%x[i,] + A[j,] %*% W[i,]
    if (Y[i,j]== 1){
      V_star[i,j] <- tnorm(1, 0, Inf, mean1, sigmaeps2)
    } else {
      V_star[i,j] <- tnorm(1, -Inf, 0, mean1, sigmaeps2)
    }
  }
}

V <- V_star %*% solve(D)^(1/2)
newList <- list("V" = V, "V_star" = V_star, "D"=D,"B_star"=B_star)
return(newList)
}

# L<-x%*%t(B) #We create the mean by multiplying B with the design matrix X
# 
# Sigma_star <-A%*%t(A)+sigmaeps2*diag(S) #We obtain Sigma. Here sigma_epsilon^2 is 0.1
# 
# 
# #We obtain the correlation matrix R
# B_star =cov2cor(Sigma_star) 
# V<-L+rmvnormRcpp(n = n, mu=rep(0,S), sigma=B_star)
