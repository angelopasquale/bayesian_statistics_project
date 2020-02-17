fit_gjam_gibbs<-function(alpha0,ndraws,burnin,N_stick,r,S,n_sites,x,Y){

eta_h <- 1/rgamma(r, shape = 1/2,  rate = 1/1e4  )
Dz<-riwish( 2 + r - 1, 4 * diag(1/eta_h))
W<-matrix(0,nrow=n_sites,ncol=r)
for (j in seq(1,n_sites,1)) {
  W[j,] <- rmvnorm ( n = 1, mean = rep(0,times=r), sigma = diag(r) )
}
sigmaeps2 <- 1e-2
k<-rep(2,times=N_stick)
Cardinality_S = rep ( 0, times=N_stick )
V<-matrix(0,nrow=n_sites,ncol=S)
p = rep(0,times=N_stick)
pl=matrix(0,nrow=S,ncol=N_stick)
nu = 0
G = 1

Q=matrix(0,nrow=S,ncol=N_stick)
Z<-rmvnormRcpp(N_stick,rep(0,r),1/S*diag(r))
logpl<-matrix(0,nrow=S,ncol=N_stick)

#Inizializzazione di xi, p, k, pl, logpl per lo step 3 ( inizializzati come zeri dava too few
#positive probabilities)
xi <- compute_GD_prior(N_stick,alpha0,k)
# given the vector k = [1,..., 1] of length = S
p[1] <- xi[1]
p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
p[N_stick] <- 1 - sum( p[1:N_stick-1] )
k = sample(1:N_stick, size = S, replace=TRUE, prob = p)
pl <- matrix(0, nrow = S, ncol = N_stick)  
logpl <- matrix(0, nrow = S, ncol = N_stick)



n_cov = ncol(x) # number of covariates 
muBeta = rep ( 0, times=n_cov ) # Prior mean of the beta coefficients
sigmaBeta = diag ( n_cov ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = n_cov) # coefficient matrix
for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 10*sigmaBeta )
}
sigmaB = 10
eta_h<-1/100

list_B <- list()
list_A <- list()
list_Z <- list()
list_k <- list()
list_R <- list()
list_sigmaeps2 <- list()

#Gibbs sampler
for(i in 1:ndraws){
  #1
  #A<-STEP1(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B)
  #A<-L$A
  #Cardinality_S<-L$Cardinality_S
  A<-STEP1RCPP(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B,eta_h,Y)
  
  #2
  #W<-STEP2(sigmaeps2,A,r,V,B,x,n_sites,W)
  W<-STEP2RCPP(sigmaeps2,A,r,V,B,x,n_sites,W)
  #3
  #STEP3(S,N_stick,sigmaeps2,p,pl,logpl)
  k<-STEP3RCPP(p,V,Z,x,B,W,sigmaeps2,N_stick)
  #if(i==19){k_1<-k}
  #4
  #p<-STEP4(p,N_stick,S,alpha0)
  p<-STEP4RCPP(p,N_stick,S,alpha0,k)
  #5
 #sigmaeps2<-STEP5(n_sites,V,B,x,A,W,S,nu,G)
 sigmaeps2<-STEP5RCPP(n_sites,V,B,x,A,W,S,nu,G,sigmaeps2,Z,k,r,Y)
  #6
Dz<-STEP6(r,Dz,Z,N_stick)
 #7
 lista<-STEP7RCPP(A,sigmaeps2,S,B,n_sites,W,Y,x)
  V<-lista$V
 # V1<-STEP7RCPP(A,sigmaeps2,S,B,n_sites,W,Y,x)
  V_star<-lista$V_star
  B_star=lista$B_star
  D=lista$D
  #8
  #B<-STEP8TaylorRCPP(S,x,sigmaeps2,V,W,A,sigmaB,n_cov,D,B)
  B<-STEP8RCPP(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star)
  list_B[[i]] <- B
  list_A[[i]] <- A
  list_Z[[i]] <- Z
  list_k[[i]] <- k
 # list_R[[niter]] <- R
  list_sigmaeps2[[i]] <- sigmaeps2
}

#Here we get the chain for one element of matrix A
chain<-list()
for(i in (burnin:ndraws)){
  chain[[i]]<-list_A[[i]][1,2]
chain}
#Cutting away burnin NULL values in the chain
chain<-chain[burnin:ndraws]


#Here is created mean matrix A of the chain
A_media=matrix(0,nrow=S,ncol=r)
for(i in (1:S)){
  for(j in (1:r)){
    tot=0
    for(k in (burnin:ndraws)){
      tot=tot+list_A[[k]][i,j]
    }
    A_media[i,j]=tot/(ndraws-burnin)
  }
}

#A = return_list$A
bp<-list()
A_sup<<-apply(simplify2array(list_A),1:2,quantile,0.95)
A_inf<<-apply(simplify2array(list_A),1:2,quantile,0.05)
EE<-matrix(0,nrow=S,ncol=r)
for(r in(burnin:ndraws)){
  EE<-list_A[[r]]
  bp[[r]]<-dim(uniquecombs(EE))[1]
}
bp<-bp[burnin:ndraws]

#bp<-bp[burnin:ndraws]
#What you take back in the main
#h<-list("A_media"=A_media,"A_true"=data$A_true,"A_inf"=A_inf,"A_sup"=A_sup,"x"=x,"Y"=Y,"bp"=bp)
return(chain)
#return(h)
}
  