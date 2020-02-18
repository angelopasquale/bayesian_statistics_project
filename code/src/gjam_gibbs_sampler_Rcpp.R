gjam_gibbs_sampler_Rcpp<-function(alpha0,ndraws,burnin,N_stick,r,S,n_sites,x,Y){

tic("Parameters initialization from prior distributions:")
  
#eta_h parameter for Dz
eta_h <- 1/rgamma(r, shape = 1/2,  rate = 1/1e4  )
#Dz variance of Z
Dz<-riwish( 2 + r - 1, 4 * diag(1/eta_h))
#Construction of the matrix W
W<-matrix(0,nrow=n_sites,ncol=r)
for (j in seq(1,n_sites,1)) {
  W[j,] <- rmvnorm ( n = 1, mean = rep(0,times=r), sigma = diag(r) )
}
#sigmaeps2 variance of residuals
sigmaeps2 <- 1e-2
#k vector of labels
k<-rep(2,times=N_stick)
#Cardinality of species
Cardinality_S = rep ( 0, times=N_stick )
#V matrix of latent variables
V<-matrix(0,nrow=n_sites,ncol=S)
#p probabilities needed for stick breaking
p = rep(0,times=N_stick)
#nu and G needed for full conditional of sigmaeps2
nu = 0
G = 1
#Q matrix needed for computation of A
Q=matrix(0,nrow=S,ncol=N_stick)
#prior matrix Z
Z<-rmvnormRcpp(N_stick,rep(0,r),1/S*diag(r))


#Generalized dirichlet prior and stick breaking
xi <- compute_GD_prior(N_stick,alpha0,k)
p[1] <- xi[1]
p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
p[N_stick] <- 1 - sum( p[1:N_stick-1] )
k = sample(1:N_stick, size = S, replace=TRUE, prob = p)
pl <- matrix(0, nrow = S, ncol = N_stick)  
logpl <- matrix(0, nrow = S, ncol = N_stick)
# number of covariates
n_cov = ncol(x)  
# Prior mean of the beta coefficients
muBeta = rep ( 0, times=n_cov ) 
# Prior variance-covariance matrix of the beta coefficients
sigmaBeta = diag ( n_cov ) 
# B coefficient matrix
B <- matrix(data = 0, nrow = S, ncol = n_cov) 
for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 10*sigmaBeta )
}
#Varianca for full conditional of B
sigmaB = 10


#Lists for chains
list_B <- list()
list_A <- list()
list_Z <- list()
list_k <- list()
list_R <- list()
list_sigmaeps2 <- list()
toc()
timer1<-0
timer2<-0
timer3<-0
timer4<-0
timer5<-0
timer6<-0
timer7<-0
timer8<-0
tic("Gibbs loop")
#Gibbs sampler
for(i in 1:ndraws){
  
  ######################### STEP 1 ##########################
  ####################### posterior of A ####################
  
  step1result<-STEP1RCPP(N_stick,r,sigmaeps2,Z,Q,k,x,Dz,W,V,B,eta_h,Y)
  A<-step1result$A
  timer1<-step1result$timer+timer1
  
  
  
  ######################### STEP 2 ##########################
  ####################### posterior of W ####################
  
  step2result<-STEP2RCPP(sigmaeps2,A,r,V,B,x,n_sites,W)
  W<-step2result$W
  timer2<-timer2+step2result$timer
  
  
  
  ######################### STEP 3 ##########################
  ####################### posterior of A ####################
 
  
  step3result<-STEP3RCPP(p,V,Z,x,B,W,sigmaeps2,N_stick)
  k<-step3result$k
  timer3<-timer3+step3result$timer
  
  
  
  ######################### STEP 4 ##########################
  ####################### posterior of p ####################
  step4result<-STEP4RCPP(p,N_stick,S,alpha0,k)
  p<-step4result$p
  timer4<-step4result$timer4+timer4
  
  
  
  ######################### STEP 5 ##########################
  ####################### posterior of sigmaeps2 ####################
  step5result<-STEP5RCPP(n_sites,V,B,x,A,W,S,nu,G,sigmaeps2,Z,k,r,Y)
  sigmaeps2<-step5result$sigmaeps2
  timer5<-timer5+step5result$timer
  
  
  
  ######################### STEP 6 ##########################
  ####################### posterior of Dz ####################
  step6result<-STEP6(r,Dz,Z,N_stick)
  Dz<-step6result$Dz
  timer6<-timer6+step6result$timer
  
  
  
  ######################### STEP 7 ##########################
  ####################### posterior of V ####################
  step7result<-STEP7RCPP(A,sigmaeps2,S,B,n_sites,W,Y,x)
  V<-step7result$V
  V_star<-step7result$V_star
  B_star=step7result$B_star
  D=step7result$D
  timer7<-step7result$timer+timer7
  
  
  ######################### STEP 8 ##########################
  ####################### posterior of B ####################
  step8result<-STEP8RCPP(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star)
  B<-step8result$B
  timer8<-timer8+step8result$timer
  
  
  #Lists of posteriors
  list_B[[i]] <- B
  list_A[[i]] <- A
  list_Z[[i]] <- Z
  list_k[[i]] <- k
 # list_R[[niter]] <- R
  list_sigmaeps2[[i]] <- sigmaeps2
}
toc()

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
A_sup<-apply(simplify2array(list_A),1:2,quantile,0.95)
A_inf<-apply(simplify2array(list_A),1:2,quantile,0.05)
# EE<-matrix(0,nrow=S,ncol=r)
# for(r in(burnin:ndraws)){
#   EE<-list_A[[r]]
#   if(dim(uniquecombs(EE))[1]=="NULL"){bp[[r]]<-0}
#   else{
#   bp[[r]]<-dim(uniquecombs(EE))[1]}
# }
# bp<-bp[burnin:ndraws]

#bp<-bp[burnin:ndraws]
#What you take back in the main
#h<-list("A_media"=A_media,"A_true"=data$A_true,"A_inf"=A_inf,"A_sup"=A_sup,"x"=x,"Y"=Y,"bp"=bp)
#return(chain)
#return(h)
time_list<-list("timer1"=timer1,"timer2"=timer2,"timer3"=timer3,"timer4"=timer4,"timer5"=timer5,"timer6"=timer6,"timer7"=timer7,"timer8"=timer8)
return(list("time_list"=time_list,"chain"=chain,"A_sup"=A_sup,"A_inf"=A_inf,"A_media"=A_media))
}
  