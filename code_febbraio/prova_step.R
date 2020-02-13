prova_step<-function(alpha0,niter,N_stick,r,S){
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

#Questa funzione genera la matrice x delle covariate e la matrice B dei beta
simulation_fun<-function(Sp=S,nsamples=n_sites, r=4, K_t=4){
  S<-Sp
  n<- nsamples
  #iterations<-it
  
  #create the design matrix: here we use the intercept and one covariate (scaled!)
  X<-cbind(rep(1,n),scale(runif(0,100,n=n))) 
  
  #create the coefficient matrix B (Sx2)
  idx<-sample(S)
  B_0<-scale(seq(0,100,length.out=S)[idx]) #intercept
  idx<-sample(S)
  B_1<-scale(seq(0,100,length.out=S)[idx]) #covariate coefficient
  B<-cbind(B_0,B_1) #Coefficient matrix
  L<-X%*%t(B) #We create the mean by multiplying B with the design matrix X
  
  
  #We create A and then Sigma
  A<-matrix(NA,nrow=K_t,ncol=r) #A is the matrix with the atoms only. Initialization.
  sig=matrix(runif(n=r*r),ncol=r) #sig is the variance covariance matrix of Z.
  for(i in 1:K_t){
    A[i,]<-mvrnorm(n = 1, rep(0,r), Sigma=1*diag(r)) # We sample the unique values of A
  }
  idx<-sample((1:K_t),S,replace=T) #idx represents the labels
  Lambda<-A[idx,] #We expand A using the labels to obtain Lambda (what we used to call A)
  Sigma_true<-Lambda%*%t(Lambda)+0.1*diag(S) #We obtain Sigma. Here sigma_epsilon^2 is 0.1
  
  
  #We obtain the correlation matrix R
  R_true=cov2cor(Sigma_true) 
  
  #We sample V from the model
  Y_cont<-L+rmvnorm(n = n, mean=rep(0,S), sigma=R_true)
  # We obtain our binary dataset. Secondo me vi conviene prima testare il modello con Y_cont e poi su Y
  Y<- ifelse(Y_cont>0,1,0)
  
  give_back<-list(A_true = Lambda, B_true=B, Xdesign=X, mu_true=L, R_true=R_true, V=Y_cont, Y=Y)
  return(give_back)
}
data<-simulation_fun()
x<-data$Xdesign
#initialization of B
#B<-data$B_true
n_cov = ncol(x) # number of covariates (no intercept)
muBeta = rep ( 0, times=n_cov ) # Prior mean of the beta coefficients
sigmaBeta = diag ( n_cov ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = n_cov) # coefficient matrix
for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 10*sigmaBeta )
}
sigmaB = 10
Y<-data$Y
eta_h<-1/100

list_B <- list()
list_A <- list()
list_Z <- list()
list_k <- list()
list_R <- list()
list_sigmaeps2 <- list()

#Lancio gli steps
for(i in 1:niter){
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
 lista<-STEP7(A,sigmaeps2,S,B,n_sites,W,Y,x)
  V<-lista$V
 # V1<-STEP7RCPP(A,sigmaeps2,S,B,n_sites,W,Y,x)
  V_star<-lista$V_star
  B_star=lista$B_star
  D=lista$D
  #8
 B<-STEP8(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star)
  #B1<-STEP8RCPP(S,x,sigmaeps2,V_star,W,A,sigmaB,n_cov,D,B_star)
  list_B[[i]] <- B
  list_A[[i]] <- A
  list_Z[[i]] <- Z
  list_k[[i]] <- k
 # list_R[[niter]] <- R
  list_sigmaeps2[[i]] <- sigmaeps2
}
h<-list()
for(i in (1:niter)){
  h[[i]]<-list_A[[i]][2,2]
}
h<-as.numeric(h)
h<-as.data.frame(h)
h<-as.mcmc(h)
h<-ggs(h)
x11()
ggs_traceplot(h) + 
  geom_area(colour="blue") + xlab('Year') + ylab('Activity ') + theme(text = element_text(size=15, family="LM Roman 10")) +
  labs(fill = 'Technologies')+ ggtitle('Baseline - Activity') 
x11()
ggs_running(h)
h
}
  