# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())
# Clear console
cat("\014") 

needed_packages <- c("MASS","coda","ggmcmc","extrafont","mgcv","mvtnorm", "matlib", "devtools", "MCMCpack", "gjam", "invgamma", "MixMatrix", "tictoc", "corpcor") 
new_packages <- needed_packages[!(needed_packages %in% installed.packages()[, "Package"])] 
if (length(new_packages)) install.packages(new_packages) 
lapply(needed_packages, require, character.only = TRUE) 

# Install **TTF** Latin Modern Roman fonts from www.fontsquirrel.com/fonts/latin-modern-roman
# Import the newly installed LModern fonts, change the pattern according to the 
# filename of the lmodern ttf files in your fonts folder
# font_import(pattern = "lmodern*")
# loadfonts(device = "win")
# par(family = "LM Roman 10")

# required source functions
Rcpp::sourceCpp('src/cpp/cppFns.cpp')
for (i in seq(1,8)) {
  source(paste("steps_R/STEP",i,".R", sep = ""))
}
for (i in seq(1,8)) {
  source(paste("steps_Rcpp/STEP",i,"RCPP.R", sep = ""))
}
source("steps_Rcpp/STEP8TaylorRCPP.R")
source("src/compute_GD_prior.R")
source("src/tnorm.R")
source("src/gjamHfunctions.R")
source("src/check_CR.R")
source("src/gjam_gibbs_sampler_R.R")
source("src/gjam_gibbs_sampler_Rcpp.R")

# Initialization of function parameters
alpha0<-1 #Dirichlet mass parameter
ndraws=3000 #number of iterations
burnin=500 #number of discarded iterations
S<-150 #number of species
N_stick=min(150,S) #level of truncation of the Dirichlet process
r=4 #number of latent factors
n_sites=100 #number of locations


# This function generates matrix x of covariates and matrix Y of absence/presence
simulation_fun<-function(Sp=S,nsamples=n_sites, r=4, K_t=3){
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
  Ybin<- ifelse(Y_cont>0,1,0)
  Ycont<-Y_cont
  give_back<-list(A_true = Lambda, B_true=B, Xdesign=X, mu_true=L, R_true=R_true, V=Ycont, Y=Ybin)
  
  return(give_back)
}

data<-simulation_fun()
x<-data$Xdesign #x matrix of covariates
Y<-data$Y #Binary matrix Y

#Call of the function (R or Rcpp)
#result<-gjam_gibbs_sampler_Rcpp(alpha0,ndraws,burnin,N_stick,r,S, n_sites,x,Y)
result<-gjam_gibbs_sampler_Rcpp(alpha0,ndraws,burnin,N_stick,r,S, n_sites,x,Y)

# 
# 
# #Analysis of output
chain<-as.vector(result$chain)
chain<-as.numeric(chain)
chain<-as.data.frame(chain)


# # #h<-as.mcmc.list(h1,h2,h3)

# #Intervalli di confidenza per A
mat_bin=matrix(0,nrow=S,ncol=r)
check_CR(result$A_inf, result$A_sup, data$A_true, mat_bin)
