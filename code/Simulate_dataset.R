
#function that generate a dataset by sampling from the model with given parameters

# Sp is the number of species
# nsamples is the number of sites
# r is the number of latent factors
# K_t is the number of unique lines of A
library(MASS)
library(repmis)
library(rlist)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(NLRoot)
library(Rcpp)

needed_packages <- c("mvtnorm", "matlib", "devtools", "MCMCpack", "invgamma", "MixMatrix", "tictoc", "corpcor") 
new_packages <- needed_packages[!(needed_packages %in% installed.packages()[, "Package"])] 
if (length(new_packages)) install.packages(new_packages) 
lapply(needed_packages, require, character.only = TRUE) 

library(tictoc)
library(corpcor) # for make positive definite

library("mvtnorm")
library("matlib")
library("devtools")
#install_github("cran/MCMCpack")#, for inverse wishart, but just for the moment
library("MCMCpack")
library("invgamma")
library("MixMatrix")

setwd("/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/implementation/gjam/")
Rcpp::sourceCpp('src/cppFns.cpp') #in gjam sources

simulation_fun<-function(Sp=5,nsamples=500, r=5, K_t=40){
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
    A[i,]<-mvrnorm(n = 1, rep(20,r), Sigma=3*diag(r)) # We sample the unique values of A
  }
  idx<-sample((1:K_t),S,replace=T) #idx represents the labels
  Lambda<-A[idx,] #We expand A using the labels to obtain Lambda (what we used to call A)
  Sigma_true<-Lambda%*%t(Lambda)+0.1*diag(S) #We obtain Sigma. Here sigma_epsilon^2 is 0.1
  
  
  #We obtain the correlation matrix R
  R_true=cov2cor(Sigma_true) 
  
  #We sample V from the model
  Y_cont<-L+rmvnormRcpp(n = n, mu=rep(0,S), sigma=R_true)
  # We obtain our binary dataset. Secondo me vi conviene prima testare il modello con Y_cont e poi su Y
  Y<- ifelse(Y_cont>0,1,0)
  
  give_back<-list(B_true=B, Xdesign=X, mu_true=L, R_true=R_true, V=Y_cont, Y=Y)
  return(give_back)
}

data<-simulation_fun()

source(file = "/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/Bayesian_Statistics_Project/code/GJAM_Gibbs_Sampler.R")
return_list <- GJAM_Gibbs_Sampler(data$Xdesign, data$Y, 5, 40, 1e0, 1e3, 1)

Q <- apply(simplify2array(return_list$B), 1:2, quantile, prob = c(0.05, 0.95))
M <- apply(simplify2array(return_list$B), 1:2, mean)
