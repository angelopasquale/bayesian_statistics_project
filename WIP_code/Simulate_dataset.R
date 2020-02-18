# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())
# Clear console
cat("\014") 
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

#setwd("/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/implementation/gjam/")
setwd("C:/Users/loren/OneDrive/Desktop/bayesian_project/gjam_2.2.7/gjam")
#setwd("~/Desktop/Polimi5anno/Bayesiana/Progetto_Bayes/gjam_2.2.7")
Rcpp::sourceCpp("src/cppFns.cpp") #in gjam sources

simulation_fun<-function(Sp=15,nsamples=400, r=6, K_t=4){
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
  Y_cont<-L+rmvnormRcpp(n = n, mu=rep(0,S), sigma=R_true)
  # We obtain our binary dataset. Secondo me vi conviene prima testare il modello con Y_cont e poi su Y
  Y<- ifelse(Y_cont>0,1,0)
  
  give_back<-list(A_true = Lambda, B_true=B, Xdesign=X, mu_true=L, R_true=R_true, V=Y_cont, Y=Y)
  return(give_back)
}

data<-simulation_fun()
source(file = "C:/Users/loren/OneDrive/Desktop/bayesian_project/bayesian_statistics_project/code/GJAM_Gibbs_Sampler_R_last.R")
#source(file = "/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/Bayesian_Statistics_Project/code/GJAM_Gibbs_Sampler_R_last.R")
#source(file = "/Users/lombardata/Desktop/Polimi5anno/Bayesiana/Progetto_Bayes/bayesian_statistics_project/code/GJAM_Gibbs_sampler_R_last.R")

return_list <- GJAM_Gibbs_Sampler_R_last(data$Xdesign, data$Y, 6, 80, 1e1, 1e2, 1)

B = return_list$B

apply(simplify2array(B),1:2,quantile,0.95)
apply(simplify2array(B),1:2,quantile,0.05)
apply(simplify2array(B),1:2, mean)

### DIAGNOSTIC IN CODA
#X_mc <- mcmc(data = X, start = burnin + 1, end = niter, thin = thin)
#plot(X_mc)
#summary(X_mc)

# Autocovariances
#acfplot(X_mc, lag.max = 30)

# (0.025; 0.5, 0.975) quantiles 
#cumuplot(X_mc)

# Effective sample size:
#effectiveSize(X_mc)
#dim(X_mc)[1]

Q <- apply(simplify2array(return_list$B), 1:2, quantile, prob = c(0.05, 0.95))
M <- apply(simplify2array(return_list$B), 1:2, mean)


########### Simulation from rewritten .gjamReduct ###########

library(gjam)

f <- gjamSimData(n = 100, S = 5, Q = 4, typeNames = 'PA')

source(file = "/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/Bayesian_Statistics_Project/code/gjam.R")

rl   <- list(r = 3, N = 4)
ml   <- list(ng = 3000, burnin = 500, typeNames = 'PA', reductList = rl)

out  <- .gjamReduct(f$formula, xdata = f$xdata, ydata = f$ydata, modelList = ml)

out$parameters$betaMu         # S by M coefficient matrix alpha
out$parameters$betaStandXmu   # S by M standardized for X
out$parameters$betaStandXWmu  # (S-F) by M standardized for W/X, centered factors

out$parameters$betaTable        # SM by stats posterior summary
out$parameters$betaStandXtable  # SM by stats posterior summary
out$parameters$betaStandXWtable # (S-F)M by stats posterior summary

out$parameters$sigMu         # S by S covariance matrix omega
out$parameters$sigSe         # S by S covariance std errors

pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T)
gjamPlot(output = out, plotPars = pl)
