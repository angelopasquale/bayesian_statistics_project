#Librerie
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())
# Clear console
cat("\014") 
#function that generate a dataset by sampling from the model with given parameters
library("mvtnorm")
library("MASS")
library("corpcor")
library("MCMCpack")
library("coda")
library("ggmcmc")
library("extrafont")
# Install **TTF** Latin Modern Roman fonts from www.fontsquirrel.com/fonts/latin-modern-roman
# Import the newly installed LModern fonts, change the pattern according to the 
# filename of the lmodern ttf files in your fonts folder
#font_import(pattern = "lmodern*")
#loadfonts(device = "win")
#par(family = "LM Roman 10")

Rcpp::sourceCpp('cppFns.cpp')
#Inizializzazione variabili
for (i in seq(1,8)) {
  source(paste("STEP",i,".R", sep = ""))
}
source("compute_GD_prior.R")
source("STEP3RCPP.R")
source("tnorm.R")
source("STEP1RCPP.R")
source("STEP2RCPP.R")
source("STEP4RCPP.R")
source("STEP5RCPP.R")
source("STEP7RCPP.R")
source("STEP8RCPP.R")
source("STEP8TaylorRCPP.R")
source("gjamHfunctions.R")
source("check_IC.R")
source("prova_step.R")

#alpha0=c(1e-3,1e-2,1e-1,1,1e1,1e2,1e3)
alpha0<-1e2
niter=5
burnin=0
N_stick=13
r=4
S<-5 #number of species
n_sites=180
#q<-list()
#for(i in 1:length(alpha0)){
#nn<-prova_step(alpha0[i],niter)
#q[[i]]<-nn
#}
#l<-list()
#for(i in 1:length(alpha0)){
#  l[[i]]<-sum(q[[i]]==1)/length(q[[i]])
#}
#l
h<-prova_step(alpha0,niter,burnin,N_stick,r,S)
#h2<-prova_step(alpha0,niter,N_stick,r,S)
#h3<-prova_step(alpha0,niter,N_stick,r,S)
#h<-list(h1,h2,h3)
#h<-as.mcmc(h)
#h2<-as.mcmc(h2)
#h3<-as.mcmc(h3)
#h<-as.mcmc.list(h1,h2,h3)
#h<-ggs(h)
#GRAFICI DELLA CATENA

#Traceplot
#x11()
#ggs_traceplot(h) + 
#  geom_area(colour="blue") + xlab('Iterations') + ylab('A[2,2] ') + theme(text = element_text(size=15, family="LM Roman 10")) +
 # labs(fill = 'Technologies')+ ggtitle('Traceplot') 
#Running Mean
#x11()
#ggs_running(h)
mat_bin=matrix(0,nrow=S,ncol=r)
check_IC(h$A_inf, h$A_sup, h$A_true, mat_bin)

#Autocorrelation
#x11()
#ggs_autocorrelation(h1)

