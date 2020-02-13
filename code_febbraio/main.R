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
#library("extrafont")
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
source(".tnorm.R")
source("STEP1RCPP.R")
source("STEP2RCPP.R")
source("STEP4RCPP.R")
source("STEP5RCPP.R")
source("STEP7RCPP.R")
source("STEP8RCPP.R")
source("gjamHfunctions.R")
source("check_IC.R")
source("prova_step.R")
#alpha0=c(1e-3,1e-2,1e-1,1,1e1,1e2,1e3)
alpha0<-1e2
niter=200
N_stick=13
r=4
S<-5 #number of species
n_sites=20
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
h<-prova_step(alpha0,niter,N_stick,r,S)