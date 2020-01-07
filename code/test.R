rm ( list = ls() )
setwd("/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/Bayesian_Statistics_Project/code")
library("mvtnorm")
library("matlib")
library("devtools")
#install_github("cran/MCMCpack")#, for inverse wishart, but just for the moment
library("MCMCpack")
library("invgamma")
library("MixMatrix")
library("gdirmn")
library(gjam)
library(tictoc)

n = 15
S = 20
n_cov = 3

f <- gjamSimData(n = n, S = S, Q = n_cov, typeNames = 'PA')
# The object `f` includes elements needed to analyze the simulated data set.  
# `f$typeNames` is a length-$S$ `character vector`. The `formula` follows 
# standard R syntax. It does not start with `y ~`, because gjam is multivariate.
# The multivariate response is supplied as a $n \times S$ `matrix` or
# `data.frame ydata`.
summary(f)
#ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
#out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
#summary(out)

x = as.matrix(f$xdata) # matrix of measured (environmental) covariates (n * k)
Y = as.matrix(f$ydata) # matrix of n_species presence/absence data (in a continuous framework) (n * S)

# Number of iterations
posteriorDraws = 1e4
burnInIterations = 1

# Number of latent factors
r = 5
N_stick = 20
alpha_0 = 1e2

# burnInIterations: butto via i valori fino al burnInIterations per mantenere
# solo la stazionare
# A, sigmaeps2, B, k

source(file = "GJAM_Gibbs_Sampler.R")
return_list <- GJAM_Gibbs_Sampler(x, Y, r, N_stick, alpha_0, posteriorDraws, burnInIterations)

# voglio la distribuzione di un nuovo Y -> Y_new sim N(mu_new, sigma_new).
# dalla catena ho la distribuzione a posteriori di mu_new e di sigma_new.
# per ogni elemento della lista, per ogni B, mi calcolo mu e sigma (media e varianza
# a posteriori). Se aggiungo un nuovo sito, avrò un nuovo elemento e voglio la predittiva
# 

# 1. ANALISI DELLA CONVERGENZA CATENE 
# plot del numero di clusters in funzione di alpha
number_of_clusters=c()
count=1
alpha=c(1e-3, 1e-2, 1e-1, 1e-0, 1e1, 20, 30, 40, 50, 60, 70, 80, 90, 100)
for(i in alpha){
  alpha_0=abs(i)
  return_list <- GJAM_Gibbs_Sampler(x, Y, r, N_stick, alpha_0, posteriorDraws, burnInIterations)
  number_of_clusters[count]=length(unique((return_list$k[[1]])))
  count=count+1
}
par(family = "LM Roman 10")

plot(alpha, number_of_clusters)



