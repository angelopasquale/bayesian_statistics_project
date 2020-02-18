########### Simulation from rewritten .gjamReduct ###########
source(file = "gjamHfunctions.R")
library(gjam)

f <- gjamSimData(n = 30, S = 6, Q = 2, typeNames = 'PA')

#source(file = "/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/Bayesian_Statistics_Project/code/gjam.R")

rl   <- list(r = 4, N = 4)
ml   <- list(ng = 2500, burnin = 500, typeNames = 'PA', reductList = rl)

out  <- gjamReduct(f$formula, xdata = f$xdata, ydata = f$ydata, modelList = ml)

out$parameters$betaMu         # S by M coefficient matrix alpha
out$parameters$betaStandXmu   # S by M standardized for X
out$parameters$betaStandXWmu  # (S-F) by M standardized for W/X, centered factors

out$parameters$betaTable        # SM by stats posterior summary
out$parameters$betaStandXtable  # SM by stats posterior summary
out$parameters$betaStandXWtable # (S-F)M by stats posterior summary

out$parameters$sigMu         # S by S covariance matrix omega
out$parameters$sigSe         # S by S covariance std errors

pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T)
x11()
g <-gjamPlot(output = out, plotPars = pl)
