########### Simulation from rewritten .gjamReduct ###########
source(file = "gjamHfunctions.R")
library(gjam)

f <- gjamSimData(n = 100, S = 5, Q = 2, typeNames = 'PA')

rl   <- list(r = 4, N = 13)
ml   <- list(ng = 1000, burnin = 100, typeNames = 'PA', reductList = rl)

#Here we call the original gjam function
tic("GJAM")
out  <- .gjam(f$formula, xdata = f$xdata, ydata = f$ydata, modelList = ml)
toc()
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
