needed_packages <- c("MASS","coda","ggmcmc","extrafont","mgcv","mvtnorm", "matlib", "devtools", "MCMCpack", "gjam", "invgamma", "MixMatrix", "tictoc", "corpcor") 
new_packages <- needed_packages[!(needed_packages %in% installed.packages()[, "Package"])] 
if (length(new_packages)) install.packages(new_packages) 
lapply(needed_packages, require, character.only = TRUE) 