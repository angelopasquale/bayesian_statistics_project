# Joint Species Distribution Models - Gibbs sampling implementation in RCpp

This project aims at implementing a Gibbs Sampler for the GJAM model presented by [Taylor and Rodriguez](https://projecteuclid.org/euclid.ba/1478073617). Here below the complete model:

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites and installing

### Installing

The script **install.R** provides automatical installation of required **R** packages as well as **Rcpp** interface. This is done simply by the following check: 

```
needed_packages  <- c("MASS","coda","ggmcmc","extrafont","mgcv","mvtnorm", "matlib", "devtools", "MCMCpack", "gjam", "invgamma","MixMatrix", "tictoc", "corpcor")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])] if (length(new_packages))
install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)}
```

## Running the tests

Run the script **main.R** for a testcase.  

### Break down into end to end tests

Data and model true parameters are simulated through **simulation_fun** contained in **main.R**. Here, dimensions of the model need to be fixed:
* **Sp** = number of species
* **nsamples** = number of sites 
* **r** = number of latent factor 
* **K_t** = number of clusters in the matrix **A_true**.

Then, a call to the Gibbs sampler function **gjam_gibbs_sampler** is performed and posterior draws for the model are obtained. 

Afterwards, the **check_CR** function is called in order to check credible regions for the chains with respect to the true (benchmark) values.

Moreover, the analysis of convergence of the chains is performed through traceplots, autocorrelation and running mean.

## Built With

* [Armadillo](http://arma.sourceforge.net) - for linear algebra in C+++
* [Eigen](https://eigen.tuxfamily.org/dox/) - for linear algebra in C+++
* [gjam](https://cran.r-project.org/web/packages/gjam/index.html) - for GJAM modelling
* [ggmcmc](https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html) - for analysis of Gibbs chains

## Authors

* **Angelo Pasquale** - *Politecnico di Milano - Ecole Centrale de Nantes - MSc in Computational Science and Engineering* -
* **Matteo Contini** - *Politecnico di Milano - Ecole Centrale de Lyon - MSc in Applied Statistics* -
* **Lorenzo Fiorello** - *Politecnico di Milano - MSc in Applied Statistics* -

## Acknowledgments

* Bayesian Statistics course professor and assistants
* Project tutors Doct. Riccardo Corradin and Poggiato Giovanni (PhD at Laboratoire d’Ecologie Alpine - Inria Grenoble)


