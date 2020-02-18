# Joint Species Distribution Models - Gibbs sampling implementation in RCpp

This project aims at implementing a Gibbs Sampler for the here below GJAM model, presented by [Taylor and Rodriguez](https://projecteuclid.org/euclid.ba/1478073617):

![alt text](https://github.com/angelopasquale/bayesian_statistics_project/blob/master/images/model-1.png)

For a detailed mathematical description of the model, please refer to our report, contained in this Github repository.

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
### Code and folders structure
* ***code*** contains:
  * ***main.R*** calls the main functions:
    * **simulation_fun** : data and model true parameters generation;
    * **gjam_gibbs_sampler** : Gibbs sampler using data previously simulated;
    * **check_CR** : confidence regions for the chains;
    * **chain_convergence** : analysis of Gibbs chain.
  * ***simulation_gjam.R*** is used for a preliminary analysis based on the library [**gjam**](https://cran.r-project.org/web/packages/gjam/index.html);
  * ***src*** folder contains source files:
    * ***cpp/cppFns.cpp*** contains **c++** functions;
    * ***chain_convergence.R***;
    * ***check_CR***;
    * ***compute_GD_prior***, which computes the Generalized Dirichlet prior;
    * ***gjam_gibbs_sampler_Rcpp.R***, contains our Gibbs sampler in Rcpp;
    * ***gjam_gibbs_sampler_R.R***, contains our Gibbs sampler in Rcpp;
    * ***tnorm.R***, function for sampling from a truncated normal;
    * ***gjam.R***, which contains an extrapolation from **gjam** library source code (needed from ***simulation_gjam.R***);
    * ***gjamHfunctions.R***, which contains **gjam** library source functions (needed from ***simulation_gjam.R***);
  * ***steps_R***, which contains all steps of the Gibbs sampler in **R** as separated functions;
  * ***steps_Rcpp***, which contains all steps of the Gibbs sampler in **Rcpp** as separated functions;
* ***references*** folder contains main bibliography exploited in our work;
* ***WIP_code*** contains work-in-progress code;
* ***report.pdf*** contains a report of the project.

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

### Some results

Considering for instance an element of matrix **A**, we get the chain:
* autocorrelation
![alt text](https://github.com/angelopasquale/bayesian_statistics_project/blob/master/images/autocorr10000.png){:height="50%" width="50%"}
* values
![alt text](https://github.com/angelopasquale/bayesian_statistics_project/blob/master/images/freq10000.png){:height="50%" width="50%"}
* running mean
![alt text](https://github.com/angelopasquale/bayesian_statistics_project/blob/master/images/rm10000iter.png){:height="50%" width="50%"}
* traceplot
![alt text](https://github.com/angelopasquale/bayesian_statistics_project/blob/master/images/ts.png){:height="50%" width="50%"}

## Built With

* [Armadillo](http://arma.sourceforge.net) - for linear algebra in C++
* [Eigen](https://eigen.tuxfamily.org/dox/) - for linear algebra in C++
* [gjam](https://cran.r-project.org/web/packages/gjam/index.html) - for GJAM modelling
* [ggmcmc](https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html) - for analysis of Gibbs chains

## Authors

* **Angelo Pasquale** - *Politecnico di Milano - Ecole Centrale de Nantes - MSc in Computational Science and Engineering*
* **Matteo Contini** - *Politecnico di Milano - Ecole Centrale de Lyon - MSc in Applied Statistics*
* **Lorenzo Fiorello** - *Politecnico di Milano - MSc in Applied Statistics*

## Acknowledgments

* Bayesian Statistics course Professor Alessandra Guglielmi
* Project supervisors Dr.s Riccardo Corradin and Poggiato Giovanni (PhD at Laboratoire d’Ecologie Alpine - Inria Grenoble)
* Project tutor Dr. Mario Beraha


