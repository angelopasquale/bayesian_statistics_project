# Joint Species Distribution Models - Gibbs sampling implementation in RCpp

This project aims at implementing a Gibbs Sampler for the GJAM model presented by [Taylor and Rodriguez](https://projecteuclid.org/euclid.ba/1478073617). Here below the complete model:


<p align="center"><img alt="\begin{equation}&#10;\begin{split}&#10;    \boldsymbol{V}_i|\boldsymbol{k},\boldsymbol{Z},\boldsymbol{w}_i,\boldsymbol{B},\sigma_{\epsilon}^2 &amp; \sim \mathcal{N}_S(\boldsymbol{B}\boldsymbol{X}_i+\boldsymbol{Q}(k)\boldsymbol{Z}\boldsymbol{w}_i,\sigma_{\epsilon}^2\boldsymbol{I}_S), \quad \text{for } i=1,\dots,n \\&#10;    [\boldsymbol{B}, \sigma_{\epsilon}^2] &amp; \propto \frac{1}{\sigma_{\epsilon}^2} \\&#10;    \boldsymbol{w}_i &amp; \sim \mathcal{N}_r(\mathbf{0},\boldsymbol{I}_r) \\&#10;    k_l|\mathbf{p} &amp; \overset{ \text{iid} }{ \sim } \sum_{j=1}^N p_j\delta_j(k_l), \quad \text{for } l=1,\dots S \\&#10;    \boldsymbol{Z}_j|\boldsymbol{D_z} &amp; \overset{ \text{iid} }{ \sim } \mathcal{N}_r(\mathbf{0},\boldsymbol{D_z}), \quad \text{for } j=1,\dots N \\&#10;    \boldsymbol{p} &amp; \sim \mathcal{GD}_N(a_{\alpha},b_{\alpha})\\&#10;    \boldsymbol{D_z} &amp; \sim \mathcal{IW}(2+r-1,4\text{diag}(\frac{1}{\eta_1},\dots, \frac{1}{\eta_r})) \\&#10;    \eta_h &amp; \sim \mathcal{IG}(\frac{1}{2},\frac{1}{10^4}), \quad \text{for } h=1,\dots, r&#10;\end{split}&#10;\end{equation}" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/a95e823763e7c1b1f5eaec470207f275.svg" align="middle" width="595.33009635pt" height="288.38918625pt"/></p>

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
* ***WIP_code*** contains work-in-progress code.

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

* [Armadillo](http://arma.sourceforge.net) - for linear algebra in C++
* [Eigen](https://eigen.tuxfamily.org/dox/) - for linear algebra in C++
* [gjam](https://cran.r-project.org/web/packages/gjam/index.html) - for GJAM modelling
* [ggmcmc](https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html) - for analysis of Gibbs chains

## Authors

* **Angelo Pasquale** - *Politecnico di Milano - Ecole Centrale de Nantes - MSc in Computational Science and Engineering* -
* **Matteo Contini** - *Politecnico di Milano - Ecole Centrale de Lyon - MSc in Applied Statistics* -
* **Lorenzo Fiorello** - *Politecnico di Milano - MSc in Applied Statistics* -

## Acknowledgments

* Bayesian Statistics course Professor Alessandra Guglielmi
* Project supervisors Dr.s Riccardo Corradin and Poggiato Giovanni (PhD at Laboratoire d’Ecologie Alpine - Inria Grenoble)
* Project tutor Dr. Mario Beraha
