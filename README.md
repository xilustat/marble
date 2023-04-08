
<!-- README.md is generated from README.Rmd. Please edit that file -->

# marble

> Robust Marginal Bayesian Variable Selection 

Recently, multiple marginal variable selection methods have been developed and shown to be effective in Gene-Environment interactions studies. We propose a novel marginal Bayesian variable selection method for Gene-Environment interactions studies. In particular, our marginal Bayesian method is robust to data contamination and outliers in the outcome variables. With the incorporation of spike-and-slab priors, we have implemented the Gibbs sampler based on Markov Chain Monte Carlo (MCMC). The core algorithms of the package have been developed in C++.

## How to install

 - To install from github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("xilustat/marble")

## Examples

#### Example.1 (default method: robust sparse marginal selection)

    library(marble)
    data(dat)
    
    max.steps=10000
    fit=marble(X, Y, E, clin, max.steps=max.steps)
    selected=GxESelection(fit,sparse=TRUE)
    selected
    

#### Example.2 (alternative: non-robust sparse group selection)

    fit=marble(X, Y, E, clin, max.steps=max.steps, sparse=FALSE)
    selected=GxESelection(fit,sparse=FALSE)
    selected


## Methods

This package provides implementation for methods proposed in

  - Lu, X., Fan, K., Ren, J., and Wu, C. (2021). Identifying Gene–Environment Interactions With Robust Marginal Bayesian Variable Selection. Frontiers in Genetics. 12:667074.
