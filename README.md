# PRISM

**Partial Rho Identification through Scale Modeling: Robust Covariance Analysis**

**Version:** 0.1.0

**Author:** Maxwell Konnaris

**Maintainer:** Maxwell Konnaris <Maxwellkonnaris@gmail.com>

## Description

PRISM is an R package for running simulations and robust covariance analysis on count data using a bootstrap scale estimation approach to solve partially identified systems. This package is particularly useful for researchers working with high-dimensional data who need robust methods to estimate covariances.

## Installation

You can install the development version of PRISM from GitHub:

```r
# install.packages("devtools")
devtools::install_github("maxwellkonnaris/PRISM")
```

## Dependencies
PRISM imports the following R packages: \n
<br>
tidyverse \n
driver (downloaded through devtools::install_github("jsilve24/driver") \n
rBeta2009 \n
pbapply \n
progressr \n
foreach \n
doSNOW \n 
parallel \n
stats \n 
MCMCpack \n
ggplot2 \n
gridExtra \n
grid \n
Ensure you have these packages installed before using PRISM.

## Usage
Run Analysis on All Pairwise Taxa \n
The primary function of the PRISM package is estimate_covariance, which runs a bootstrapped analysis on the input data matrix Y using scale models which account for uncertainty due to a partially identified system where the dataset does not include information about scale. \n
See Scale Reliant Inference for more information: 

```r
results = estimate_covariance(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000, variance = FALSE, upperx = 1.0)
forest_plot(results$final_results, save="png", filename="example_dataset")
sigmaplot(results$all_inner_results, save="png", filename="example_dataset")
```

## Parameters
Y: A matrix of data with observations in columns and variables in rows. \n
alpha: A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros. \n
rhobound: A numeric value specifying the bound for the rho1 and rho2 parameters. Defaults to 0.8. \n
S: An integer specifying the number of bootstrap samples. Defaults to 1000. \n
variance: A boolean indicating whether to include diagonal pairs. Defaults to FALSE. \n
upperx: A numeric value specifying the upper bound for the x parameter. Defaults to 1. \n

## Return Value
The function returns a list containing two data frames:

final_results: A data frame with the final results of the analysis, including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances. \n
all_inner_results: A data frame with detailed results from the inner loop of the analysis for each bootstrap sample.

## Functions
estimate_covariance(): Main function for estimating the covariance of all pairwise rows \n
forest_plot(): Plotting the range and confidence interval based on the final_results from estimate_covariance function \n
sigmaplot(): Plotting the relationship between prior parameters specified to minimize or maximize the covariance objective function

## License
PRISM is licensed under the GPL-3 license. See the LICENSE file for more details.

## Bugs
Bug reports are welcome, please direct any issues to 

## Contact
For any questions or issues, please contact the maintainer: Maxwell Konnaris Maxwellkonnaris @ gmail . com

