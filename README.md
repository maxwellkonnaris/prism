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
PRISM imports the following R packages: <br>
<br>
- driver (downloaded through devtools::install_github("jsilve24/driver") <br>
- tidyverse<br>
- rBeta2009<br>
- pbapply<br>
- progressr<br>
- foreach<br>
- doSNOW<br>
- parallel<br>
- stats<br>
- MCMCpack<br>
- ggplot2<br>
- profvis<br>
- gridExtra<br>
- grid<br>
<br>
Ensure you have these packages installed before using PRISM.

## Usage
Run Analysis on All Pairwise Combinations <br>
The primary function of the PRISM package is estimate_covariance, which runs a bootstrapped analysis on the input data matrix Y using scale models which account for uncertainty due to a partially identified system where the dataset does not include information about scale. <br>
See Scale Reliant Inference for more information: 

```r
results = estimate_covariance(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000, variance = FALSE, upperx = 1.0)
forest_plot(results$final_results, save="png", filename="example_dataset")
sigmaplot(results$all_inner_results, save="png", filename="example_dataset")
```

## Parameters
**Y**: A matrix of data with observations in columns and variables in rows. <br>
**alpha**: A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros. <br>
**rhobound**: A numeric value specifying the bound for the rho1 and rho2 parameters. Defaults to 0.8. <br>
**S**: An integer specifying the number of bootstrap samples. Defaults to 1000. <br>
**variance**: A boolean indicating whether to include diagonal pairs. Defaults to FALSE. <br>
**upperx**: A numeric value specifying the upper bound for the x parameter. Defaults to 1. <br>

## Return Value
The function returns a list containing two data frames:<br>
<br>
**final_results**: A data frame with the final results of the analysis, including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances. <br>
**all_inner_results**: A data frame with detailed results from the inner loop of the analysis for each bootstrap sample.

## Functions
**estimate_covariance()**: Main function for estimating the covariance of all pairwise rows <br>
**forest_plot()**: Plotting the range and confidence interval based on the final_results from estimate_covariance function <br>
**sigmaplot()**: Plotting the relationship between prior parameters specified to minimize or maximize the covariance objective function

## License
PRISM is licensed under the GPL-3 license. See the LICENSE file for more details.

## Bugs
Bug reports are welcome, please direct any issues to 

## Contact
For any questions or issues, please contact the maintainer: Maxwell Konnaris Maxwellkonnaris @ gmail . com

