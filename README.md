# PRISM

**Partial Rho Identification through Scale Modeling: Robust Covariance and Correlation Analysis**

**Version:** 0.1.0

**Author:** Maxwell Konnaris, Michelle Nixon, and Justin Silverman

**Maintainer:** Maxwell Konnaris

**Affiliation(s):** Pennsylvania State University

## Description

PRISM is an R package for running simulations and robust covariance analysis on count data using a bootstrap scale estimation approach to solve partially identified systems. This package is particularly useful for researchers working with high-dimensional data who need robust methods to estimate covariances. PRISM utilizes parallel computing for computationally efficient estimation dependent on available CPU cores, the size of matrix, and the number of bootstrap samples specified. <br>
Local 11 core CPUs, 10x286 matrix, and default parameters: Expected 1 minutes 32 seconds. <br>
HPC 47 core CPUs, 10x286 matrix, and 1000 bootstrap samples: Expected 15 seconds. <br>

## Installation

You can install the development version of PRISM from GitHub:

```r
# install.packages("devtools")
devtools::install_github("maxwellkonnaris/PRISM")
```

PRISM is also available via CRAN:

```r
install.packages("PRISM")
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
<br>
The primary function of the PRISM package is estimate_covariance, which runs a bootstrapped analysis on the input data matrix Y using scale models which account for uncertainty due to a partially identified system where the dataset does not include information about scale. <br>
<br>
See Scale Reliant Inference for more information: https://arxiv.org/abs/2201.03616 

```r
results = estimate_covariance(Y, alpha = rep(0, nrow(Y)), rhobound = 0.9, S = 1000, upperscalevariance = 1.0)
forest_plot(results$final_results, save="png", filename="example_dataset")
sigmaplot(results$all_inner_results, save="png", filename="example_dataset")
calculate_bootstrap_summary(results$all_inner_results, group_col = "comparison", exclude_cols = c("d1", "d2"), save_as_csv = TRUE, csv_path = "my_summary_stats.csv")
```

## Parameters
- **Y**: A matrix of data with observations in columns and variables in rows. <br>
- **alpha**: A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros. <br>
- **rhobound**: A numeric value specifying the bound for the rho1 and rho2 parameters. Defaults to 0.9. <br>
- **S**: An integer specifying the number of bootstrap samples. Defaults to 1000. <br>
- **upperscalevariance**: A numeric value specifying the upper bound for the scalevariance parameter. Defaults to 1.0. <br>

## Return Value
The function returns a list containing two data frames:<br>
<br>
- **final_results**: A data frame with the final results of the analysis, including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances. <br>
- **all_inner_results**: A data frame with detailed results from the inner loop of the analysis for each bootstrap sample.

## Main Functions
- **estimate_covariance()**: Main function for estimating the covariance of all pairwise rows <br>
- **forest_plot()**: Plotting the range and confidence interval based on the final_results from estimate_covariance function <br>
- **sigmaplot()**: Plotting the relationship between prior parameters specified to minimize or maximize the covariance objective function
- **calculate_bootstrap_summary()**: Generate a table of summary stats for a csv, more specifically applied to the results$all_inner_results if youre interested in the summary stats from the bootstrap performed

## License
PRISM is licensed under the GPL-3 license. See the LICENSE file for more details.

## Contact: Bugs, Questions, or Issues
Bug reports, questions, or issues are welcome, please direct any issues to https://github.com/maxwellkonnaris/prism/issues


