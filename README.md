# PRISM 

**Partial Rho Identification through Scale Modeling: Robust Covariance and Correlation Analysis**

**Version:** 0.1.0

**Author:** Maxwell Konnaris, Michelle Nixon, and Justin Silverman

**Maintainer:** Maxwell Konnaris

**Affiliation(s):** Pennsylvania State University

## Description

PRISM is an R package for running simulations and robust covariance analysis on count data using a bootstrap scale estimation approach to solve partially identified systems. This package is particularly useful for researchers working with high-dimensional data who need robust methods to estimate covariances. PRISM utilizes parallel computing for computationally efficient estimation dependent on available CPU cores, the size of matrix, and the number of bootstrap samples specified. <br>
<br>
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
PRISM imports the following R packages:

- **driver**: For managing simulations and workflows. (Separate install below and not auto imported)
- **tidyverse**: A collection of R packages for data manipulation and visualization (includes `ggplot2`, `dplyr`, `tidyr`, etc.).
- **rBeta2009**: A package for working with beta distributions.
- **pbapply**: An alternative to `apply` functions with built-in progress bars.
- **progress**: A package for displaying progress bars in R scripts.
- **progressr**: Provides an API to track the progress of computations, used for longer running tasks.
- **foreach**: For loop constructs that allow parallel and distributed execution.
- **doSNOW**: A parallel backend for the `foreach` package, particularly for Windows users.
- **parallel**: The core R package that supports parallel computing.
- **stats**: A core R package for statistical functions.
- **MCMCpack**: For performing Markov Chain Monte Carlo (MCMC) simulations.
- **ggplot2**: Part of the tidyverse, used for advanced data visualization.
- **profvis**: A graphical profiler for R to analyze performance bottlenecks.
- **gridExtra**: A package for arranging multiple grid-based plots (e.g., ggplot2) on a single page.
- **grid**: Core R package that provides low-level functions for creating and manipulating graphical objects.
- **compositions**: Provides tools for working with compositional data, including CLR and ILR transformations.
- **GGally**: Extends `ggplot2` with additional plotting functionality, especially for data diagnostics and visualizations.
- **filelock**: Used to safely manage file locking when running parallel computations.
- **nloptr**: A package for nonlinear optimization, used in optimization problems within the package.
- **reshape2**: Used for data reshaping, particularly in preparing data for analysis.
- **plotly**: For creating interactive visualizations.
- **htmlwidgets**: For creating interactive web visualizations within R.
- **viridis**: A color palette for `ggplot2` plots, optimized for perceptual uniformity and colorblind-friendliness.
- **alphashape3d**: For 3D alpha shapes, used in visualization and analysis of 3D data.
- **phyloseq**: For microbiome data analysis, used for Banocc covariance estimation.
- **SpiecEasi**: For sparse inverse covariance estimation for ecological association inference. 
- **CCLasso**: For sparse covariance estimation using compositional data using Lasso.
- **SparCC**: For sparse covariance estimation using compositional data.
- **propr**: For analyzing proportionality relationships in compositional data.
- **BAnOCC**: Bayesian Analaysis Of Compositional Covariance.

<br>

Ensure you have these packages installed before using PRISM.

## Usage
Run Analysis on All Pairwise Combinations <br>
<br>
The primary function of the PRISM package is estimate_covariance, which runs a bootstrapped analysis on the input data matrix Y using scale models which account for uncertainty due to a partially identified system where the dataset does not include information about scale. <br>
<br>
See Scale Reliant Inference for more information: https://arxiv.org/abs/2201.03616 

```r
results = estimate_covariance(Y, alpha = rep(0, nrow(Y)), lowerrhobound = rep(0.9, nrow(Y)), upperrhobound = rep(0.9, nrow(Y)), S = 1000, lowerscalestdev = 0.45, upperscalestdev = 0.6, outputdirectory='/tables/')
forest_plot(results$final_results, save="png", filename="example_dataset")
sigmaplot(results$all_inner_results, save="png", filename="example_dataset")
calculate_bootstrap_summary(results$all_inner_results, group_col = "comparison", exclude_cols = c("d1", "d2"), save_as_csv = TRUE, csv_path = "my_summary_stats.csv")
```

## Parameters for Main Estimation Functions ( estimate_covariance() )
- **Y**: A matrix of observed counts, where rows represent variables (e.g., taxa) and columns represent observations (samples). The matrix should have at least two rows and two columns. <br>
- **alpha**: A numeric vector of Dirichlet priors with the same length as the number of rows in \code{Y}. If \code{alpha} is a scalar, it will be replicated for each row. Defaults to \code{0.5}. <br>
- **lowerrhobound**: A numeric vector of lower bounds for correlation parameters (\eqn{\rho}). Each element specifies the lower bound for a row of \code{Y}. If a scalar is provided, it will be replicated for each row. Defaults to \code{-1.0}. <br>
- **upperrhobound**: A numeric vector of upper bounds for correlation parameters (\eqn{\rho}). Each element specifies the upper bound for a row of \code{Y}. If a scalar is provided, it will be replicated for each row. Defaults to \code{1.0}. <br>
- **S**: An integer specifying the number of bootstrap samples. Larger values reduce Monte Carlo error but increase computation time. Defaults to \code{1000}. <br>
- **lowerscalestdev**: A numeric value specifying the lower bound of the standard deviation of the scale. Defaults to \code{0.49}. <br>
- **upperscalestdev**: A numeric value specifying the upper bound of the standard deviation of the scale. Defaults to \code{0.51}. <br>
- **algorithm**: A character string specifying the optimization algorithm to be used. Can be one of \code{"COBYLA"}, \code{"MMA"}, \code{"AUGLAG_COBYLA"}, \code{"AUGLAG_MMA"}, or \code{"GRID_SEARCH"}. Defaults to \code{"COBYLA"}. <br>
- **outputdirectory**: A character string specifying the directory to save results for grid search. If \code{NULL}, the current working directory is used. Defaults to \code{NULL}.

## Return Value
The function returns a list containing two data frames:<br>
<br>
- **final_results**: A data frame with the final results of the analysis, including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances. <br>
- **all_inner_results**: A data frame with detailed results from the inner loop of the analysis for each bootstrap sample.

## Main Estimation Functions
- **estimate_covariance()**: Main function for estimating the covariance of all pairwise rows using a Multinomial Dirichlet Bootstrap <br>
- **estimate_covariance_MLN()**: Main function for estimating the covariance of all pairwise rows using a Multinomial Logistic Normal Bootstrap <br>
- **estimate_covariance_convergence()**: Main function for estimating the covariance of all pairwise rows using a Multinomial Dirichlet Bootstrap across several different iterations of bootstrap samples <br>

## Main Plotting and Diagnostic Functions
- **forest_plot()**: Plotting the range and confidence interval based on the final_results from estimate_covariance function <br>
- **proportiondontcoverzerobars()**: Plotting the distribution of the proportion of intervals that do not cover zero based on the final_results from estimate_covariance function <br>
- **sigmaplot()**: Plotting the relationship between prior parameters specified to minimize or maximize the covariance objective function <br>
- **plot_bivariate_grid()**: Plotting a bivariate grid to show how each parameter contributes to the optimization of the covariance interval for a given grid search result between two taxa <br>
- **plot_rpars_3d_scatter()**: Plotting a 3D scatter of the grid search optimization of the covariance between two taxa <br>
- **calculate_bootstrap_summary()**: Generate a table of summary stats for a csv, more specifically applied to the results$all_inner_results if youre interested in the summary stats from the bootstrap performed <br>


## License
PRISM is licensed under the GPL-3 license. See the LICENSE file for more details.

## Contact: Bugs, Questions, or Issues
Bug reports, questions, or issues are welcome, please direct any issues to https://github.com/maxwellkonnaris/prism/issues


