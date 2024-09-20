#' @title Package Startup Script
#' @description Installs necessary GitHub and Bioconductor packages
#' @details
#' This function installs the following packages if they are not already installed:
#' - GitHub packages: driver, propr, SpiecEasi
#' - Bioconductor package: phyloseq, banocc
#' The function is executed upon loading the package.
#' @param libname The library name.
#' @param pkgname The package name.
#' @importFrom devtools install_github
#' @importFrom BiocManager install
.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  } 
  if (!requireNamespace("driver", quietly = TRUE)) {
    devtools::install_github("jsilve24/driver")
  }
  if (!requireNamespace("propr", quietly = TRUE)) {
    devtools::install_github("tpq/propr")
  }
  if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
    devtools::install_github("zdk123/SpiecEasi")
  }
  if (!requireNamespace("banocc", quietly = TRUE)) {
    BiocManager::install("banocc", force=TRUE)
  }
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    BiocManager::install("phyloseq", force=TRUE)
  }
}
