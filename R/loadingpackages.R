#' @title Package Startup Script
#' @description Installs necessary GitHub and Bioconductor packages if not installed.
#' @param libname The library name.
#' @param pkgname The package name.
.onLoad <- function(libname, pkgname) {
  # Check and install BiocManager if not available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install GitHub packages if not installed
  if (!requireNamespace("driver", quietly = TRUE)) {
    message("Installing 'driver' package from GitHub...")
    devtools::install_github("jsilve24/driver")
  }
  
  if (!requireNamespace("propr", quietly = TRUE)) {
    message("Installing 'propr' package from GitHub...")
    devtools::install_github("tpq/propr")
  }
  
  if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
    message("Installing 'SpiecEasi' package from GitHub...")
    devtools::install_github("zdk123/SpiecEasi")
  }
  
  # Install Bioconductor packages if not installed
  if (!requireNamespace("banocc", quietly = TRUE)) {
    message("Installing 'banocc' package from Bioconductor...")
    BiocManager::install("banocc")
  }
  
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    message("Installing 'phyloseq' package from Bioconductor...")
    BiocManager::install("phyloseq")
  }
}
