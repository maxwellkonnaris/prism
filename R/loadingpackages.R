#' @title Package Startup Script
#' @description Installs necessary GitHub and Bioconductor packages if not installed.
#' @param libname The library name.
#' @param pkgname The package name.
.onLoad <- function(libname, pkgname) {  
  # Function to install GitHub packages
  install_github_package <- function(pkg, repo) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing '", pkg, "' package from GitHub...", sep = ""))
      tryCatch({
        devtools::install_github(repo)
      }, error = function(e) {
        message(paste("Failed to install '", pkg, "' from GitHub. Error: ", e$message, sep = ""))
      })
    }
  }

  # Function to install Bioconductor packages
  install_bioconductor_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing '", pkg, "' package from Bioconductor...", sep = ""))
      tryCatch({
        BiocManager::install(pkg, dependencies = TRUE)
      }, error = function(e) {
        message(paste("Failed to install '", pkg, "' from Bioconductor. Error: ", e$message, sep = ""))
      })
    }
  }

  # Install GitHub packages
  install_github_package("driver", "jsilve24/driver")
  install_github_package("propr", "tpq/propr")
  install_github_package("SpiecEasi", "zdk123/SpiecEasi")
  
  # Install Bioconductor packages
  install_bioconductor_package("banocc")
  install_bioconductor_package("phyloseq")
}
