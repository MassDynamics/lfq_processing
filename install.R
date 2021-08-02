#!/usr/bin/Rscript


install.packages("devtools")
library(devtools)

ensure_package_installed <- function (package, repos = repos) {
  if(!require(package, character.only=TRUE)) {
    install.packages(package, repos = repos)
    library(package, character.only=TRUE)
  }
}

ensure_package_installed_with_version <- function (package, version, repos = repos) {
  if(!require(package, character.only=TRUE)) {
    install_version(package, version = version, repos = repos)
    library(package, character.only=TRUE)
  } else if (packageVersion(package) != version) {
    install_version(package, version = version, repos = repos)
    library(package, character.only=TRUE)
  }

  if(packageVersion(package) != version) {
    stop("issue with: ", package, " ", version)
  }

}

# 
# ensure_package_installed("grDevices", repos = list("http://cran.rstudio.com/", "https://cran.ms.unimelb.edu.au/"))
# ensure_package_installed("devtools", repos = list("http://cran.rstudio.com/", "https://cran.ms.unimelb.edu.au/"))
# 
# install.packages("githubinstall", repos = "http://cran.rstudio.com/")
# githubinstall::gh_install_packages("sjmgarnier/viridis", ask = FALSE)
# 
# ensure_package_installed_with_version("rlang", "0.4.10",repos = "http://cran.rstudio.com/")
# 
# ensure_package_installed("parallel", repos = "http://cran.rstudio.com/")
# ensure_package_installed("doParallel", repos = "http://cran.rstudio.com/")
# ensure_package_installed("GGally", repos = "http://cran.rstudio.com/")
# ensure_package_installed("gridExtra", repos = "http://cran.rstudio.com/")
# ensure_package_installed("rmarkdown", repos = "http://cran.rstudio.com/")
# ensure_package_installed("Hmisc", repos = "http://cran.rstudio.com/")
# ensure_package_installed("broom", repos = "http://cran.rstudio.com/")
# ensure_package_installed("FNN", repos = "http://cran.rstudio.com/")
# ensure_package_installed("R.utils", repos = "http://cran.rstudio.com/")
# ensure_package_installed("rmarkdown", repos = "http://cran.rstudio.com/")
# ensure_package_installed("corrplot", repos = "http://cran.rstudio.com/")
# ensure_package_installed("umap", repos = "http://cran.rstudio.com/")
# ensure_package_installed("ggpubr", repos = "http://cran.rstudio.com/")
# ensure_package_installed("stringi", repos = "http://cran.rstudio.com/")
# 
# ensure_package_installed("bit64", repos = "http://cran.rstudio.com/")
# ensure_package_installed("snow", repos = "http://cran.rstudio.com/")
# ensure_package_installed("tinytex", repos = "http://cran.rstudio.com/")
# 
# ensure_package_installed_with_version("data.table", "1.14.0",repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("foreach", "1.5.1", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("stringr", "1.4.0", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("jsonlite", "1.7.2", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("plotly", "4.9.4.1", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("ggplot2", "3.3.5", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("knitr", "1.33", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("rmarkdown","2.9", repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("FactoMineR", '2.4', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("factoextra", '1.0.7', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("testthat", '3.0.4', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("ggrepel", '0.9.1', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("Hmisc", '4.5-0', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("corrplot", '0.90', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("Rcpp", '1.0.7', repos = "http://cran.rstudio.com/")
# ensure_package_installed_with_version("RcppArmadillo", '0.10.6.0.0', repos = "http://cran.rstudio.com/")
# 
# 
# ensure_package_installed_with_version("readr", '2.0.0', repos = "http://cran.rstudio.com/")
# 
# #### Install packages Bioconductor
# ensure_package_installed_with_version("BiocManager", '1.30.16', repos = "http://cran.rstudio.com/")
# BiocManager::install(c("Biobase", "limma", "qvalue"))
# 
# biocManager_valid <- BiocManager::valid()
# if (typeof(biocManager_valid) == 'list'){
#   BiocManager::valid()$out_of_date
#   print(BiocManager::valid()$out_of_date)
#   sessionInfo()
#   stop("BiocManager is not valid")
# }


ensure_package_installed_with_version("BiocManager", '1.30.16', repos = "http://cran.rstudio.com/")
install.packages("renv")

renv::restore(library = '/usr/local/lib64/R/library/')


sessionInfo()