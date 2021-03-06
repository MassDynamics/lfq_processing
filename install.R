#!/usr/bin/Rscript

install.packages("devtools", repos = "http://cran.rstudio.com/")
install.packages("githubinstall", repos = "http://cran.rstudio.com/")
githubinstall::gh_install_packages("sjmgarnier/viridis", ask = FALSE)
install.packages("data.table", repos = "http://cran.rstudio.com/")
install.packages("foreach", repos = "http://cran.rstudio.com/")
install.packages("parallel", repos = "http://cran.rstudio.com/")
install.packages("doParallel", repos = "http://cran.rstudio.com/")
install.packages("ggplot2", repos = "http://cran.rstudio.com/")
install.packages("GGally", repos = "http://cran.rstudio.com/")
install.packages("stringr", repos = "http://cran.rstudio.com/")
install.packages("gridExtra", repos = "http://cran.rstudio.com/")
install.packages("rmarkdown", repos = "http://cran.rstudio.com/")
install.packages("Hmisc", repos = "http://cran.rstudio.com/")
install.packages("broom", repos = "http://cran.rstudio.com/")
install.packages("FNN", repos = "http://cran.rstudio.com/")
install.packages("R.utils", repos = "http://cran.rstudio.com/")
install.packages("rmarkdown", repos = "http://cran.rstudio.com/")
install.packages("knitr", repos = "http://cran.rstudio.com/")
install.packages("Rcpp", repos = "http://cran.rstudio.com/")
install.packages("FactoMineR", repos = "http://cran.rstudio.com/")
install.packages("factoextra", repos = "http://cran.rstudio.com/")
install.packages("corrplot", repos = "http://cran.rstudio.com/")
install.packages("umap", repos = "http://cran.rstudio.com/")
install.packages("ggrepel", repos = "http://cran.rstudio.com/")
install.packages("ggpubr", repos = "http://cran.rstudio.com/")
install.packages("jsonlite", repos = "http://cran.rstudio.com/")
install.packages("plotly", repos = "http://cran.rstudio.com/")
install.packages("statmod", repos = "http://cran.rstudio.com/")


#### Intall packages Bioconductor
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("Biobase", suppressUpdates = TRUE)
BiocInstaller::biocLite("limma", suppressUpdates = TRUE)
BiocInstaller::biocLite("qvalue", suppressUpdates = TRUE)

#######
