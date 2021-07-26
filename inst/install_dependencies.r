# production environment install script

cat("Preparing R enviroment")

list.of.packages <- c("devtools", "webshot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(devtools)

#expects to be run from .. (home dir)
devtools::install_deps(pkg = ".", upgrade = "never", repos = "http://cran.us.r-project.org")

# for pdf exports
webshot::install_phantomjs()

