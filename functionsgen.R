rm(list=ls(all=TRUE)) # removes existing data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require("devtools")
require("usethis")
# install.packages("roxygen2")
require("roxygen2")
require("nealeLabFunctions")
devtools::load_all()
devtools::document()
# devtools::load_all()

getwd()
roxygenise()

install_github("WHG1990/CCTools")

library(CCTools)


CCTools::
