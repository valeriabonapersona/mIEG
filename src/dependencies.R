## Script: "Download your dependencies"

## project: mIEG

## Author: Valeria Bonapersona
## contact: v.bonapersona-2 (at) umcutrecht.nl



# Environment preparation -------------------------------------------------
rm(list = ls())

# Install packages
## from CRAN
list_cran_packages <- c(
  "tidyverse", # data handling
  "remotes" # download from github
)


new_packages <- list_cran_packages[!(list_cran_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

## from Github
if(!require(osfr)) remotes::install_github("centerforopenscience/osfr")


# Packrat for dependenies ------------------------------------------------
#packrat::init("src/")

