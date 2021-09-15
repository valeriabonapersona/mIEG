## Script: "Download your dependencies"

## project: mIEG

## Author: Valeria Bonapersona
## contact: v.bonapersona-2 (at) umcutrecht.nl



# Environment preparation -------------------------------------------------
rm(list = ls())

# Install packages
## from CRAN
list_cran_packages <- c(
  "tidyverse", "readxl", "stringr", # data handling -> should all be in tidyverse but dont work
  "remotes", # download from github
  "metafor", # estimates and meta-analysis
  "BSDA", # t-test from summary statistics
  "multcomp", # contrasts with multiple testing correction
  "ggplot2", 
  "ggpubr",
  "svglite"# to combine multiple plots
)


new_packages <- list_cran_packages[!(list_cran_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

## from Github
if(!require(osfr)) remotes::install_github("centerforopenscience/osfr")


# library
lapply(list_cran_packages, library, character.only = TRUE)

# Packrat for dependenies ------------------------------------------------
#packrat::init("src/")

