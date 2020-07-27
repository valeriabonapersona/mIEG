## Script: "Download data"

## project: mIEG

## Author: Valeria Bonapersona
## contact: v.bonapersona-2 (at) umcutrecht.nl

## DOES NOT WORK

# Environment preparation -------------------------------------------------
source("src/utilities.R")




# Connection to open science framework ------------------------------------

OSF_PAT <- "FTQ6344mqWM2ahibwpfi9sWezdt8KK8FH7P4saUN2AsxnsQ3cZETNFXzfqWHXQVmyKJX7m" #V: Change PAT.
osf_auth(OSF_PAT)

project <- osf_retrieve_node("t7xcu")

project %>% 
  osf_ls_files() %>%
  filter(name == "mIEG_outcomes.xlsx") %>%
  osf_download(overwrite = TRUE)

pubs <- read_excel("mIEG_outcomes.xlsx", sheet = "Publications")
exps <- read_excel("mIEG_outcomes.xlsx", sheet = "Experiments")
outs <- read_excel("mIEG_outcomes.xlsx", sheet = "Outcomes")