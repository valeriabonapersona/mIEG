### mIEG - PREPROCESSING SCRIPT

#Prepare environment. 
rm(list = ls())
set.seed(3351)

#Libraries.
library(tidyverse)
library(metafor)
library(osfr)
library(readxl)
library(reshape2)
library(wesanderson)
library(ggplot2)

# LOAD DATA ---------------------------------------------------------------

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

#Check dataset IDs for coherence.
identical(unique(pubs$ID), unique(exps$ID))
identical(unique(exps$nest), unique(outs$nest))


# RE-COMBINE DATA ---------------------------------------------------------

#Combine sheets into one dataset from which to create sub-datasets.
dat <- data.frame(pubs[match(exps$ID, pubs$ID),], exps)
dat <- data.frame(dat[match(outs$nest, dat$nest),], outs)

### DATASET 1: Meta-Analysis Males Only

#Extract data that we will used for MA only. The other data will only be 
#evaluated on SR level (see below).

dat_ma_m <- dat[dat$meta == 1,] #cFos and stressors(var meta).
dat_ma_m <- dat_ma_m[dat_ma_m$sex == 'M',] #male only.
#areasLevel2: PFC, TH, HPF, HYP, PFC
dat_ma_m <- dat_ma_m[dat_ma_m$areaLevel2 %in% c('PFC','TH','HPF','HYP','AMY'),] 
#areasLevel1: ACA, PVT, CA1, CA3, DG, PVN, CEA, BLA
#dat_ma_m <- dat_ma_m[dat_ma_m$areaLevel1 %in% c('ACA','PL','IL',
#                                                 'CA1','CA2','CA3','DG',
#                                                 'vCA1','vCA2','vCA3','vDG',
#                                                 'dCA1','dCA2','dCA3','dDG',
#                                                 'PVH','PVT',
#                                                 'CEA','BLA','MEA'),] 
#rename CA1, CA3, DG consistently (remove d, v)
dat_ma_m$areaLevel1[dat_ma_m$areaLevel1 %in% c('vCA1','vCA2','vCA3','vDG','dCA1','dCA2','dCA3','dDG')] <-
  substring(dat_ma_m$areaLevel1[dat_ma_m$areaLevel1 %in% c('vCA1','vCA2','vCA3','vDG','dCA1','dCA2','dCA3','dDG')],2)
#select variables of interest only. 
dat_ma_m <- dat_ma_m[,c('ID','nest','each', #id vars
                        'authors','year', #for labs
                        'species','strain','origin', #animal info
                        'model','mCage','mTimeStart','mTimeLength','mHoursTotal', #model info
                        'hit2', #second hit
                        'tAcuteStressor','tStressorType','tNovel','tWaitPeriod', #acute stressor info
                        'outMeasure','outTechnique','outUnit', #outcome measures
                        'areaLevel1','areaLevel2', #brain areas
                        'nC','avgC','avgCtype','varC','varCtype', #control data
                        'nE','avgE','avgEtype','varE','varEtype', #experimental data
                        'dataOrigin','cohensD','bias','blindRand','sigEffect')]
dat_ma_m$authors <- gsub("^(.*?),.*", "\\1 et al.",dat_ma_m$authors)
rownames(dat_ma_m) <- seq(1,nrow(dat_ma_m))
#save file. Upload to osf.
write.csv(dat_ma_m, 'dat_ma_m.csv')
osf_upload(project, 'dat_ma_m.csv', overwrite = T)


### DATASET 2: Females, with same data as above

dat_ma_f <- dat[dat$meta == 1,] #cFos and stressors(var meta).
dat_ma_f <- dat_ma_f[dat_ma_f$sex == 'F',] #male only.
#areasLevel2: PFC, TH, HPF, HYP, PFC
dat_ma_f <- dat_ma_f[dat_ma_f$areaLevel2 %in% c('PFC','TH','HPF','HYP','PFC','AMY'),] 
#areasLevel1: ACA, PVT, CA1, CA3, DG, PVN, CEA, BLA
#dat_ma_f <- dat_ma_f[dat_ma_f$areaLevel1 %in% c('ACA','PL','IL',
#                                                'CA1','CA2','CA3','DG',
#                                                'vCA1','vCA2','vCA3','vDG',
#                                                'dCA1','dCA2','dCA3','dDG',
#                                                'PVH','PVT',
#                                                'CEA','BLA','MEA'),] 
#rename CA1, CA3, DG consistently (remove d, v)
dat_ma_f$areaLevel1[dat_ma_f$areaLevel1 %in% c('vCA1','vCA3','vDG','dCA1','dCA3','dDG')] <-
  substring(dat_ma_f$areaLevel1[dat_ma_f$areaLevel1 %in% c('vCA1','vCA3','vDG','dCA1','dCA3','dDG')],2)
#select variables of interest only. 
dat_ma_f <- dat_ma_f[,c('ID','nest','each', #id vars
                        'authors','year', #for labs
                        'species','strain','origin', #animal info
                        'model','mCage','mTimeStart','mTimeLength','mHoursTotal', #model info
                        'hit2', #second hit
                        'tAcuteStressor','tStressorType','tNovel','tWaitPeriod', #acute stressor info
                        'outMeasure','outTechnique','outUnit', #outcome measures
                        'areaLevel1','areaLevel2', #brain areas
                        'nC','avgC','avgCtype','varC','varCtype', #control data
                        'nE','avgE','avgEtype','varE','varEtype', #experimental data
                        'dataOrigin','cohensD','bias','blindRand','sigEffect')]
dat_ma_f$authors <- gsub("^(.*?),.*", "\\1 et al.",dat_ma_f$authors)
rownames(dat_ma_f) <- seq(1,nrow(dat_ma_f))
#save file. Upload to osf.
write.csv(dat_ma_f, 'dat_ma_f.csv')
osf_upload(project, 'dat_ma_f.csv', overwrite = T)

### DATASET 3: SR Rest

dat_sr <- dat[dat$meta == 0,] #cFos and stressors(var meta).
#areasLevel2: PFC, TH, HPF, HYP, PFC
dat_sr <- rbind(dat_sr,dat[!(dat$areaLevel2 %in% c('PFC','TH','HPF','HYP','AMY')),] )
dat_sr <- unique(dat_sr)

#Calculate SEM into SD.
levels(as.factor(c(as.character(dat_sr$varEtype),
                   as.character(dat_sr$varCtype)))) # if NA, assumed SEM
dat_sr$sdC <- as.numeric(as.character(dat_sr$varC)) * 
   sqrt(as.numeric(as.character(dat_sr$nC)))
dat_sr$sdE <- as.numeric(as.character(dat_sr$varE)) * 
   sqrt(as.numeric(as.character(dat_sr$nE)))

dat_sr$ttest <- tsum.test(dat_sr$avgC, dat_sr$sdC, dat_sr$nC, dat_sr$avgE, dat_sr$sdE, dat_sr$nE)$p.value
dat_sr_cfos <- dat_sr[dat_sr$iegName == 'cFos',]
dat_sr_cfos_stressor <- dat_sr_cfos[dat_sr_cfos$meta == '0',]
dat_sr_cfos_ba <- dat_sr_cfos[dat_sr_cfos$meta == '1',]
write.csv(dat_sr_cfos_ba, 'dat_sr_cfos_ba.csv')
write.csv(dat_sr_cfos_stressor, 'dat_sr_cfos_stressor.csv')

dat_sr_fosB <- dat_sr[dat_sr$iegName == 'dFosB',]
dat_sr_arc <- dat_sr[dat_sr$iegName == 'Arc',]
dat_sr_egr <- dat_sr[dat_sr$iegName %in% c('Egr1','Egr2','Egr4'),]

write.csv(dat_sr_fosB, 'dat_sr_fosB.csv')
write.csv(dat_sr_arc, 'dat_sr_arc.csv')
write.csv(dat_sr_egr, 'dat_sr_egr.csv')

#Check that all data is accounted for.
nrow(dat_sr) + nrow(dat_ma_m) + nrow(dat_ma_f) == nrow(dat)
nrow(dat_sr_fosB) + nrow(dat_sr_arc) + nrow(dat_sr_cfos_ba) + nrow(dat_sr_cfos_stressor) + 
   nrow(dat_sr_egr) == nrow(dat_sr)


# DESCRPTIVES OF DATASET --------------------------------------------------

#total pubs, exps, outs
nrow(pubs)
nrow(exps)
nrow(outs)

#n animals total
nests <- subset(dat, !duplicated(nest))
sum(as.numeric(nests$nC), na.rm = TRUE) + sum(as.numeric(nests$nE), na.rm = TRUE)

#percentage significant
sum(dat$sigEffect)/nrow(dat)

#percentage rand and blind
pubs <- subset(dat, !duplicated(ID))
length(which(pubs$blindRand == 1))/35

#prepare ROB data for plots
dat_rob <- dat[,c('ID','nest', 'bias', 'blindRand', 'meta',
                 'RB_seqGeneration','RB_baseline','RB_allocation',
                 'RB_housing','RB_blindExp','RB_outAss','RB_outBlind',
                 'RB_incData','RB_control')]
dat_rob <- subset(dat_rob, !duplicated(nest))


### END ### 

