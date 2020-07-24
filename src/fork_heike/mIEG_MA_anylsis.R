### mIEG Analysis Script

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
  filter(name == "dat_ma_m.csv") %>%
  osf_download(overwrite = TRUE)

project %>% 
  osf_ls_files() %>%
  filter(name == "dat_ma_f.csv") %>%
  osf_download(overwrite = TRUE)

dat_ma_m <- read.csv('dat_ma_m.csv')
dat_ma_f <- read.csv('dat_ma_f.csv')


# META1 - MALE -----------------------------------------------------------

### DESCRIPTIVES
#Number of animals.
nCtlTot <- sum(dat_ma_m$nC[unique(dat_ma_m$nest)])
nExpTot <- sum(dat_ma_m$nE[unique(dat_ma_m$nest)])
#Number of publications.
length(unique(dat_ma_m$ID))
#Percentage rats.
dat_ma_m %>% count(species)
dat_ma_m %>% count(model)

dat_ma_m %>% count(areaLevel2,tAcuteStressor)

#Calculate SEM into SD.
levels(as.factor(c(as.character(dat_ma_m$varEtype),
                   as.character(dat_ma_m$varCtype)))) # if NA, assumed SEM
dat_ma_m$sdC <- as.numeric(as.character(dat_ma_m$varC)) * 
  sqrt(as.numeric(as.character(dat_ma_m$nC)))
dat_ma_m$sdE <- as.numeric(as.character(dat_ma_m$varE)) * 
  sqrt(as.numeric(as.character(dat_ma_m$nE)))

#Calculate ES.
dat_ma_m <- escalc(m1i = as.numeric(avgE), sd1i = sdE, 
                   n1i = as.numeric(as.character(dat_ma_m$nE)), 
                   m2i = as.numeric(avgC), sd2i = sdC, 
                   n2i = as.numeric(as.character(dat_ma_m$nC)), 
                   measure = "SMD",
                   data = dat_ma_m)

#Add cohens d from stats papers only & impute var. 
dat_ma_m$yi[is.na(dat_ma_m$yi)] <- dat_ma_m$cohensD[is.na(dat_ma_m$yi)]
dat_ma_m$vi[is.na(dat_ma_m$vi)] <- mean(dat_ma_m$vi, na.rm = T)

#Update var types.
dat_ma_m$areaLevel1     <- as.factor(dat_ma_m$areaLevel1)
dat_ma_m$areaLevel2     <- as.factor(dat_ma_m$areaLevel2)
dat_ma_m$hit2           <- as.factor(dat_ma_m$hit2)
dat_ma_m$tStressorType  <- as.factor(dat_ma_m$tStressorType)
dat_ma_m$tAcuteStressor <- as.factor(dat_ma_m$tAcuteStressor)

#Create z-score variable.
dat_ma_m$zi <- (as.numeric(dat_ma_m$yi) - mean(as.numeric(dat_ma_m$yi), na.rm = T)) / 
  sd(as.numeric(dat_ma_m$yi), na.rm = T)
dat_ma_m <- dat_ma_m[dat_ma_m$zi < 3.29,]

#Build model.
mod_ma_m <- rma.mv(yi, vi, 
                       random = list(~1 | each, ~1 | nest),
                       method = "REML",
                       data = dat_ma_m,
                       slab=paste(authors, year, sep=", "))
summary(mod_ma_m)

#Add moderators. 
mod_ma_m <- rma.mv(yi, vi, 
                      random = list(~1 | each, ~1 | nest),
                      method = "REML",
                      data = dat_ma_m,
                      mods = ~ areaLevel2:tAcuteStressor:hit2 - 1,
                      slab = paste(authors, year, sep=", "))
summary(mod_ma_m)

#Statistical Analysis.
res_hit <- anova(mod_ma_m, L = rbind(c(0,0,0,0,0,0,0,0,
                            0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),  #2nd hit to 0
                          c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,
                            0,0,0,0,0,0,0,0,0,0))) #No 2nd hit to 0

res_AS <- anova(mod_ma_m, L = rbind(c(0.125,0.125,0.125,0,0,0,0,0,0.125,0.125,0.125,0.125,0.125,0,0,0,0,0), #BL to 0
                          c(0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1))) # AS vs 0

res_area <- anova(mod_ma_m, L = rbind(c(0.25,0,0,0.25,0,0,0,0,0.25,0,0,0,0,0.25,0,0,0,0),  #AMY
                          c(0,0.25,0,0,0.25,0,0,0,0,0.25,0,0,0,0,0.25,0,0,0),  #HPF
                          c(0,0,0.25,0,0,0.25,0,0,0,0,0.25,0,0,0,0,0.25,0,0),  #HYP
                          c(0,0,0,0,0,0,(1/3),0,0,0,0,(1/3),0,0,0,0,(1/3),0),  #PFC
                          c(0,0,0,0,0,0,0,(1/3),0,0,0,0,(1/3),0,0,0,0,(1/3)))) # PVT


# PLOTS META1 - MALE ------------------------------------------------------

#Barplot area.
dat_plt        <- cbind(res_area$Lb,res_area$se,
                        c(c('AMY','HPF','HYP','PFC','TH')))
dat_plt        <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','se','BA')
dat_plt$yi     <- as.numeric(as.character(dat_plt$yi))
dat_plt$se     <- as.numeric(as.character(dat_plt$se))

ggplot(dat_plt) + 
  geom_bar(aes(BA,yi), stat = 'identity', alpha=0.5, width=0.7, color = 'black', fill = 'white') +
  geom_errorbar(aes(BA, ymin = yi-1.96*se, ymax = yi+1.96*se), color = 'black', width = 0.3) +
  scale_x_discrete(labels = c('Amygdala','Hippocampus','Hypothalamus','Medial PFC','Thalamus')) +
  scale_y_continuous(limits = c(-0.5,1)) +
  ylab("Hedge's g") +
  geom_vline(aes(xintercept=1.5), colour = 'gray') +
  geom_vline(aes(xintercept=2.5), colour = 'gray') +
  geom_vline(aes(xintercept=3.5), colour = 'gray') +
  geom_vline(aes(xintercept=4.5), colour = 'gray') +
  geom_hline(aes(yintercept=0)) + 
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'gray'),
        panel.grid.minor.y = element_line(colour = 'gray'),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.position = 'none',
        axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20,angle = 45, hjust = 1),
        axis.title.y = element_text(size=25, vjust=5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'))

#Barplot hit.
dat_plt        <- cbind(res_hit$Lb,res_hit$se,
                        c(c('1','0')))
dat_plt        <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','se','Hit')
dat_plt$yi     <- as.numeric(as.character(dat_plt$yi))
dat_plt$se     <- as.numeric(as.character(dat_plt$se))

ggplot(dat_plt) + 
  geom_bar(aes(Hit,yi), stat = 'identity', alpha=0.5, width=0.7, color = 'black', fill = 'white') +
  geom_errorbar(aes(Hit, ymin = yi-1.96*se, ymax = yi+1.96*se), color = 'black', width = 0.3) +
  scale_x_discrete(labels = c('No Second Hit','Second Hit')) +
  scale_y_continuous(limits = c(-0.5,1)) +
  ylab("Hedge's g") +
  geom_vline(aes(xintercept=1.5), colour = 'gray') +
  geom_vline(aes(xintercept=2.5), colour = 'gray') +
  geom_vline(aes(xintercept=3.5), colour = 'gray') +
  geom_vline(aes(xintercept=4.5), colour = 'gray') +
  geom_hline(aes(yintercept=0)) + 
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'gray'),
        panel.grid.minor.y = element_line(colour = 'gray'),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.position = 'none',
        axis.text.y = element_text(size=30), 
        axis.text.x = element_text(size=20,angle = 45, hjust = 1),
        axis.title.y = element_text(size=30, vjust=5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'))

#Barplot AS.
dat_plt        <- cbind(res_AS$Lb,res_AS$se,
                        c(c('0','1')))
dat_plt        <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','se','Hit')
dat_plt$yi     <- as.numeric(as.character(dat_plt$yi))
dat_plt$se     <- as.numeric(as.character(dat_plt$se))

ggplot(dat_plt) + 
  geom_bar(aes(Hit,yi), stat = 'identity', alpha=0.5, width=0.7, color = 'black', fill = 'white') +
  geom_errorbar(aes(Hit, ymin = yi-1.96*se, ymax = yi+1.96*se), color = 'black', width = 0.3) +
  scale_x_discrete(labels = c('Baseline','Acute Stressor')) +
  scale_y_continuous(limits = c(-0.5,1)) +
  ylab("Hedge's g") +
  geom_vline(aes(xintercept=1.5), colour = 'gray') +
  geom_vline(aes(xintercept=2.5), colour = 'gray') +
  geom_vline(aes(xintercept=3.5), colour = 'gray') +
  geom_vline(aes(xintercept=4.5), colour = 'gray') +
  geom_hline(aes(yintercept=0)) + 
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'gray'),
        panel.grid.minor.y = element_line(colour = 'gray'),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.position = 'none',
        axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20,angle = 45, hjust = 1),
        axis.title.y = element_text(size=25, vjust=5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'))


# META2 - FEMALE ----------------------------------------------------------

### DESCRIPTIVES
#Number of animals.
nCtlTot <- sum(dat_ma_f$nC[unique(dat_ma_f$nest)])
nExpTot <- sum(dat_ma_f$nE[unique(dat_ma_f$nest)])
#Number of publications.
length(unique(dat_ma_f$ID))
#Percentage rats.
dat_ma_f %>% count(species)
dat_ma_f %>% count(model)

dat_ma_f %>% count(tAcuteStressor)
dat_ma_f %>% count(areaLevel2)
dat_ma_f %>% count(hit2)#Remove PFC and only after AS. 

#dat_ma_f <- dat_ma_f[dat_ma_f$tAcuteStressor == 1,]
dat_ma_f <- dat_ma_f[dat_ma_f$areaLevel2 %in% c('AMY','HPF','HYP','TH'),]

#Calculate SEM into SD.
levels(as.factor(c(as.character(dat_ma_f$varEtype),
                   as.character(dat_ma_f$varCtype)))) # if NA, assumed SEM
dat_ma_f$sdC <- as.numeric(as.character(dat_ma_f$varC)) * 
  sqrt(as.numeric(as.character(dat_ma_f$nC)))
dat_ma_f$sdE <- as.numeric(as.character(dat_ma_f$varE)) * 
  sqrt(as.numeric(as.character(dat_ma_f$nE)))

#Calculate ES.
dat_ma_f <- escalc(m1i = as.numeric(avgE), sd1i = sdE, 
                   n1i = as.numeric(as.character(dat_ma_f$nE)), 
                   m2i = as.numeric(avgC), sd2i = sdC, 
                   n2i = as.numeric(as.character(dat_ma_f$nC)), 
                   measure = "SMD",
                   data = dat_ma_f)

#Add cohens d from stats papers only & impute var. 
dat_ma_f$yi[is.na(dat_ma_f$yi)] <- dat_ma_f$cohensD[is.na(dat_ma_f$yi)]
dat_ma_f$vi[is.na(dat_ma_f$vi)] <- mean(dat_ma_f$vi, na.rm = T)

#Update var types.
dat_ma_f$areaLevel1     <- as.factor(dat_ma_f$areaLevel1)
dat_ma_f$areaLevel2     <- as.factor(dat_ma_f$areaLevel2)
dat_ma_f$hit2           <- as.factor(dat_ma_f$hit2)
dat_ma_f$tStressorType  <- as.factor(dat_ma_f$tStressorType)
dat_ma_f$tAcuteStressor <- as.factor(dat_ma_f$tAcuteStressor)

#Create z-score variable.
dat_ma_f$zi <- (as.numeric(dat_ma_f$yi) - mean(as.numeric(dat_ma_f$yi), na.rm = T)) / 
  sd(as.numeric(dat_ma_f$yi), na.rm = T)
dat_ma_f <- dat_ma_f[dat_ma_f$zi < 3.29,]

#Build model.
mod_ma_f <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_ma_f,
                   slab=paste(authors, year, sep=", "))
summary(mod_ma_f)

#Add moderators. 
mod_ma_f <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_ma_f,
                   mods = ~ hit2:tAcuteStressor:areaLevel2 -1,
                   slab = paste(authors, year, sep=", "))
summary(mod_ma_f)


# PLOTS - META2 FEMALE ----------------------------------------------------

#Barplot area.
dat_plt        <- cbind(mod_ma_f$b,mod_ma_f$ci.lb,mod_ma_f$ci.ub)
dat_plt        <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','ci.lb','ci.ub')
dat_plt$yi     <- as.numeric(as.character(dat_plt$yi))
dat_plt$ci.lb  <- as.numeric(as.character(dat_plt$ci.lb))
dat_plt$ci.ub  <- as.numeric(as.character(dat_plt$ci.ub))

ggplot(dat_plt) + 
  geom_bar(aes(1, yi), stat = 'identity', alpha=0.5, width=0.7, color = 'rosybrown', fill = 'rosybrown') +
  geom_errorbar(aes(1, ymin = ci.lb, ymax = ci.ub), color = 'rosybrown', width = 0.3) +
  scale_x_discrete(labels = c('Amygdala','Hippocampus','Hypothalamus','Thalamus')) +
  scale_y_continuous(limits = c(-1.1,0.75)) +
  ylab("Hedge's g") +
  geom_vline(aes(xintercept=1.5), colour = 'gray') +
  geom_vline(aes(xintercept=2.5), colour = 'gray') +
  geom_vline(aes(xintercept=3.5), colour = 'gray') +
  geom_vline(aes(xintercept=4.5), colour = 'gray') +
  geom_hline(aes(yintercept=0)) + 
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'gray'),
        panel.grid.minor.y = element_line(colour = 'gray'),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        legend.position = 'none',
        axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20,vjust=6),
        axis.title.y = element_text(size=25, vjust=5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'))


### END ###
