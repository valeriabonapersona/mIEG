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


# EXPLORE HIT2 ------------------------------------------------------------

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


#Build model - Strain.
dat_new <- dat_ma_m[dat_ma_m$strain %in% c('SD','WI'),]
mod_ma_m <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_new,
                   mods = ~ strain -1,
                   slab=paste(authors, year, sep=", "))
summary(mod_ma_m)
dat_new %>% count(strain)

dat_plt <- cbind(mod_ma_m$b,mod_ma_m$ci.lb,mod_ma_m$ci.ub,c('SD','WI'))
dat_plt <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','ci.lb','ci.ub','strain')
dat_plt$yi <- as.numeric(as.character(dat_plt$yi))
dat_plt$ci.lb <- as.numeric(as.character(dat_plt$ci.lb))
dat_plt$ci.ub <- as.numeric(as.character(dat_plt$ci.ub))

ggplot(dat_plt) + 
  geom_bar(aes(as.factor(strain),yi), color='black',fill='white', stat = 'identity', 
           position = position_dodge(width=0.8), alpha=0.5, width=0.7) + 
  geom_errorbar(aes(as.factor(strain), ymin = ci.lb, ymax = ci.ub), 
                position = position_dodge(width=0.8), width = 0.5) +
  scale_x_discrete(labels = c('Sprague-Dawley','Wistar')) +
    scale_y_continuous(limits = c(-0.5,1.5)) +
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

#Build model - Outcome Measure.
dat_new <- dat_ma_m[dat_ma_m$outMeasure %in% c('P','M'),]
mod_ma_m <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_new,
                   mods = ~ outMeasure -1,
                   slab=paste(authors, year, sep=", "))
summary(mod_ma_m)
dat_new %>% count(strain)

dat_plt <- cbind(mod_ma_m$b,mod_ma_m$ci.lb,mod_ma_m$ci.ub,c('mRNA','Protein'))
dat_plt <- as.data.frame(dat_plt)
names(dat_plt) <- c('yi','ci.lb','ci.ub','outMeasure')
dat_plt$yi <- as.numeric(as.character(dat_plt$yi))
dat_plt$ci.lb <- as.numeric(as.character(dat_plt$ci.lb))
dat_plt$ci.ub <- as.numeric(as.character(dat_plt$ci.ub))

ggplot(dat_plt) + 
  geom_bar(aes(as.factor(outMeasure),yi), color='black',fill='white', stat = 'identity', 
           position = position_dodge(width=0.8), alpha=0.5, width=0.7) + 
  geom_errorbar(aes(as.factor(outMeasure), ymin = ci.lb, ymax = ci.ub), 
                position = position_dodge(width=0.8), width = 0.5) +
  scale_x_discrete(labels = c('mRNA','Protein')) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
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




# EXPLORE BIAS ------------------------------------------------------------

par(mar=c(5,5,1,2))

dat_ma_bias <- rbind(dat_ma_m, dat_ma_f)
dat_ma_bias$yi <- abs(as.numeric(dat_ma_bias$yi))
mod_ma_m <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_ma_bias,
                   mods = ~ bias,
                   slab=paste(authors, year, sep=", "))
summary(mod_ma_m)
preds <- predict(mod_ma_m,newmods=c(0:10))

plot(NA, NA, xlim=c(0,6), ylim=c(0,4),
     xlab="Risk of Bias", ylab="Effect Size", 
     bty = 'n')

size <- 1 / sqrt(na.pass(dat_ma_bias$vi)) 
size <- size / max(size/0.05)

### add points
symbols(dat_ma_bias$bias, dat_ma_bias$yi, circles=size, inches=FALSE, add=TRUE, bg='black')

### add predicted values (and corresponding CI bounds)
lines(0:10, preds$pred)
lines(0:10, preds$ci.lb, lty="dashed")
lines(0:10, preds$ci.ub, lty="dashed")

cor(dat_ma_bias$yi, dat_ma_bias$bias)


par(mar=c(5,5,1,2))

dat_new <- dat_ma_m
dat_new$yi <- abs(as.numeric(dat_new$yi))
mod_ma_m <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = dat_new,
                   mods = ~ year,
                   slab=paste(authors, year, sep=", "))
summary(mod_ma_m)
preds <- predict(mod_ma_m,newmods=c(2000:2020))

plot(NA, NA, xlim=c(2000,2020), ylim=c(0,4),
     xlab="Year", ylab="Effect Size", 
     bty = 'n')

size <- 1 / sqrt(na.pass(dat_new$vi)) 
size <- size / max(size/0.2)

### add points
symbols(dat_new$year, dat_new$yi, circles=size, inches=FALSE, add=TRUE, bg='black')

### add predicted values (and corresponding CI bounds)
lines(2000:2020, preds$pred)
lines(2000:2020, preds$ci.lb, lty="dashed")
lines(2000:2020, preds$ci.ub, lty="dashed")

