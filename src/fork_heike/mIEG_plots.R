
# ROB PLOT ----------------------------------------------------------------

#Prepare labels.
text_fill <- c('No','Unclear','Yes','Computerized')
text_x    <- c('0','25','50','75','100')
text_y    <- c('Sequence Generation',
               'Baseline Characteristics',
               'Allocation Concealment',
               'Random Housing',
               'Blinded Experimenter',
               'Random Allocation',
               'Binded Assessment',
               'Incomplete Data Addressed',
               'Adequate Control Group')
col <- c('#1a9641', '#a6d96a', '#fdae61', '#d7191c')

#Melt data for ggplot.
dat_rob %>%
  melt(id.vars = c('ID', 'nest', 'bias', 'blindRand', 'meta')) -> dat_rob
dat_rob$value <- factor(dat_rob$value, levels = rev(c('C','L','UC','H')))
dat_rob$variable <- factor(dat_rob$variable, levels = rev(levels(dat_rob$variable)))

dat_rob %>% count(value) -> UC
UC[2,2] /  72

#Plot.
ggplot() +
  geom_bar(data = dat_rob,
           aes(x = variable, fill = value),
           position = 'fill', width = 0.5, color = 'black')+
  coord_flip()+
  guides(fill=guide_legend(title=element_blank(), reverse = T))+
  scale_fill_manual(values=rev(col),labels=text_fill)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=text_x,
                     limits=c(0,1),
                     expand=c(0,0))+
  ylab('Percentage')+
  scale_x_discrete(labels=rev(text_y))+
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(color='black'),
        axis.title.x=element_text(size=26),
        axis.text=element_text(size=26,color='black'),
        legend.text=element_text(size=26),
        legend.key.size=unit(9,'mm'),
        panel.background=element_blank(),
        legend.position = 'bottom',
        plot.margin = margin(1,1,1,1,'cm'))


# STUDY CHARACTERISTICS PLOT ----------------------------------------------

list(ls = rm())
# @ heike - I did not understan what the levels of the first plot are, so I changed it

pubs <- read_excel(paste0(raw,"outcomes.xlsx"), sheet = "Publications")
exps <- read_excel(paste0(raw,"outcomes.xlsx"), sheet = "Experiments")
outs <- read_excel(paste0(raw,"outcomes.xlsx"), sheet = "Outcomes")

#Prepare labels.
text_x    <- c('0','25','50','75','100')
text_y    <- rev(c('Species','Sex','Model','2nd Hit',
                   'Stressor',
                   'IEG','Measure'))


#SC plot data.
SC <- data.frame(exps[match(outs$nest, exps$nest),], outs)
SC <- SC[,c(1,2,39,3,7,8,23,28,40,41,47)]

SC[,c(seq(4,8))] <- lapply(SC[,c(seq(4,8))], as.factor)

#animal
SC$species[which(SC$species == "rats")] <- "rat"
SC$species        <- factor(SC$species, levels = c('rat','mice'))
SC$sex            <- factor(SC$sex, levels = c('M','F','P','NS'))

#model
SC$model          <- factor(SC$model, levels = c('MS','LG','LBN'))
SC$hit2           <- factor(SC$hit2, levels = c('0','1'))

#stressor
SC$tAcuteStressor <- factor(SC$tAcuteStressor, levels = c('0','1'))

#outcome
SC$iegName[which(SC$iegName %in% c('Egr1','Egr2','Egr4'))] <- 'Egr'
SC$iegName        <- factor(SC$iegName, levels = c('cFos','Arc','dFosB','Egr'))
SC$outMeasure     <- factor(SC$outMeasure, levels = c('P','M'))
SC$areaLevel2[which(SC$areaLevel2 %in% c('CB','CTX','HB','MB','STR', 'PD'))] <- 'OTH'
SC$areaLevel2     <- factor(SC$areaLevel2, levels = c('AMY','HPF','PFC','HYP','TH','OTH'))

SC[,c(seq(4,11))] <- lapply(SC[,c(seq(4,11))], as.numeric)

SC <- melt(SC, id.vars = c('ID','nest','each'))


SC_area <- SC[SC$variable == 'areaLevel2',]
SC_area %>% count(variable,value) -> SC_area
#percent(SC_area$n/sum(SC_area$n)) # does not work
text = rev(c('Amygdala','Hippocampus','Medial PFC','Hypothalamus','Thalamus','Other'))

ggplot(SC_area, aes(x = '', y = n, fill = forcats::fct_rev(as.factor(value))))+
  geom_bar(stat = 'identity', color = 'black') +
  coord_polar('y', start = 0) + 
  scale_fill_brewer(palette = 'Spectral', labels = text)+
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"),
        panel.background=element_blank(),
        axis.text = element_blank(),
        legend.text = element_text (size = 20),
        legend.key.size=unit(10,'mm'),
        legend.title = element_blank(),
        legend.position = 'right')+
  guides(fill=guide_legend(title=element_blank(), reverse = T))+
  geom_text(aes(y = n/6 + (c(0, cumsum(n)[-length(n)]) + (n/3)), 
              label = scales::percent(n/sum(n))), size=8, 
              position = position_nudge(x = 0.2))
  
SC <- SC[!(SC$variable == 'areaLevel2'),]


ggplot() +
  geom_bar(data = SC,
           aes(x = rev(variable), fill = forcats::fct_rev(as.factor(value))),
           position = 'fill', width = 0.5, color = 'black')+
  coord_flip()+
  guides(fill=guide_legend(title=element_blank(), reverse = T))+
  scale_fill_brewer(palette = 'YlGnBu')+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=text_x,
                     limits=c(0,1),
                     expand=c(0,0))+
  ylab('Percentage')+
  scale_x_discrete(labels=text_y)+
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(color='black'),
        axis.title.x=element_text(size=20),
        axis.text=element_text(size=20,color='black'),
        legend.text=element_text(size=20),
        legend.key.size=unit(9,'mm'),
        panel.background=element_blank(),
        legend.position = 'bottom',
        plot.margin = margin (1,1,1,1,'cm'))


# SUMMARY SR --------------------------------------------------------------
## Table with brain area groupings
ba_grouped <- outs %>% 
  
  # get required info
  dplyr::select(areaLevel2, areaPub) %>% 
  rename(ba_grouped = areaLevel2, ba_publication = areaPub) %>% 
  
  # get unique
  unique() %>% 
  
  # grouping as used in analysis
  mutate(
    ba_grouped = ifelse(ba_grouped %in% c('CB','CTX','HB','MB','STR', 'PD'), "Other", ba_grouped),
    
    ba_grouped = case_when(
      ba_grouped %in% c('CB','CTX','HB','MB','STR', 'PD') ~ "Other",
      
      ba_grouped == "STR" ~ "Striatum", 
      ba_grouped == "TH" ~ "Thalamus", 
      ba_grouped == "PFC" ~ "Prefrontal cx", 
      ba_grouped == "HYP" ~ "Hypothalamus", 
      ba_grouped == "AMY" ~ "Amygdala",
      ba_grouped == "HPF" ~ "Hippocampus",
      TRUE ~ ba_grouped
    )) %>% 
  
  # wide format
  group_by(ba_grouped) %>%
  summarize(ba_publication = paste(ba_publication, collapse = "; "))

write.csv(ba_grouped, "results/ba_grouped.csv")


# table acute stressors ---------------------------------------------------
acute_stress_grouped <- exp %>% 
  
  # get required info
  dplyr::select(contains("Stressor")) %>%
  
  unique() %>%
  
  # manual corrections heike
  mutate(
    tAcuteStressor = ifelse(tStressorType == "NA", 0, 1), 
    stressor_intensity = case_when(
      tStressorType %in% c('OFT','DLB','NE','EPM') ~ "mild", # @ heike --> shouldnt this be mild? >> below you had 1
      tStressorType %in% c('CRD','FST','MWM','RS',
                           'FS in inhibitory avoidance task',
                           'Shock in shock-probe burial task') ~"severe", 
      TRUE ~ "mild"
    )
  ) %>%
  
  # remove rest condition
  filter(tAcuteStressor == 1) %>%
  
  # wide format
  group_by(stressor_intensity) %>% 
  summarize(stressor_type = paste(tStressorType, collapse = "; "))
  
  
write.csv(acute_stress_grouped, "results/acute_stress_grouped.csv")


# PREPARATION META --------------------------------------------------------
# @ heike - I do not know what data_MA is, so I couldnt run it

levels(factor(c(data_MA$varEtype,data_MA$varCtype))) # if NA, assumed SEM
data_MA$sdC <- as.numeric(as.character(data_MA$varC)) * 
  sqrt(as.numeric(as.character(data_MA$nC)))
data_MA$sdE <- as.numeric(as.character(data_MA$varE)) * 
  sqrt(as.numeric(as.character(data_MA$nE)))

data_MA <- escalc(m1i = as.numeric(avgE), sd1i = sdE, 
                  n1i = as.numeric(as.character(data_MA$nE)), 
                  m2i = as.numeric(avgC), sd2i = sdC, 
                  n2i = as.numeric(as.character(data_MA$nC)), 
                  measure = "SMD",
                  data = data_MA)

# add variable for stressor intensity
data_MA$stressorIntensity <- 0
data_MA$stressorIntensity[which(data_MA$tStressorType %in% c('OFT','DLB','NE','EPM'))] <- 1
data_MA$stressorIntensity[which(data_MA$tStressorType %in% c('CRD','FST','MWM','RS',
                                                             'FS in inhibitory avoidance task',
                                                             'Shock in shock-probe burial task'))] <- 1
data_MA$stressorIntensity <- as.factor(data_MA$stressorIntensity)

# zscore for outliers
data_MA$zi <- (as.numeric(data_MA$yi) - mean(as.numeric(data_MA$yi), na.rm = T)) / 
  sd(as.numeric(data_MA$yi), na.rm = T)

#direction of effect?

data_MA$areaLevel2    <- as.factor(data_MA$areaLevel2)
data_MA$hit2          <- as.factor(data_MA$hit2)
data_MA$tStressorType <- as.factor(data_MA$tStressorType)

data_MA_F <- data_MA[which(data_MA$sex == 'F'),]
data_MA_M <- data_MA[which(data_MA$sex == 'M'),]

mod_F <- rma.mv(yi, vi, 
                #subset = (areaLevel2 == 'AMY'),
                random = list(~1 | each, ~1 | nest),
                mods   = ~ areaLevel2 : hit2 : stressorIntensity -1,
                method = "REML",
                data = data_MA_F)
summary(mod_F)

#Effect of Hit
anova(mod_F, L = c(-0.5,-0.5,(1/3),(1/3),(1/3)))
#Effect of Amygdala
anova(mod_F, L = rbind(c(0.5,0,0.5,0,0),
                       c(0,0.5,0,0,0.5),
                       c(0,0,0,1,0)))

# with V's processed data
meta <- readRDS("~/surfdrive/Work/PhD/mELA/mIEG/mIEG/data/processed/meta.RDS")
data_MA_M <- meta %>% filter(sex == "M")
mod_M <- rma.mv(yi, vi, 
                #subset = (areaLevel2 == 'AMY'),
                random = list(~1 | each, ~1 | nest),
                mods   = ~ areaLevel2 : hit2 -1,
                method = "REML",
                data = data_MA_M,
                slab = paste(each))
summary(mod_M)

#Effect of Hit
anova(mod_M, L = c(-(1/3),-(1/3),-(1/3),(1/3),(1/3),(1/3)))
#Effect of Amygdala
anova(mod_M, L = rbind(c(0.5,0,0,0.5,0,0),
                       c(0,0.5,0,0,0.5,0),
                       c(0,0,0.5,0,0,0.5)))



# EXPLORE -----------------------------------------------------------------

# Build model for plotting purposes only. No moderators necessary. 
mod_M <- rma.mv(yi, vi, 
                random = list(~1 | each, ~1 | nest),
                method = "REML",
                data = data_MA_M,
                slab = paste(each))

#


forest(mod_M, xlim=c(-16, 8), at=c(-3, -2, -1,  0, 1, 2, 3, 4, 5),
       ilab=cbind(data_MA_M$areaLevel1, data_MA_M$stressorIntensity, data_MA_M$tStressorType),
       ilab.xpos=c(-9.5,-6,-4.5), cex=0.75, ylim=c(-1,105),
       order=order(data_MA_M$areaLevel2,data_MA_M$areaLevel1,data_MA_M$stressorIntensity), 
       rows=c(3:20,23:58,61:88,91:104),
       xlab="hedge's g", mlab="", psize=1)



# EXPLORE CORRELATIONS ----------------------------------------------------

#Calculate SEM into SD.
levels(as.factor(c(as.character(dat$varEtype),
                   as.character(dat$varCtype)))) # if NA, assumed SEM
dat$sdC <- as.numeric(as.character(dat$varC)) * 
  sqrt(as.numeric(as.character(dat$nC)))
dat$sdE <- as.numeric(as.character(dat$varE)) * 
  sqrt(as.numeric(as.character(dat$nE)))

#Calculate ES.
dat <- escalc(m1i = as.numeric(avgE), sd1i = sdE, 
                   n1i = as.numeric(as.character(dat$nE)), 
                   m2i = as.numeric(avgC), sd2i = sdC, 
                   n2i = as.numeric(as.character(dat$nC)), 
                   measure = "SMD",
                   data = dat)

# Scatterplot Hours of MS. 
par(mar=c(5,5,1,2))

dat <- data.frame(exps[match(pubs$ID, exps$ID),], pubs)
dat <- data.frame(outs[match(dat$nest, outs$nest),], dat)

res <- rma(yi, vi, mods = ~ as.numeric(as.character(bias)), data=dat)
preds <- predict(res, newmods=c(0:10))

### set up plot (risk ratios on y-axis, absolute latitude on x-axis)
plot(NA, NA, xlim=c(0,7), ylim=c(-2,3),
     xlab="Year of Publication", ylab="Effect Size", 
     bty = 'n')

size <- 1 / sqrt(na.pass(dat$vi)) 
size <- size / max(size/0.3)

### add points
symbols(dat$bias, dat$yi, circles=rep(0.1,332), inches=FALSE, add=TRUE, bg='black')

### add predicted values (and corresponding CI bounds)
lines(0:60, preds$pred)
lines(0:60, preds$ci.lb, lty="dashed")
lines(0:60, preds$ci.ub, lty="dashed")

### dotted line at RR=1 (no difference between groups)
abline(h=0, lty="dotted")

# Scatterplot Year.
par(mar=c(5,5,1,2))

res <- rma(yi, vi, mods = ~ year, data=dat)
preds <- predict(res, newmods=c(2000:2020))

### set up plot (risk ratios on y-axis, absolute latitude on x-axis)
plot(NA, NA, xlim=c(2000,2020), ylim=c(-2,3),
     xlab="Year of Publication", ylab="Effect Size", 
     bty = 'n')

size <- 1 / sqrt(na.exclude(dat$vi)) 
size <- size / max(size/0.3)

### add points
symbols(dat$year, dat$yi, circles=rep(0.1,332), inches=FALSE, add=TRUE, bg='black')

### add predicted values (and corresponding CI bounds)
lines(2000:2020, preds$pred)
lines(2000:2020, preds$ci.lb, lty="dashed")
lines(2000:2020, preds$ci.ub, lty="dashed")

### dotted line at RR=1 (no difference between groups)
abline(h=0, lty="dotted")


### END ###




# NEW PLOTS ---------------------------------------------------------------


# BAR PLOT MALES ----------------------------------------------------------

dat <- dat_ma_m[dat_ma_m$areaLevel1 == 'DG',]
dat <- cbind(dat$yi, dat$vi, dat$tAcuteStressor, dat$hit2)
dat <- as.data.frame(dat)
names(dat) <- c('yi','vi','tAcuteStressor','hit2')

dat %>% group_by(tAcuteStressor, hit2) %>% summarise(mean(yi), mean(vi)) -> dat_DG



# SR PLOT ARC -------------------------------------------------------------

arc_dat <- dat[dat$iegName == 'Arc',]
egr_dat <- dat[dat$iegName %in% c('Egr1','Egr3','Egr4'),]
fosb_dat <- dat[dat$iegName %in% c('dFosB'),]

#Calculate SEM into SD.
levels(as.factor(c(as.character(arc_dat$varEtype),
                   as.character(arc_dat$varCtype)))) # if NA, assumed SEM
arc_dat$sdC <- as.numeric(as.character(arc_dat$varC)) * 
  sqrt(as.numeric(as.character(arc_dat$nC)))
arc_dat$sdE <- as.numeric(as.character(arc_dat$varE)) * 
  sqrt(as.numeric(as.character(arc_dat$nE)))

arc_dat <- escalc(m1i = as.numeric(avgE), sd1i = sdE, 
                   n1i = as.numeric(as.character(arc_dat$nE)), 
                   m2i = as.numeric(avgC), sd2i = sdC, 
                   n2i = as.numeric(as.character(arc_dat$nC)), 
                   measure = "SMD",
                   data = arc_dat)

arc_dat %>% group_by(areaLevel2,hit2) %>% summarise(mean(yi), mean(vi)) -> arc_dat
