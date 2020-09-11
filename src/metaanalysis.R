## Script: "Meta-analysis cFos"

## project: mIEG

## Input: preprocessed data from data_preparation.R adjusted by report_descriptives.rmd
## Output: meta-analytic results

## Author: Valeria Bonapersona & Heike Schuler
## contact: v.bonapersona-2 (at) umcutrecht.nl


# Environment preparation -------------------------------------------------
source("src/utilities.R")
df <- readRDS(paste0(processed, "meta.RDS"))


# Outliers ----------------------------------------------------------------
# Create z-score variable.
df$zi <- (df$yi - mean(df$yi, na.rm = T)) / 
  sd(df$yi, na.rm = T)

df_filtered <- df %>% filter(sex == "M")

#Check for influential outliers and remove.
res <- rma(yi, vi, measure="ZCOR", data=df_filtered)
inf <- influence(res)
df_filtered <- df_filtered[-which(inf$inf$inf == '*'),]


# Model -------------------------------------------------------------------
mod <- rma.mv(yi, vi, 
              random = list(~1 | each, ~1 | nest),
              method = "REML",
              data = df_filtered,
              mods = ~ tAcuteStressor:areaLevel2 -1, 
              slab = paste(authors, year, sep=", "))


# Main analysis -----------------------------------------------------------
## Create contrasts per RQ
interactions <- names(data.frame(mod$X))

# RQ1: effects of ELA on presence acute stressor
acutestressor <- ifelse(str_detect(interactions, "tAcuteStressor1"), 1, 0)
contr_acutestressor <- rbind(
  ifelse(acutestressor == 1, 1/sum(acutestressor == 1), 0), # acute vs 0
  ifelse(acutestressor == 0, 1/sum(acutestressor == 0), 0) # baseline vs 0
)

anova(mod, L = contr_acutestressor) 
summary(glht(mod, linfct = contr_acutestressor, df=df.residual(mod)), test=adjusted('holm'))


# RQ2: effects of ELA on brain area
brainarea <- sub(".*2", "", interactions)
contr_brainarea <- NULL
for (area in unique(brainarea)) {
  x <- ifelse(brainarea == area, 1/sum(brainarea == area), 0)
  contr_brainarea <- rbind(contr_brainarea, x)
}

anova(mod, L = contr_brainarea) 
summary(glht(mod, linfct = contr_brainarea, df=df.residual(mod)), test = adjusted('holm'))


# Subgroup: hit2 ----------------------------------------------------------
# estimates for each subgroup
subgroup_hit0 <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ tAcuteStressor:areaLevel2 -1,
                        subset = (hit2==0),
                        slab = paste(authors, year, sep=", "))
subgroup_hit1 <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ tAcuteStressor:areaLevel2 -1,
                        subset = (hit2==1),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_hit <- data.frame(estimate = c(coef(subgroup_hit0), coef(subgroup_hit1)), 
                             stderror = c(subgroup_hit0$se, subgroup_hit1$se),
                             meta = c(rep("hit0", length(subgroup_hit0$se)),
                                      rep("hit1", length(subgroup_hit1$se))), 
                             stressor = c('BL','AS','BL','AS','BL','AS','BL','BL',
                                          'BL','AS','BL','AS','BL','AS','BL','AS','BL','AS'),
                             tau2 = c(rep(subgroup_hit0$tau2, length(subgroup_hit0$se)),
                                      rep(subgroup_hit1$tau2, length(subgroup_hit1$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_hit, digits=3)

# Wald test: to test the difference between the two
# with(comparison_hit, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))

# meta-analyze interaction 2nd hit by stressor
rma(estimate, sei=stderror, mods = ~ meta:stressor -1, method="FE", data=comparison_hit, digits=3)

# Sensitivity: Species ----------------------------------------------------
# model compared to complete model
mod_sens <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = df_filtered,
                   subset = (species=='rat'),
                   mods = ~ tAcuteStressor:areaLevel2 -1, 
                   slab = paste(authors, year, sep=", "))

# summarize together
comparison_sens <- data.frame(estimate = c(coef(mod), coef(mod_sens)), 
                              stderror = c(mod$se, mod_sens$se),
                              sens = c(rep("all", length(mod$se)),
                                       rep("wo_mice", length(mod_sens$se))), 
                              tau2 = c(rep(mod$tau2, length(mod$se)),
                                       rep(mod_sens$tau2, length(mod_sens$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_sens, digits=3)



# Sensitivity: ELA model --------------------------------------------------
# model compared to complete model
mod_sens <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = df_filtered,
                   subset = (model=='MS'),
                   mods = ~ tAcuteStressor:areaLevel2 -1, 
                   slab = paste(authors, year, sep=", "))

# summarize together
comparison_sens <- data.frame(estimate = c(coef(mod), coef(mod_sens)), 
                              stderror = c(mod$se, mod_sens$se),
                              sens = c(rep("all", length(mod$se)),
                                       rep("MS_only", length(mod_sens$se))), 
                              tau2 = c(rep(mod$tau2, length(mod$se)),
                                       rep(mod_sens$tau2, length(mod_sens$se))))

# meta-analyze with fixed effect
res <- rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_sens, digits=3)


# Exploratory: Stressor Type ----------------------------------------------
#categorize into groups. 
phys <- c('CRD','FS in inhibitory avoidance task','Shock in shock-probe burial task','RS')
psyc <- c('EPM','FST','NE','OFT')

df_filtered$tStressorTypeGroup[df_filtered$tStressorType %in% phys] <- 'physical'
df_filtered$tStressorTypeGroup[df_filtered$tStressorType %in% psyc] <- 'psychological'
df_filtered$tStressorTypeGroup[df_filtered$tStressorType == 'NA'] <- 'none'

# estimates for each subgroup
subgroup_phys <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ areaLevel2 -1,
                        subset = (tStressorTypeGroup=='physical'),
                        slab = paste(authors, year, sep=", "))
subgroup_psyc <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ areaLevel2 -1,
                        subset = (tStressorTypeGroup=='psychological'),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_type <- data.frame(estimate = c(coef(subgroup_phys), coef(subgroup_psyc)), 
                              stderror = c(subgroup_phys$se, subgroup_psyc$se),
                              meta = c(rep("physical", length(subgroup_phys$se)),
                                       rep("psychological", length(subgroup_psyc$se))), 
                              tau2 = c(rep(subgroup_phys$tau2, length(subgroup_phys$se)),
                                       rep(subgroup_psyc$tau2, length(subgroup_psyc$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_type, digits=3)


# Exploratory: Novelty ----------------------------------------------------
# estimates for each subgroup
subgroup_novel <- rma.mv(yi, vi, 
                         random = list(~1 | each, ~1 | nest),
                         method = "REML",
                         data = df_filtered,
                         mods = ~ areaLevel2 -1,
                         subset = (tNovel == 1),
                         slab = paste(authors, year, sep=", "))
subgroup_old <- rma.mv(yi, vi, 
                       random = list(~1 | each, ~1 | nest),
                       method = "REML",
                       data = df_filtered,
                       mods = ~ areaLevel2 -1,
                       subset = (tNovel == 0),
                       slab = paste(authors, year, sep=", "))

# summarize together
comparison_novelty <- data.frame(estimate = c(coef(subgroup_old), coef(subgroup_novel)), 
                                 stderror = c(subgroup_old$se, subgroup_novel$se),
                                 meta = c(rep("old", length(subgroup_old$se)),
                                          rep("novel", length(subgroup_novel$se))), 
                                 tau2 = c(rep(subgroup_old$tau2, length(subgroup_old$se)),
                                          rep(subgroup_novel$tau2, length(subgroup_novel$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_novelty, digits=3)


# Exploratory: Time -------------------------------------------------------
metaregression_time <- rma.mv(yi, vi, 
                              random = list(~1 | each, ~1 | nest),
                              method = "REML",
                              data = df_filtered,
                              mods = ~ tWaitPeriod -1,
                              subset = (tAcuteStressor == 1 & outMeasure == 'M'),
                              slab = paste(authors, year, sep=", "))
metaregression_time

### OR ###

# estimates for each subgroup
subgroup_early <- rma.mv(yi, vi, 
                         random = list(~1 | each, ~1 | nest),
                         method = "REML",
                         data = df_filtered,
                         mods = ~ areaLevel2 -1,
                         subset = (tAcuteStressor == 1 & outMeasure == 'P'& tWaitPeriod <= 60),
                         slab = paste(authors, year, sep=", "))
subgroup_late <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ areaLevel2 -1,
                        subset = (tAcuteStressor == 1 & outMeasure == 'P' & tWaitPeriod > 60),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_time <- data.frame(estimate = c(coef(subgroup_early), coef(subgroup_late)), 
                              stderror = c(subgroup_early$se, subgroup_late$se),
                              meta = c(rep("early", length(subgroup_early$se)),
                                       rep("late", length(subgroup_late$se))), 
                              tau2 = c(rep(subgroup_early$tau2, length(subgroup_early$se)),
                                       rep(subgroup_late$tau2, length(subgroup_late$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_time, digits=3)


# Exploratory: Out Measure ------------------------------------------------
# estimates for each subgroup
subgroup_prot <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ tAcuteStressor:areaLevel2 -1,
                        subset = (outMeasure == 'P'),
                        slab = paste(authors, year, sep=", "))
subgroup_mrna <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ tAcuteStressor:areaLevel2 -1,
                        subset = (outMeasure == 'M'),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_measure <- data.frame(estimate = c(coef(subgroup_prot), coef(subgroup_mrna)), 
                                 stderror = c(subgroup_prot$se, subgroup_mrna$se),
                                 meta = c(rep("protein", length(subgroup_prot$se)),
                                          rep("mRNA", length(subgroup_mrna$se))), 
                                 tau2 = c(rep(subgroup_prot$tau2, length(subgroup_prot$se)),
                                          rep(subgroup_mrna$tau2, length(subgroup_mrna$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_measure, digits=3)


### END ###