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


# Model M -----------------------------------------------------------------
df_filtered <- df %>% filter(sex == "M", species =='rat', areaLevel2 != 'PFC')
mod <- rma.mv(yi, vi, 
              random = list(~1 | each, ~1 | nest),
              method = "REML",
              data = df_filtered,
              mods = ~ tAcuteStressor:areaLevel2 -1, # ideally my model would also include hit2 and species since those are my subgroups
              slab = paste(authors, year, sep=", "))

## http://www.metafor-project.org/doku.php/tips:multiple_factors_interactions


# Main analysis ------------------------------------------------
## Create contrasts per RQ
interactions <- names(data.frame(mod$X))

# RQ1: effects of ELA on presence acute stressor
acutestressor <- ifelse(str_detect(interactions, "tAcuteStressor1"), 1, 0)
contr_acutestressor <- rbind(
  ifelse(acutestressor == 1, 1/sum(acutestressor == 1), 0), # acute vs 0
  ifelse(acutestressor == 0, 1/sum(acutestressor == 0), 0) # baseline vs 0
)

anova(mod, L = contr_acutestressor) # add multiple testing correction?
summary(glht(mod, linfct = contr_acutestressor, df=df.residual(mod)), test=adjusted('holm'))


# RQ2: effects of ELA on brain area
brainarea <- sub(".*2", "", interactions)
contr_brainarea <- NULL
for (area in unique(brainarea)) {
  x <- ifelse(brainarea == area, 1/sum(brainarea == area), 0)
  contr_brainarea <- rbind(contr_brainarea, x)
}

anova(mod, L = contr_brainarea) # add multiple testing correction?
summary(glht(mod, linfct = contr_brainarea, df=df.residual(mod)), test = adjusted('holm'))


# Subgroup: hit2 ----------------------------------------------------------
## tutorial from: http://www.metafor-project.org/doku.php/tips:comp_two_independent_estimates?s[]=subset

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
                             tau2 = c(rep(subgroup_hit0$tau2, length(subgroup_hit0$se)),
                                      rep(subgroup_hit1$tau2, length(subgroup_hit1$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_hit, digits=3)

# Wald test: to test the difference between the two
with(comparison_hit, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))


# Model F -----------------------------------------------------------------
df_filtered <- df %>% filter(sex == "F", areaLevel2 != 'PFC', outMeasure == 'P',species='rat')
mod <- rma.mv(yi, vi, 
              random = list(~1 | each, ~1 | nest),
              method = "REML",
              data = df_filtered,
              mods = ~ tAcuteStressor:areaLevel2 -1, 
              slab = paste(authors, year, sep=", "))


# Main analysis ------------------------------------------------
## Create contrasts per RQ
interactions <- names(data.frame(mod$X))

# RQ1: effects of ELA on presence acute stressor
acutestressor <- ifelse(str_detect(interactions, "tAcuteStressor1"), 1, 0)
contr_acutestressor <- rbind(
  ifelse(acutestressor == 1, 1/sum(acutestressor == 1), 0), # acute vs 0
  ifelse(acutestressor == 0, 1/sum(acutestressor == 0), 0) # baseline vs 0
)

anova(mod, L = contr_acutestressor) # add multiple testing correction?
summary(glht(mod, linfct = contr_acutestressor, df=df.residual(mod)), test=adjusted('holm'))


# RQ2: effects of ELA on brain area
brainarea <- sub(".*2", "", interactions)
contr_brainarea <- NULL
for (area in unique(brainarea)) {
  x <- ifelse(brainarea == area, 1/sum(brainarea == area), 0)
  contr_brainarea <- rbind(contr_brainarea, x)
}

anova(mod, L = contr_brainarea) # add multiple testing correction?
summary(glht(mod, linfct = contr_brainarea, df=df.residual(mod)), test = adjusted('holm'))


# Subgroup: hit2 ----------------------------------------------------------
## tutorial from: http://www.metafor-project.org/doku.php/tips:comp_two_independent_estimates?s[]=subset

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
                             tau2 = c(rep(subgroup_hit0$tau2, length(subgroup_hit0$se)),
                                      rep(subgroup_hit1$tau2, length(subgroup_hit1$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_hit, digits=3)

# Wald test: to test the difference between the two
with(comparison_hit, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))


