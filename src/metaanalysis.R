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
df_filtered <- df_filtered[-which(inf$inf$inf == '*'),] %>%
  mutate(
    tAcuteStressor = ifelse(tStressorType == "NA", 0, 1), 
    stressor_intensity = case_when(
      tStressorType %in% c('OFT','DLB','NE','EPM') ~ "mild", # @ heike --> shouldnt this be mild? >> below you had 1
      tStressorType %in% c('CRD','FST','MWM','RS',
                           'FS in inhibitory avoidance task',
                           'Shock in shock-probe burial task') ~"severe", 
      TRUE ~ "mild"
    )
  )


# Model -------------------------------------------------------------------
mod <- rma.mv(yi, vi, 
              random = list(~1 | each, ~1 | nest),
              method = "REML",
              data = df_filtered,
              
              mods = ~ areaLevel2:factor(tAcuteStressor):hit2 -1, 
           #   mods = ~ tAcuteStressor:areaLevel2 -1, 
              slab = paste(authors, year, sep=", "))


# Main analysis -----------------------------------------------------------
## Create contrasts per RQ
interactions <- names(data.frame(mod$X))

# RQ1: effects of ELA on presence acute stressor
acutestressor <- ifelse(str_detect(interactions, "tAcuteStressor.1"), 1, 0)
contr_acutestressor <- rbind(
  ifelse(acutestressor == 1, 1/sum(acutestressor == 1), 0), # acute vs 0
  ifelse(acutestressor == 0, 1/sum(acutestressor == 0), 0) # baseline vs 0
)

anova(mod, L = contr_acutestressor) 
main_res <- summary(glht(mod, linfct = contr_acutestressor, df=df.residual(mod)), test=adjusted('holm'))


# RQ2: effects of ELA on brain area
brainarea <- gsub('(.*)\\.\\w+', '\\1', interactions)

contr_brainarea <- NULL
for (area in unique(brainarea)) {
  x <- ifelse(brainarea == area, 1/sum(brainarea == area), 0)
  contr_brainarea <- rbind(contr_brainarea, x)
}

anova(mod, L = contr_brainarea) 
res_ba <- summary(glht(mod, linfct = contr_brainarea, df=df.residual(mod)), test = adjusted('holm'))

# Subgroup: hit2 ----------------------------------------------------------
# estimates for each subgroup
subgroup_hit0 <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ factor(tAcuteStressor):areaLevel2 -1,
                        subset = (hit2==0),
                        slab = paste(authors, year, sep=", "))
subgroup_hit1 <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ factor(tAcuteStressor):areaLevel2 -1,
                        subset = (hit2==1),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_hit <- data.frame(estimate = c(coef(subgroup_hit0), coef(subgroup_hit1)), 
                             stderror = c(subgroup_hit0$se, subgroup_hit1$se),
                             meta = c(rep("hit0", length(subgroup_hit0$se)),
                                      rep("hit1", length(subgroup_hit1$se))), 
                             stressor = c('BL','AS','BL','AS','BL','AS','AS','AS',
                                          'BL','AS','BL','AS','BL','AS','BL','AS','BL','AS'),
                             tau2 = c(rep(subgroup_hit0$tau2, length(subgroup_hit0$se)),
                                      rep(subgroup_hit1$tau2, length(subgroup_hit1$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_hit, digits=3)

# Wald test: to test the difference between the two
# with(comparison_hit, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))

# meta-analyze interaction 2nd hit by stressor
hit_mod <- rma(estimate, sei=stderror, mods = ~ meta:stressor -1, method="FE", data=comparison_hit, digits=3)



# Exploratory: severity ----------------------------------------------------------
df_filtered %>% 
  filter(tAcuteStressor == "1") %>%
  group_by(hit2, stressor_intensity) %>% 
  count()

severity_res <- rma.mv(yi, vi, 
                  random = list(~1 | each, ~1 | nest),
                  method = "REML",
                  data = df_filtered,
                  mods = ~stressor_intensity:hit2 -1,
                  subset = (tAcuteStressor=="1"),
                  slab = paste(authors, year, sep=", "))

# Visualizations ----------------------------------------------------------
## baseline & acute stress
df_fig_main <- data.frame(
  type = c("Acute stress", "Baseline"), 
  facet = rep("Main analysis", 2),
  g = main_res$test$coefficients, 
  sem = main_res$test$sigma,
  pval = main_res$test$pvalues, 
  n = df_filtered %>% group_by(tAcuteStressor) %>% count() %>% arrange(-tAcuteStressor) %>% pull(n)
)

g_main <- df_fig_main %>%
  mutate(type = factor(type, levels = c("Baseline", "Acute stress"))) %>%
  ggplot(aes(type, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  
  # beautiful
  ylim(c(-0.2,1.1)) + 
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  facet_grid(~facet) + 
  geom_text(x = "Baseline", y = max(df_fig_main$g) + 0.4, label = "*", size = 20) + 
  geom_text(aes(y = 0.05, label = paste("n = ", n), size = 20)) + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        legend.position = "none")

## second hit
df_fig_hit <- data.frame(
  hit = c("no hit", "2nd hit", "no hit", "2nd hit"),
  type = c("Acute stress", "Acute stress", "Baseline", "Baseline"), 
  g = hit_mod$b, 
  sem = hit_mod$se,
  pval = hit_mod$pval,
  n = df_filtered %>% group_by(tAcuteStressor, hit2) %>% count() %>% arrange(-tAcuteStressor) %>% pull(n)
)

g_hit <- df_fig_hit %>%
  mutate(type = factor(type, levels = c("Baseline", "Acute stress"))) %>%
  mutate(hit = factor(hit, levels = c("no hit", "2nd hit"))) %>%
    mutate(text = ifelse(pval < 0.05, "*", "")) %>%
    
  ggplot(aes(hit, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  facet_grid(~type) + 
  
  # beautiful
  ylim(c(-0.3,1.3)) + 
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  geom_text(aes(label = text, y = g+0.3),  size = 20) + 
  geom_text(aes(y = 0.05, label = paste("n = ", n), size = 20)) + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = "none")


ggarrange(
  g_main, g_hit, nrow = 1, widths = c(1,2),
  labels = c("A)", "B)")
)

## brain area
df_fig_ba <- data.frame(
  ba = c("AMY", "HPF", "HYP", "AMY", "HPF", "HYP", "PFC", "TH", "PFC", "TH"),
  type = c("Baseline", "Baseline", "Baseline", "Acute stress", "Acute stress", 
           "Acute stress", "Acute stress","Acute stress", "Baseline", "Baseline"), 
  g = res_ba$test$coefficients, 
  sem = res_ba$test$sigma,
  pval = res_ba$test$pvalues
) %>% 
  left_join(
    df_filtered %>% 
      group_by(tAcuteStressor, areaLevel2) %>% 
      count() %>% 
      rename(ba = areaLevel2, type = tAcuteStressor) %>% 
      mutate(type = ifelse(type == 1, "Acute stress", "Baseline"))
  )

g_ba <- df_fig_ba %>%
  mutate(
    type = factor(type, levels = c("Baseline", "Acute stress")),
    ba = case_when(
      ba == "AMY" ~ "Amygdala", 
      ba == "PFC" ~ "Prefrontal cx", 
      ba == "HYP" ~ "Hypothalamus", 
      ba == "HPF" ~ "Hippocampus", 
      ba == "TH" ~ "Thalamus"
    )) %>%
  ggplot(aes(ba, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  facet_grid(~type) + 
  
  # beautiful
  ylim(c(-0.5,1.1)) + 
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  geom_text(aes(y = -0.5, label = paste("n = ", n), size = 20)) + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = "none")

g_ba


# Cartoon discussion ------------------------------------------------------
discussion_a <- df_fig_main %>%
  mutate(type = factor(type, levels = c("Baseline", "Acute stress"))) %>%
  ggplot(aes(type, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  
  # beautiful
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  facet_grid(~type, scales = "free_x") + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        legend.position = "none")

df_discussion <- data.frame(
  type = rep(c("Baseline", "Acute stress"), each = 4),
  group = rep(c("Control", "ELA"), 4),
  hit = rep(rep(c("no hit", "2nd hit"), each = 2),2),
  g = c(1,4,3,5,9,8,10, 13),
  sem = rep(0.5, 8)
)


df_discussion %>%
  mutate(type = factor(type, levels = c("Baseline", "Acute stress"))) %>%
  mutate(hit = factor(hit, levels = c("no hit", "2nd hit"))) %>%
  mutate(
    label = ifelse(type == "Acute stress" & hit == "no hit", "", "*"),
    y_label = ifelse(type == "Acute stress", 14, 8)) %>%
  
  mutate(each_group = paste(group, hit, sep = "_")) %>%
  ggplot(aes(hit, g, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.9), colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), position =position_dodge(width=0.9), width = 0.2) +
  
  # beautiful
  theme_bw() + 
  labs(fill = "Exp group:",
       y = "cFos expression",
       x = "") +
  
  scale_fill_manual(values = c("white", "grey")) +
  facet_grid(~type, scales = "free_x") + 
  geom_text(aes(y = y_label, label = label), size = 20) + 
  theme(text = element_text(size = 20), 
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

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
rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_sens, digits=3)


# Exploratory: Time -------------------------------------------------------
metaregression_time <- rma.mv(yi, vi, 
                              random = list(~1 | each, ~1 | nest),
                              method = "REML",
                              data = df_filtered,
                              mods = ~ tWaitPeriod -1,
                              subset = (tAcuteStressor == 1 & outMeasure == 'M'),
                              slab = paste(authors, year, sep=", "))
metaregression_time

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

# Publication bias --------------------------------------------------------
funnel(mod, xlab = "Hedge's g")
mod_uni <- rma(yi, vi, 
              method = "REML",
              data = df_filtered,
              mods = ~ tAcuteStressor:areaLevel2 -1, 
              slab = paste(authors, year, sep=", "))
regtest(mod_uni)

high_es <- c(26452320, 24513388)


# Forest plot -------------------------------------------------------------

forest(mod)
