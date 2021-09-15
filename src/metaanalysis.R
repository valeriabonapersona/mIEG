## Script: "Meta-analysis cFos"

## project: mIEG

## Input: preprocessed data from data_preparation.R adjusted by report_descriptives.rmd
## Output: meta-analytic results

## Author: Valeria Bonapersona & Heike Schuler
## contact: v.bonapersona-2 (at) umcutrecht.nl


# Environment preparation -------------------------------------------------
source("src/utilities.R")
df <- readRDS(paste0(processed, "meta.RDS"))

# sys review
length(unique(df$ID))

df %>% 
  mutate(sig = ifelse(p_val <= 0.05, "yes", "no")) %>% 
  group_by(sig) %>% 
  count()

df %>% 
  filter(sigEffect == 1) %>% 
  mutate(pos = ifelse(yi > 0, "yes", "no")) %>% 
  group_by(pos) %>% 
  count()

# Outliers ----------------------------------------------------------------
# Create z-score variable.
df$zi <- (df$yi - mean(df$yi, na.rm = T)) / 
  sd(df$yi, na.rm = T)

df_filtered <- df %>% filter(sex == "M")

#Check for influential outliers and remove.
res <- rma(yi, vi, measure="ZCOR", data=df_filtered)
inf <- influence(res)
outliers <- df_filtered[which(inf$inf$inf == '*'),]
df_filtered <- df_filtered[-which(inf$inf$inf == '*'),] %>%
  mutate(
    tAcuteStressor = ifelse(tStressorType == "NA", 0, 1), # change 0 and 1 to acute and rest
    stressor_intensity = case_when(
      tStressorType %in% c('OFT','DLB','NE','EPM') ~ "mild",
      tStressorType %in% c('CRD','FST','MWM','RS',
                           'FS in inhibitory avoidance task',
                           'Shock in shock-probe burial task') ~"severe", 
      tAcuteStressor == 0 ~ "not_applicable",
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

res_anova_ba <- anova(mod, L = contr_brainarea) 
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
                             stressor = c('BL','AS','BL','AS','BL','AS','BL','AS','BL','AS','AS',
                                          rep(c('BL', 'AS'), 6)),
                             # stressor = c('BL','AS','BL','AS','BL','AS','AS','AS',
                             #              'BL','AS','BL','AS','BL','AS','BL','AS','BL','AS'),
                             tau2 = c(rep(subgroup_hit0$tau2, length(subgroup_hit0$se)),
                                      rep(subgroup_hit1$tau2, length(subgroup_hit1$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ meta -1, method="FE", data=comparison_hit, digits=3)

# Wald test: to test the difference between the two
# with(comparison_hit, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))

# meta-analyze interaction 2nd hit by stressor
hit_mod <- rma(estimate, sei=stderror, mods = ~ meta:stressor -1, method="FE", data=comparison_hit, digits=3)

p.adjust(hit_mod$pval, method = "holm", n = length(hit_mod$pval))

 

# Sensitivity: species -----------------------------------------------------
subgroup_rat <- rma.mv(yi, vi, 
                        random = list(~1 | each, ~1 | nest),
                        method = "REML",
                        data = df_filtered,
                        mods = ~ factor(tAcuteStressor):areaLevel2 -1,
                        subset = (species=="rat"),
                        slab = paste(authors, year, sep=", "))

# summarize together
comparison_species <- data.frame(estimate = c(coef(mod), coef(subgroup_rat)), 
                              stderror = c(mod$se, subgroup_rat$se),
                              sens = c(rep("all", length(mod$se)),
                                       rep("rat_only", length(subgroup_rat$se))), 
                              tau2 = c(rep(mod$tau2, length(mod$se)),
                                       rep(subgroup_rat$tau2, length(subgroup_rat$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_species, digits=3)


# Sensitivity: species only acute -----------------------------------------------------
subgroup_rat_acute <- rma.mv(yi, vi, 
                       random = list(~1 | each, ~1 | nest),
                       method = "REML",
                       data = df_filtered %>% filter(tAcuteStressor==1),
                       mods = ~ areaLevel2 -1,
                       subset = (species=="rat"),
                       slab = paste(authors, year, sep=", "))
acute <- rma.mv(yi, vi, 
               random = list(~1 | each, ~1 | nest),
               method = "REML",
               data = df_filtered %>% filter(tAcuteStressor==1),
               mods = ~ areaLevel2 -1,
               slab = paste(authors, year, sep=", "))

# summarize together
comparison_species_acute <- data.frame(estimate = c(coef(acute), coef(subgroup_rat_acute)), 
                                 stderror = c(acute$se, subgroup_rat_acute$se),
                                 sens = c(rep("all", length(acute$se)),
                                          rep("rat_only", length(subgroup_rat_acute$se))), 
                                 tau2 = c(rep(acute$tau2, length(acute$se)),
                                          rep(subgroup_rat_acute$tau2, length(subgroup_rat_acute$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_species_acute, digits=3)

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
  mutate(type = ifelse(type == "Baseline", "At rest", "Acute stress")) %>%
  mutate(type = factor(type, levels = c("At rest", "Acute stress"))) %>%  ggplot(aes(type, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  
  # beautiful
  ylim(c(-0.2,0.8)) + 
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  facet_grid(~facet) + 
  geom_text(x = "At rest", y = max(df_fig_main$g) + 0.2, label = "*", size = 20) + 
  geom_text(aes(y = 0.05, label = paste("n = ", n), size = 20)) + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        legend.position = "none")

## second hit
df_fig_hit <- data.frame(
  hit = rep(c("No additional hits", "Additional hits"),2),
  type = c("Acute stress", "Acute stress", "Baseline", "Baseline"), 
  g = hit_mod$b, 
  sem = hit_mod$se,
  pval = hit_mod$pval,
  n = df_filtered %>% group_by(tAcuteStressor, hit2) %>% count() %>% arrange(-tAcuteStressor) %>% pull(n)
)

g_hit <- df_fig_hit %>%
  mutate(type = ifelse(type == "Baseline", "At rest", "Acute stress")) %>%
  mutate(type = factor(type, levels = c("At rest", "Acute stress"))) %>%
  mutate(hit = factor(hit, levels = c("No additional hits", "Additional hits"))) %>%
    mutate(text = ifelse(pval < 0.05, "*", "")) %>%
    
  ggplot(aes(hit, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  facet_grid(~type) + 
  
  # beautiful
  ylim(c(-0.5,1)) + 
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

ggsave(paste0(figs_svg, "main_res.eps"), width = 9, height = 6)

## brain area
df_fig_ba <- data.frame(
  ba = str_remove_all(res_anova_ba$hyp[,1], "\\:.*") %>% str_remove_all(".*2"),
  type = str_remove_all(res_anova_ba$hyp[,1], ".*\\)") %>% str_remove_all("\\:.*"), 
  g = res_ba$test$coefficients, 
  sem = res_ba$test$sigma,
  pval = res_ba$test$pvalues
) %>% 
  mutate(type = ifelse(type == "0", "Baseline", "Acute stress")) %>%
  left_join(
    df_filtered %>% 
      group_by(tAcuteStressor, areaLevel2) %>% 
      count() %>% 
      rename(ba = areaLevel2, type = tAcuteStressor) %>% 
      mutate(type = ifelse(type == 1, "Acute stress", "Baseline"))
  )

g_ba <- df_fig_ba %>%
  mutate(
    type = ifelse(type == "Baseline", "At rest", type),
    type = factor(type, levels = c("At rest", "Acute stress")),
    ba = case_when(
      ba == "AMY" ~ "Amygdala", 
      ba == "PFC" ~ "Prefrontal cx", 
      ba == "HYP" ~ "Hypothalamus", 
      ba == "HPF" ~ "Hippocampus", 
      ba == "TH" ~ "Thalamus", 
      ba == "MB" ~ "Midbrain"
    )) %>%
  ggplot(aes(ba, g)) + 
  geom_bar(stat = "identity", fill = "white", colour = "black") + 
  geom_errorbar(aes(ymin = g-sem, ymax = g+sem), width = 0.2) +
  facet_grid(~type) + 
  
  # beautiful
  ylim(c(-0.6,1.1)) + 
  theme_bw() + 
  labs(x="", y = "Hedge's g") + 
  geom_text(aes(y = -0.6, label = paste("n = ", n), size = 20)) + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = "none")

g_ba

ggsave(paste0(figs_svg, "ba.eps"), width = 10, height = 5)


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
  g = c(1,4,3,5,9,7,11, 15),
  sem = c(rep(0.5, 4), rep(1.5,4))
)


df_discussion %>%
  mutate(type = ifelse(type == "Baseline", "At rest", type)) %>%
  mutate(type = factor(type, levels = c("At rest", "Acute stress"))) %>%
  mutate(hit = factor(hit, levels = c("no hit", "2nd hit"))) %>%
  mutate(
    label = ifelse(type == "Acute stress" & hit == "no hit", "", "*"),
    y_label = ifelse(type == "Acute stress", 16, 8)) %>%
  
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

ggsave(paste0(figs_svg, "discussion.eps"), width = 10, height = 6)

# Sensitivity: RNA/protein ----------------------------------------------------
df_filtered %>% 
  group_by(outTechnique) %>% 
  count()

# model compared to complete model
mod_sens <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = df_filtered,
                   subset = (outTechnique=='IHC'),
                   mods = ~ areaLevel2:factor(tAcuteStressor):hit2 -1, 
                   slab = paste(authors, year, sep=", "))

# summarize together
comparison_sens <- data.frame(estimate = c(coef(mod), coef(mod_sens)), 
                              stderror = c(mod$se, mod_sens$se),
                              sens = c(rep("all", length(mod$se)),
                                       rep("IHC_only", length(mod_sens$se))), 
                              tau2 = c(rep(mod$tau2, length(mod$se)),
                                       rep(mod_sens$tau2, length(mod_sens$se))))

# meta-analyze with fixed effect
rma(estimate, sei=stderror, mods = ~ sens -1, method="FE", data=comparison_sens, digits=3)

# Sensitivity: ELA model --------------------------------------------------
df_filtered %>% 
  group_by(model) %>% 
  count()

# model compared to complete model
mod_sens <- rma.mv(yi, vi, 
                   random = list(~1 | each, ~1 | nest),
                   method = "REML",
                   data = df_filtered,
                   subset = (model=='MS'),
                   mods = ~ areaLevel2:factor(tAcuteStressor):hit2 -1, 
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
## adapted from metafor tutorial
### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

df_filtered$stressor_intensity <- factor(df_filtered$stressor_intensity, levels=c('severe','mild','not_applicable'))
df_filtered$species[df_filtered$ID == '29063642'] = 'rat'
dat <- df_filtered %>% 
  mutate(
    acute_stressor = case_when(
      str_detect(tStressorType, "FS") ~ "FS", 
      str_detect(tStressorType, "Shock") ~ "Shock burial", 
      tStressorType == "OFT" ~ "OF",
      tStressorType == "NA" ~ "",
      T ~ as.character(tStressorType)
    ), 
    hit2 = ifelse(hit2=="0", "no", "yes"),
  ) %>%
  arrange(species, model, stressor_intensity, hit2, desc(yi)) %>% 
  mutate(my_order = row_number())
  
### set up forest plot (with 2x2 table counts added; the 'rows' argument is
### used to specify in which rows the outcomes will be plotted)
start = 3
break_small_size = 3
break_medium_size = 5
break_large_size = 7

# Mice
break_1 <- dat %>% filter(species == 'mice' & model == 'LBN') %>% nrow() + start -1
break_2 <- dat %>% filter(species == 'mice' & model == 'MS') %>% nrow() + break_1 + break_small_size -1

# Rats
break_3 <- dat %>% filter(species == 'rat' & model == 'LG') %>% nrow() + break_2 + break_large_size -1
break_4 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'severe' & hit2 == 'no') %>% nrow() + break_3 + break_large_size -1
break_5 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'severe' & hit2 == 'yes') %>% nrow() + break_4 + break_small_size -1
break_6 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'mild' & hit2 == 'no') %>% nrow() + break_5 + break_medium_size -1
break_7 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'mild' & hit2 == 'yes') %>% nrow() + break_6 + break_small_size -1
break_8 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'not_applicable' & hit2 == 'no') %>% nrow() + break_7 + break_medium_size -1
break_9 <- dat %>% filter(species == 'rat' & model == 'MS' & stressor_intensity == 'not_applicable' & hit2 == 'yes') %>% nrow() + break_8 + break_small_size -1

end <- break_9 + 8

mod <- rma.mv(yi, vi, 
              random = list(~1 | each, ~1 | nest),
              method = "REML",
              data = dat,
              
              mods = ~ areaLevel2:factor(tAcuteStressor):hit2 -1, 
              #   mods = ~ tAcuteStressor:areaLevel2 -1, 
              slab = paste(authors, year, sep=", "))

op <- par(cex=0.5, font=1)
forest(mod, xlim=c(-12.5, 6), at=c(-4, -2, 0, 2, 4), 
       #atransf=exp,
       ilab=cbind(dat$model, dat$hit2, as.character(dat$acute_stressor), as.character(dat$areaLevel2)),
       ilab.xpos=c(-9.5,-8,-6.5,-5), 
       cex=1.8, ylim=c(-1, end),
       order=dat$my_order, rows=c(start:break_1, #mouse LG mild nohit2
                                  (break_1+break_small_size):break_2, # mouse MS rest nohit2
                                  (break_2+break_large_size):break_3, # rats LG shock nohit2
                                  (break_3+break_large_size):break_4,
                                  (break_4+break_small_size):break_5,
                                  (break_5+break_medium_size):break_6,
                                  (break_6+break_small_size):break_7,
                                  (break_7+break_medium_size):break_8,
                                  (break_8+break_small_size):break_9),
                                 
       mlab=mlabfun("RE Model for All Studies", mod),
       addfit=FALSE,
       psize=1)

### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex=1, font=2)

### add additional column headings to the plot
text(c(-12.5), (end), c("Author(s) and Year"), adj=0)
text(c(-9.5,-8,-6.5,-5), (end), c("ELA model", "Additional hit", "Acute stress","Area"))
text(-12.5, c(break_2+4, break_3+4, break_9+4), c("Mice", "Rats, other ELA", "Rats, MS"), adj=0)


### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex=0.8, font=2)

### add additional column headings to the plot
text(-12.5, c(break_1+1.5, break_2+1.5), c("LBN in mice, mild stress, no additional hits", 
                                       "MS in mice, at rest, no additional hits"), adj=0)
text(-12.5, c(break_3+1.5), c("LG in rats, severe stress, no additional hits"),adj=0)
text(-12.5, c(break_4+1.5,break_5+1.5), c("Severe stress, additional hits",
                                      "Severe stress, no additional hits"),adj=0)
text(-12.5, c(break_6+1.5,break_7+1.5), c("Mild stress, additional hits",
                                          "Mild stress, no additional hits"),adj=0)
text(-12.5, c(break_8+1.5,break_9+1.5), c("At rest, additional hits",
                                          "At rest, no additional hits"),adj=0)


## save width = 1800, height = 2400; save with cex = 0.5 for error bars
