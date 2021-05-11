## Script: "Overall data preparation"

## project: mIEG

## Input: excel data from osf (or downloaded)
## Output:
## 1) preprocessed data in data/temp/ folder for report "descriptives"
## 2) preprocessed data in data/processed/ folder for meta-analysis

## Author: Valeria Bonapersona & Heike Schuler
## contact: v.bonapersona-2 (at) umcutrecht.nl

'
  Script: "Risk of bias assessment"
  Project: mIEG
  
  Input: risk of bias assessment table
  Output: figure
'


# Environment preparation -------------------------------------------------
source("src/utilities.R")

rob <- readRDS(paste0(processed, "rob.RDS"))



# Visualization -----------------------------------------------------------
g <- rob %>%
  mutate(
    value_bias = case_when(
      value_bias %in% c("C", "L") ~ "Low",
      value_bias %in% c("H") ~ "High", 
      T ~ "Unclear"),
    value_bias = factor(value_bias, levels = c("Low", "Unclear", "High")), 
    
    type_bias = case_when(
      type_bias %in% c("RB_seqGeneration") ~ "Sequence generation",
      type_bias %in% c("RB_outBlind") ~ "Blind outcome assessment",
      type_bias %in% c("RB_outAss") ~ "Random outcome assessment",
      type_bias %in% c("RB_incData") ~ "Incomplete outcome data",
      type_bias %in% c("RB_housing") ~ "Random housing",
      type_bias %in% c("RB_control") ~ "Appropriate control group",
      type_bias %in% c("RB_blindExp") ~ "Blind exp performance",
      type_bias %in% c("RB_baseline") ~ "Baseline characteristics",
      type_bias %in% c("RB_allocation") ~ "Allocation concealment"
    ), 
    type_bias = factor(type_bias, levels = c(
      "Appropriate control group", "Baseline characteristics", "Sequence generation",
      "Random housing", "Allocation concealment", "Random outcome assessment",
      "Blind exp performance", "Incomplete outcome data", "Blind outcome assessment"
      ))
    ) %>%
  ggplot(aes(type_bias, fill = value_bias)) + 
  geom_bar(position = "fill") + 
  coord_flip() + 
  
  # beautify
  labs(x = "",
       y = "Percentage of studies", 
       fill = "The Risk of Bias was:") +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + 
  scale_fill_manual(values = c("#32a852", "#d9c750", "#d42020")) + 
  theme(legend.position = "bottom", 
        text = element_text(size = 18))

g

saveRDS(g, paste0(figs_rds, "rob.RDS"))
