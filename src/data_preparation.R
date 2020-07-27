## Script: "Relative paths"

## project: mIEG

## Author: Valeria Bonapersona & Heike Schuler
## contact: v.bonapersona-2 (at) umcutrecht.nl



# Environment preparation -------------------------------------------------
source("src/utilities.R")

# import data
df_list <- read_sheets(paste0(raw, "outcomes.xlsx"), 
                       c("Publications", "Experiments", "Outcomes"))


# Create dataframe --------------------------------------------------------
# Check data IDs for coherence.
identical(unique(df_list$Publications$ID), unique(df_list$Experiments$ID))
identical(unique(df_list$Experiments$nest), unique(df_list$Outcomes$nest))


# Combine sheets into one dataset from which to create sub-datasets.
df <- data.frame(df_list$Publications[match(df_list$Experiments$ID, df_list$Publications$ID),], 
                 df_list$Experiments)
df <- data.frame(df[match(df_list$Outcomes$nest, df$nest),], df_list$Outcomes)


# Checks ------------------------------------------------------------------
levels(as.factor(c(as.character(df$varEtype),
                   as.character(df$varCtype)))) # if NA, assumed SEM


# Filtering for meta-analysis ---------------------------------------------
# rename CA1, CA3, DG consistently (remove d, v)
df$areaLevel1[df$areaLevel1 %in% c('vCA1','vCA2','vCA3','vDG','dCA1','dCA2','dCA3','dDG')] <-
  substring(df$areaLevel1[df$areaLevel1 %in% c('vCA1','vCA2','vCA3','vDG','dCA1','dCA2','dCA3','dDG')],2)

# SEX NOT SELECTED HERE!
df_clean <- df %>%
  
  # studies with meta-analytic information 
  filter(
 #   meta == 1, # @ heike - how did you create this variable?, @ valeria - do it on all?
    
    areaLevel2 %in% c('PFC','TH','HPF','HYP','AMY'), # specific regions - due to frequency?
  ) %>%
  
  # select variables of interest
  select(
    ID, nest, each, #id vars
    authors, year, #for labs
    species, strain, origin, sex, #animal info
    model, mCage, mTimeStart, mTimeLength, mHoursTotal, #model info
    hit2, #second hit
    tAcuteStressor, tStressorType, tNovel, tWaitPeriod, #acute stressor info
    outMeasure, outTechnique, outUnit, #outcome measures
    areaLevel1, areaLevel2, #brain areas
    nC, avgC, avgCtype, varC, varCtype, #control data
    nE, avgE, avgEtype, varE, varEtype, #experimental data
    dataOrigin, cohensD, bias, blindRand, sigEffect
  ) %>%
  
  # preprocess variables
  mutate(
    
    # make necessary vars numeric
    across(c(varC, nC, avgC, varE, nE, avgE), as.numeric), # CANCEL WARNING MESSAGE NAs!
    
    # authors for citation
    authors = gsub("^(.*?),.*", "\\1 et al.", authors),
    
    # change variance type to sd >> here all either SEM or not specified (NS). NS considered sem - can be checked in varCtype and varEtype
    sdC = ifelse(is.na(varC) | is.na(nC), NA, varC * sqrt(nC)),
    sdE = ifelse(is.na(varE) | is.na(nE), NA, varE * sqrt(nE)), 
    
  ) %>%
  
  # cleaning
  ungroup() %>% droplevels()
  

# perform t-test for systematic review
df_clean$p_val <- tsum.test(df_clean$avgC, df_clean$sdC,df_clean$nC, 
                           df_clean$avgE, df_clean$sdE, df_clean$nE)$p.value

# estimate effect sizes for meta-analysis
df_clean <- escalc(m1i = avgE, sd1i = sdE, n1i = nE, 
                  m2i = avgC, sd2i = sdC, n2i = nC, 
                  measure = "SMD",
                  data = df_clean)

# blinding
blind_me <- sample(1:nrow(df_clean), sample(1:nrow(df_clean), 1), replace = TRUE)
df_clean$yi[blind_me] <- df_clean$yi[blind_me] * -1


# Save ----------------------------------------------------------------
saveRDS(df_clean, paste0(processed, "data_onefile.RDS"))
