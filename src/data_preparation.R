## Script: "Overall data preparation"

## project: mIEG

## Input: excel data from osf (or downloaded)
## Output:
    ## 1) preprocessed data in data/temp/ folder for report "descriptives"
    ## 2) preprocessed data in data/processed/ folder for meta-analysis
    ## 3) risk of bias assessment data

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

# Meta-analytic filter in the end. 
## @ HEIKE: does not work
df_clean <- df %>%

  # select variables of interest
  dplyr::select(
    ID, nest, each, #id vars
    authors, year, #for labs
    species, strain, origin, sex, #animal info
    model, mCage, mTimeStart, mTimeLength, mHoursTotal, #model info
    hit2, #second hit
    tAcuteStressor, tStressorType, tNovel, tWaitPeriod, #acute stressor info
    iegName, outMeasure, outTechnique, outUnit, #outcome measures
    areaLevel1, areaLevel2, #brain areas
    nC, avgC, avgCtype, varC, varCtype, #control data
    nE, avgE, avgEtype, varE, varEtype, #experimental data
    meta, # inclusion in meta-analysis
    dataOrigin, cohensD, bias, blindRand, sigEffect
  ) %>%

  # preprocess variables
  mutate(

    # make necessary vars numeric
    across(c(varC, nC, avgC, varE, nE, avgE, cohensD, tWaitPeriod), as.numeric),

    # authors for citation
    authors = gsub("^(.*?),.*", "\\1 et al.", authors),

    # experiment ID
    exp_ID = paste(ID, nest, sep = "_"),

    # change variance type to sd >> here all either SEM or not specified (NS). NS considered sem - can be checked in varCtype and varEtype
    sdC = ifelse(is.na(varC) | is.na(nC), NA, varC * sqrt(nC)),
    sdE = ifelse(is.na(varE) | is.na(nE), NA, varE * sqrt(nE)),

    # update var factor type
    across(c(areaLevel1, areaLevel2, hit2, tStressorType, tAcuteStressor), as.factor)

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

#Add cohens d from stats papers only & impute var.
df_clean$yi[is.na(df_clean$yi)] <- df_clean$cohensD[is.na(df_clean$yi)]
df_clean$vi[is.na(df_clean$vi)] <- median(df_clean$vi, na.rm = T)

# blinding
#blind_me <- sample(1:nrow(df_clean), sample(1:nrow(df_clean), 1), replace = TRUE)
#df_clean$yi[blind_me] <- df_clean$yi[blind_me] * -1


# Save 1) ----------------------------------------------------------------
saveRDS(df_clean, paste0(temp, "df_report_datapreparation.RDS"))


# Save 2) -----------------------------------------------------------------
df_meta <-
  df_clean %>%
  filter(
    sex %in% c("M", "F"),
    iegName == "cFos",
    meta == 1,
    areaLevel2 %in% c('PFC','TH','HPF','HYP','AMY', 'MB')
  ) %>%
  ungroup() %>% droplevels()


saveRDS(df_meta, paste0(processed, "meta.RDS"))



# RoB ---------------------------------------------------------------------

rob <- df_list[["Outcomes"]] %>% 
  
  # get information
  dplyr::select(each, starts_with("RB_")) %>% 
  mutate(each = str_replace_all(each, "-.*", "")) %>% 
  unique() %>% 
  filter(each %in% unique(df$ID)) %>%
  
  # transform to long
  pivot_longer(cols = starts_with("RB_"), names_to = "type_bias", values_to = "value_bias") 

saveRDS(rob, paste0(processed, "rob.RDS"))

