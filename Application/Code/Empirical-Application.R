################################################################################
######################### QP Empirical Application #############################
################################################################################

############################ Script Description ################################
#
# Author: Cameron
# 
# Date Created: 04/16/24
#
#
# Script Description:
#
#
# Last Updated: 07/31/2024 
#
#
# Notes:
#   To-Do:
#     + Clean up code comments (especially in data cleaning sections) 
# 
#   Next: 
#     + 
# 
#   Done: 
# 
#     + get weights & figure out how to incorporate them properly into the analysis 
#     + create rough draft of code to run analyses 
#     + 
#     + Note. pg 14 of 21600-User_guide.pdf provides weights 
#     + Multilevel Model sample weight example on pg 43 of 21600-User_guide.pdf
# 
################################################################################


# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  tidyverse, 
  ggplot2, 
  flextable, 
  stringr, 
  mice, 
  cobalt, 
  WeightIt, 
  boot, 
  utils, 
  lme4, 
  WeMix, 
  parallel
  # ggdag, 
  # dagitty, 
  # huxtable
)

# Load Functions 
source("Application/Functions/bootstrap_ci_paral_2.R")
source("Application/Functions/bootstrap_ci_re_paral_2.R")




# Import & Clean Data -----------------------------------------------------

## Wave I ------------------------------------------------------------------

# Import student data from Wave I of the Add Health study
# The Add Health study is a longitudinal study of adolescents in the United States
w1 <- readr::read_tsv(file = "Application/Data/ICPSR_21600/DS0001/21600-0001-Data.tsv")

# Select relevant variables from Wave I student data
# This includes demographic information, feelings scale (CES-D items), and self-esteem scale
w1 <- w1 %>% 
  select("AID":"SCH_YR", "S1":"S7", "S12", "S18", "S11", "S17", "PC22",
         "H1FS1":"H1FS19", # Feelings Scale (contains mostly CES-D items for depression assessment)
         "H1NF12B", "S44A14", "S44A18":"S44A29", 
         "H1PF30", "H1PF32":"H1PF34", "H1PF36") # Self-esteem scale items

# Rename key demographic variables for clarity
# Variable names are based on the In School Questionnaire Code Book Public Use Sample (pg. 697 of Questionnaire.pdf)
data <- w1 %>% 
  rename(age = S1, 
         sex = S2, 
         ethnicity = S4, # Indicates Hispanic or Spanish origin
         white = S6A, 
         black = S6B, 
         asian = S6C, 
         nativeAmerican = S6D, 
         raceOther = S6E)

# Data cleaning and scoring for self-esteem scale

# Replace invalid responses (6=refused; 8=don't know) with NA for self-esteem items
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]] <- 
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]], 2,
        function(x)
          ifelse(x >= 6, NA, x))

# Reverse code all 5-point scale items for self-esteem
# This ensures that higher scores consistently indicate higher self-esteem
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]] <- 
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]], 2, function(x)
    6 - x) 

# Calculate total self-esteem score by summing individual item scores
data$selfEst <- rowSums(data[, grep(pattern = "^H1PF\\d{2}", colnames(data))])


# # will those indicating an 8 (for S4 or ethnicity variable) mess things up? (same with sex)
# summary(data$sex) #summary(data$S2)
# summary(data$ethnicity) #summary(data$S4)
# # Change 8 to NA
# data$ethnicity <- ifelse(data$ethnicity == 8, NA, data$ethnicity)
# # sum(I(data$ethnicity == 8), na.rm = TRUE)

# Recode parental education
# 1 = High school or less, 2 = Some college, 3 = College graduate or more
data$momEdu <- ifelse(data$S12 %in% c(1, 2, 3, 10), 1, 
                      ifelse(data$S12 >= 4 & data$S12 <= 6, 2, 
                             ifelse(data$S12 %in% c(7, 8), 3, 
                                    NA)))

data$dadEdu <- ifelse(data$S18 %in% c(1, 2, 3, 10), 1, 
                      ifelse(data$S18 >= 4 & data$S18 <= 6, 2, 
                             ifelse(data$S18 %in% c(7, 8), 3, 
                                    NA)))

table(data$dadEdu, data$momEdu)
table(data$momEdu); sum(is.na(data$momEdu))
table(data$dadEdu); sum(is.na(data$dadEdu))

# Create a parental education variable (highest level of either parent)
data$parentalEdu <- apply(data[, stringr::str_detect(colnames(data), pattern = "Edu$")], 1, function(x) max(x, na.rm = FALSE))

# data %>% 
#   select("dadEdu", "momEdu", "parentalEdu") %>% 
#   head()

# Recode family structure
# S11: Living with mother figure, S17: Living with father figure
data$familyStruct <- data$S11 + data$S17
data$familyStruct <- ifelse(data$familyStruct == 10, NA, data$familyStruct) # Change "multiple response" to missing

table(data$S11); sum(is.na(data$S11))
table(data$S17); sum(is.na(data$S17))
table(data$S11, data$S17)

table(data$familyStruct)

# Recode health insurance gap (wave 1) [PC22: In the past 12 months, has there been a time when {NAME} had no health insurance?]
data <- data %>%
  rename(healthInsur = PC22) %>%
  mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know; 99=not applicable

# Calculate sports participation (H1NF12B)
# Sum number of sports and create a binary participation variable
data$sport <- rowSums(data[, colnames(data)[stringr::str_starts(colnames(data), pattern = "^S44A")]])
data$sportPartic <- ifelse(data$sport >= 1, 1, 
                           ifelse(data$sport == 0, 0, 
                                  NA))

table(data$sportPartic)


# Feelings Scale 
data %>% 
  select("H1FS1":"H1FS19")
# Change values for Feelings Scale items to NA
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]] <-
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]], 2, function(x)
    ifelse(x >= 6, NA, x)) # 6=refused; 8=don't know

# Reverse code items 4, 8, 11, & 15 on the Feelings Scale (4-point scale items; max value = 3)
data[, c("H1FS4", "H1FS8", "H1FS11", "H1FS15")] <-
  apply(data[, c("H1FS4", "H1FS8", "H1FS11", "H1FS15")], 2, function(x)
    3 - x)

# Score Feelings Scale 
data$feelings <- rowSums(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]])

# Import wave I sampling weights (to get CLUSTER2 to match to school weights later)
w1.w <- readr::read_tsv(file = "Application/Data/ICPSR_21600/DS0004/21600-0004-Data.tsv")

# Add cluster indicator (CLUSTER2)
data <- merge(data, w1.w[, c("AID", "CLUSTER2")], by = "AID")

# Reorder variables & add wave number to variable names
data <- data %>%
  select(AID, CLUSTER2, everything()) %>%
  select("AID", "CLUSTER2", "age", "sex",
         "ethnicity", "white", "black", "asian", "nativeAmerican", "raceOther",
         "healthInsur", "parentalEdu", "familyStruct", "sport", "sportPartic", "feelings", "selfEst")

colnames(data)[-c(1:2)] <- paste0(colnames(data)[-c(1:2)], "_w1")



## Wave III ----------------------------------------------------------------
# Wave III (21600-0008-Codebook) pg. 395

# Import student data
w3 <- readr::read_tsv(file = "Application/Data/ICPSR_21600/DS0008/21600-0008-Data.tsv")

# Drop variables in wave III student data
w3 <- w3 %>% 
  select("AID":"BIO_SEX3", 
         "H3SP7", "H3SP19":"H3SP22", # self-esteem scale [NOTE: STILL MISSING AN ITEM]
         "H3HS1") # Health insurance item 

# # Add wave III student data to data set
# data <- merge(data, w3, by = "AID")

# self-esteem scale
# Easterlin et al. (2019) article: 
# Uses responses to 5 statements to create a scale measure of self-esteem: “Do you agree or disagree with the following statement?” 
## “You have a lot of good qualities,” (H3SP19)
## “you have a lot to be proud of,” (H3SP20)
## “you like yourself just the way you are,” (H3SP21)
## and “you feel loved and wanted,”                     [CANNOT FIND COLUMN NAME FOR THIS QUESTION]
# used "Do you agree or disagree that you feel you are doing things just about right?" (H3SP22) as replacement 
### with responses on a 5-point scale ranging from 1 (strongly agree) to 5 (strongly disagree) (reverse coded), 
## and “How often was the following been true over the past week? You felt that you were just as good as other people,” (H3SP7)
### with responses ranging from 0 (never or rarely) to 3 (most of the time or all of the time). Higher scores indicate greater self-esteem.28

# change values for self-esteem items to NA 
w3$H3SP7 <- ifelse(w3$H3SP7 >= 6, NA, w3$H3SP7) # 6=refused; 8=don't know; 9=not applicable 

w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]] <-
  apply(w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]], 2,
        function(x)
          ifelse(x >= 96, NA, x)) # 96=refused; 98=don't know; 99=not applicable 

# reverse code the four 5-point scale items
w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]] <-
  apply(w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]], 2, function(x)
    6 - x) 

# Score self-esteem scale 
w3$selfEst <- rowSums(w3[, grep(pattern = "^H3SP", colnames(w3))])

# check 
# hist(w3$selfEst)

# # gap in health insurance (wave 3) 
# # 1. Over the past 12 months, how many months did you have health insurance? (H3HS1)
# w3 <- w3 %>%
#   rename(healthInsur = H3HS1) %>%
#   mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know; 99=not applicable

# Add wave III student data to data set
data <- merge(data, w3[, c("AID", "selfEst")], by = "AID")
# data <- merge(data, w3[, c("AID", "selfEst", "healthInsur")], by = "AID")

# Import Wave III weights 
# w3.w <- readr::read_tsv(file = "Empirical-Application/Data/ICPSR_21600/DS0018/21600-0018-Data.tsv") # GSWGT3_2 for level 1 weight 
# w3.sch.w <- readr::read_tsv(file = "Empirical-Application/Data/ICPSR_21600/DS0019/21600-0019-Data.tsv") # SCHWT1 W3_2_WC for level 2 weight 

# w3.w %>%
#   count(PTWGT3_2) %>%
#   dim() # 

# Add wave III weights to data set
# data <- merge(data, w3.sch.w, by = "CLUSTER2")
# data <- merge(data, w3.w, by = "AID") 

# Add wave number to variable names
colnames(data)[!grepl(pattern = paste(c("_w1$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))] <-
  paste0(colnames(data)[!grepl(pattern = paste(c("_w1$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))], "_w3")



## Wave IV -----------------------------------------------------------------

# Wave IV (21600-0022-Codebook) pg. 184 
w4 <- readr::read_tsv(file = "Application/Data/ICPSR_21600/DS0022/21600-0022-Data.tsv")

# Drop variables in wave I student data
w4 <- w4 %>% 
  select("AID":"BIO_SEX4", 
         "H4MH18":"H4MH27") # CES-D-10 scale [NOTE: DOUBLE CHECK THESE ARE CORRECT ITEMS]
         # "H4HS3") # Health insurance item 


# CES-D-10

# THESE ITEMS ARE FROM ONLINE 
# 1. I was bothered by things that usually don't bother me. (H4MH18)
# 2. I had trouble keeping my mind on what I was doing. (H4MH21)
# 3. I felt depressed. (H4MH22)
# 4. I felt that everything I did was an effort. (H4MH23)
# 5. I felt hopeful about the future. 
# 6. I felt fearful. 
# 7. My sleep was restless. 
# 8. I was happy. (H4MH24)
# 9. I felt lonely.
# 10. I could not "get going." 

# Add Health report (https://addhealth.cpc.unc.edu/wp-content/uploads/docs/data_briefs/Depressive_Symptoms_Data_Brief.pdf)
# These items ask respondents how many days in the last week they felt: “depressed,” “sad,” “like they could not shake off the blues,” “happy
# (reverse coded),” and “like life was not worth living.” The response scale for each item ranges from 0 (rarely, 0 days) to 3 (severe, 5-7 days).

#       [I BELIEVE THESE ARE ALL 10 OF CES-D-10 BUT NOT CERTAIN ABOUT THREE OF THE ITEMS]
# 18. You were bothered by things that usually don't bother you. (H4MH18)
# 19. (During the past seven days:) You could not shake off the blues, even with help from your family and your friends. (H4MH19)
# 20. (During the past seven days:) You felt you were just as good as other people. (H4MH20)
# 21. (During the past seven days:) You had trouble keeping your mind on what you were doing. (H4MH21)
# 22. (During the past seven days:) You felt depressed. (H4MH22)
# 23. (During the past seven days:) You felt that you were too tired to do things. (H4MH23)
# 24. (During the past seven days:) You felt happy. (H4MH24)
# 25. (During the past seven days:) You enjoyed life. (H4MH25)
# 26. (During the past seven days:) You felt sad. (H4MH26)
# 27. (During the past seven days:) You felt that people disliked you, during the past seven days. (H4MH27)

# Change values for CES-D-10 items to NA
w4[, colnames(w4)[str_detect(colnames(w4), pattern = "^H4MH")]] <-
  apply(w4[, colnames(w4)[str_detect(colnames(w4), pattern = "^H4MH")]], 2, function(x)
    ifelse(x >= 6, NA, x)) # 6=refused; 8=don't know

# Reverse code items 20, 24, & 25 on the CES-D-10 (4-point scale items; max value = 3)
w4[, c("H4MH20", "H4MH24", "H4MH25")] <-
  apply(w4[, c("H4MH20", "H4MH24", "H4MH25")], 2, function(x)
    3 - x)

# Score CES-D-10 (Depression) scale 
w4$depress <- rowSums(w4[, colnames(w4)[str_detect(colnames(w4), pattern = "^H4MH")]])

# # check 
# hist(w4$depress)
# table(I(w4$depress > 10)); sum(I(w4$depress > 10), na.rm = T) / sum(!is.na(w4$depress)) # 16% would be classified with depression in Easterlin study


# # gap in health insurance (wave 4) 
# # 3. Over the past 12 months, how many months did you have health insurance? (H4HS3)
# w4 <- w4 %>%
#   rename(healthInsur = H4HS3) %>%
#   mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know

# Add wave IV student data to data set
data <- merge(data, w4[, c("AID", "depress")], by = "AID")
# data <- merge(data, w4[, c("AID", "healthInsur", "depress")], by = "AID")

# GSWGT4_2 for level 1 weights
# SCHWT1 W4_2_WC for level 2 weights

# Import Wave IV weights 
# w4.w <- readr::read_tsv(file = "Empirical-Application/Data/ICPSR_21600/DS0031/21600-0031-Data.tsv") # 

# Add wave IV weights to data set
# data <- merge(data, w4.w[, -2], by = "AID")

# Add wave number to variable names
colnames(data)[!grepl(pattern = paste(c("_w\\d{1}$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))] <- 
  paste0(colnames(data)[!grepl(pattern = paste(c("_w\\d{1}$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))], "_w4")


# Calculated number of months not covered with health insurance 
# data$mnthsNotCovered <- 
#   rowSums(data[, colnames(data)[str_detect(colnames(data), pattern = "^healthInsur_")]])
# data$insuranceGap <- ifelse(data$mnthsNotCovered > 0, 1, 0)
# head(
#   data[, c(colnames(data)[str_detect(colnames(data), pattern = "^healthInsur_")], "mnthsNotCovered")]
# )

# Create insurance gap variable based only on wave I
data$healthInsur_w1 <- ifelse(data$healthInsur_w1 > 0, 1, 0)
# data$healthInsur_w3 <- ifelse(data$healthInsur_w3 > 0, 1, 0)
# data$healthInsur_w4 <- ifelse(data$healthInsur_w4 > 0, 1, 0)
# data$healthInsur_gap <- ifelse(data$healthInsur_w1 == 0 | data$healthInsur_w3 == 0 | data$healthInsur_w4 == 0,
#                                1, 0)

# Add cluster (school) size variable 
data <- data %>% 
  group_by(CLUSTER2) %>% 
  mutate(n = n()) %>% 
  ungroup()



# Drop small clusters -----------------------------------------------------

# First drop those missing sport participation (treatment)
data <- data %>% 
  filter(is.na(sportPartic_w1) != TRUE) # 4208 - 3155 = 1,053 dropped 

# Calculate cluster sizes  
data <- data %>% 
  group_by(CLUSTER2) %>% 
  mutate(n = n()) %>% 
  ungroup()

# Check clusters in order of size 
data %>% 
  group_by(CLUSTER2) %>% 
  summarize(n = max(n)) %>% 
  arrange(desc(n)) %>% 
  print(n = 132) # 123 schools (cluster sizes ranged from 1 to 81 students)

# Drop cluster sizes < 5
data <- data %>%
  filter(n >= 5) # To avoid convergence issues, we restricted our analysis to the schools with 5 or more students (121 schools)

# Look at cluster sizes
data %>% 
  group_by(CLUSTER2) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n = 132)  # 121 schools (cluster sizes ranged from 5 to 81 students)



# Missing Data ------------------------------------------------------------

# Check proportions of categorical variables & Scale variables (maybe create the scaled version of the mediator (selfEst_w3) after imputing data)
head(data)
colnames(data)

# head(scale(data))
summary(data[, c(
  "sex_w1",
  "ethnicity_w1",
  "white_w1",
  "black_w1",
  "asian_w1",
  "nativeAmerican_w1",
  "raceOther_w1"
)])   # Might need to collapse asian, native american, & other b/c the proportion of students in each group is only 0.0466, 0.058, & 0.0682 


summary(data[, c(
  "healthInsur_w1",
  "parentalEdu_w1",
  "familyStruct_w1", # Maybe drop this variable or look into its meaning more (has 3 distinct values: 0, 1, 2)
  "sport_w1",
  "sportPartic_w1",
  "feelings_w1",
  "depress_w4"
)])

# PS model includes: "sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 + parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2)"

# Need to scale: 
# Scale, rename, & drop non-scaled variable
data <- data %>%
  mutate(
    age_w1_sc = as.vector(scale(age_w1)),
    parentalEdu_w1_sc = as.vector(scale(parentalEdu_w1)),
    feelings_w1_sc = as.vector(scale(feelings_w1)), 
    selfEst_w1_sc = as.vector(scale(selfEst_w1))
  ) %>%
  select(-age_w1, -parentalEdu_w1, -feelings_w1)  # Remove original columns if desired

# Select variables in PS model 
t <- data %>% 
  select(!c(AID, CLUSTER2, sport_w1, n)) %>% 
  # select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_gap))
  # healthInsur_w3, healthInsur_w4, 
  # healthInsur_w1))
  select(c("sex_w1", "white_w1", "black_w1", 
           "sportPartic_w1", "selfEst_w3", "depress_w4", 
           "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_w1", 
           "age_w1_sc", "feelings_w1_sc", "selfEst_w1_sc"))

# Missing pattern 
md.pattern(t, rotate.names = TRUE)



# md.pairs(t)

# number of missing outcome variable & cluster size by cluster 
# data %>% 
#   group_by(CLUSTER2) %>% 
#   summarize(missing = sum(is.na(sportPartic_w1)), 
#             n = n()) %>% 
#   arrange(desc(missing)) %>% 
#   print(n = 132)

# Imputation 
imp <- mice(t, m = 1, seed = 875421)

# Create complete data 
t_imp <- complete(imp, "long", include = TRUE)

# head(t); head(t_imp[t_imp$.imp == 1, ])

# replace missing values in original data set
var <- colnames(t)
for (i in 1:length(var)) {
  
  data[, var[i]] <- t_imp[t_imp$.imp == 1, var[i]]
  
}

# Checking data after imputing values 
# mice::md.pattern(data[, var])
# summary(data$sportPartic_w1)
# head(data[, var])


# Add a scaled version of the mediator to use in outcome model 
data <- data %>%
  mutate(selfEst_w3_sc = as.vector(scale(selfEst_w3)))

# Re-check 

# Check proportions of categorical variables & Scale variables (maybe create the scaled version of the mediator (selfEst_w3) after imputing data)
head(data)
colnames(data)

head(scale(data))
summary(data[, c(
  "sex_w1",
  "ethnicity_w1",
  "white_w1",
  "black_w1",
  "asian_w1",
  "nativeAmerican_w1",
  "raceOther_w1"
)])   # Might need to collapse asian, native american, & other b/c the proportion of students in each group is only 0.0466, 0.058, & 0.0682 

# 
# summary(data[, c(
#   "healthInsur_w1",
#   "parentalEdu_w1",
#   "familyStruct_w1", # Maybe drop this variable or look into its meaning more (has 3 distinct values: 0, 1, 2)
#   "sport_w1",
#   "sportPartic_w1",
#   "feelings_w1",
#   "healthInsur_w3",
#   "healthInsur_w4",
#   "depress_w4"
# )])


# Clean up environment 
rm(imp, t, t_imp, 
   w1, w3, w4)



## ICC ---------------------------------------------------------------------

# ICC for mediator 
med_unconditional <- lme4::lmer(selfEst_w3 ~ (1 | CLUSTER2), data = data)
summary(med_unconditional)

med_var <- data.frame(lme4::VarCorr(med_unconditional))[1, 4] # variance due to level-2
med_res <- data.frame(lme4::VarCorr(med_unconditional))[2, 4] # residual variance 
med_icc <- med_var / (med_var + med_res)
med_icc # with 121 schools (5+ in size) med icc = 0.006954964

# ICC for outcome 
out_unconditional <- lme4::lmer(depress_w4 ~ (1 | CLUSTER2), data = data)
summary(out_unconditional)

out_var <- data.frame(lme4::VarCorr(out_unconditional))[1, 4] # variance due to level-2
out_res <- data.frame(lme4::VarCorr(out_unconditional))[2, 4] # residual variance 
out_icc <- out_var / (out_var + out_res)
out_icc # with 121 schools (5+ in size) outcome icc = 0.01868923

# ICC for some of the variables in PS models 
# ## Parent education
# parent_uncond <- lme4::lmer(parentalEdu_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(parent_uncond))[1, 4] / ( data.frame(lme4::VarCorr(parent_uncond))[1, 4] + data.frame(lme4::VarCorr(parent_uncond))[2, 4])
# 
# ## Feeling scale 
# feeling_uncond <- lme4::lmer(feelings_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(feeling_uncond))[1, 4] / ( data.frame(lme4::VarCorr(feeling_uncond))[1, 4] + data.frame(lme4::VarCorr(feeling_uncond))[2, 4])
# 
# ## age
# age_uncond <- lme4::lmer(age_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(age_uncond))[1, 4] / ( data.frame(lme4::VarCorr(age_uncond))[1, 4] + data.frame(lme4::VarCorr(age_uncond))[2, 4])
# 
# ## self-esteem
# selfEst_uncond <- lme4::lmer(selfEst_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(selfEst_uncond))[1, 4] / ( data.frame(lme4::VarCorr(selfEst_uncond))[1, 4] + data.frame(lme4::VarCorr(selfEst_uncond))[2, 4])
# 
# ## white 
# white_uncod <- lme4::lmer(white_w1 ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(white_uncod))[1, 4] / ( data.frame(lme4::VarCorr(white_uncod))[1, 4] + data.frame(lme4::VarCorr(white_uncod))[2, 4])



# PS & IPTW Calculation ---------------------------------------------------------------
# This script calculates propensity scores and Inverse Probability of Treatment Weights (IPTW)
# for three models: SL (Standard Logistic), FE (Fixed Effects), and RE (Random Effects).

### SL Propensity Score Model -----------------------------------------------
# Standard logistic regression model for propensity score estimation
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc",
                family = "binomial", 
                data = data)

# Predict propensity scores and log-odds for SL model
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")

# Calculate IPTW for SL model
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))

### FE Propensity Score Model -----------------------------------------------
# Fixed effects logistic regression model for propensity score estimation
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)

# Predict propensity scores and log-odds for FE model
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")

# Calculate IPTW for FE model
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))

### RE Propensity Score Model -----------------------------------------------
# Random effects logistic regression model for propensity score estimation
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc + 
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) 

# Predict propensity scores and log-odds for RE model
data$ps_re <- predict(psmod_re, type = "response")
data$ps_re_logit <- predict(psmod_re, type = "link")

# Calculate IPTW for RE model
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))



# Truncate Non-overlap Cases ---------------------------------------------------------

### SL (Single Level) ------------------------------------------------------

# Count instances where propensity scores (PS) are outside the valid range
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_sl < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_sl > 0.999)))

# Count cases below the 1st percentile and above the 99th percentile of IPTW
first_percentile <- quantile(data$iptw_sl, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_sl, probs = 0.99)

paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_sl < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_sl > ninety_ninth_percentile))
)

# Adjust IPTW values to the 1st and 99th percentile thresholds
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < first_percentile, first_percentile, 
                          ifelse(iptw_sl > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_sl)))

### FE (Fixed Effects) -------------------------------------------------------------

# Count instances where propensity scores (PS) are outside the valid range
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_fe < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_fe > 0.999)))

# Count cases below the 1st percentile and above the 99th percentile of IPTW
first_percentile <- quantile(data$iptw_fe, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_fe, probs = 0.99)

paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_fe < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_fe > ninety_ninth_percentile))
)

# Adjust IPTW values to the 1st and 99th percentile thresholds
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < first_percentile, first_percentile, 
                          ifelse(iptw_fe > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_fe)))

### RE (Random Effects) -------------------------------------------------------------

# Count instances where propensity scores (PS) are outside the valid range
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_re < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_re > 0.999)))

# Count cases below the 1st percentile and above the 99th percentile of IPTW
first_percentile <- quantile(data$iptw_re, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_re, probs = 0.99)

paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_re < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_re > ninety_ninth_percentile))
)

# Adjust IPTW values to the 1st and 99th percentile thresholds
data <- data %>% 
  mutate(iptw_re = ifelse(iptw_re < first_percentile, first_percentile, 
                          ifelse(iptw_re > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_re)))



# Covariate Balance Visuals -----------------------------------------------

# Custom function to calculate weighted variance
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    # Remove NA values
    valid <- !is.na(x) & !is.na(w)
    x <- x[valid]
    w <- w[valid]
  }
  
  # Calculate the weighted mean
  wm <- weighted.mean(x, w)
  
  # Calculate the weighted variance
  variance <- sum(w * (x - wm)^2) / sum(w)
  
  return(variance)
}

# Function to calculate Standardized Mean Difference (SMD)
calculate_smd <- function(data, treatment, covariate) {
  # Calculate the mean for the treatment group
  mean_treatment <- mean(data[[covariate]][data[[treatment]] == 1], na.rm = TRUE)
  
  # Calculate the mean for the control group
  mean_control <- mean(data[[covariate]][data[[treatment]] == 0], na.rm = TRUE)
  
  # Calculate the standard deviation for the treatment group
  sd_treatment <- sd(data[[covariate]][data[[treatment]] == 1], na.rm = TRUE)
  
  # Calculate the standard deviation for the control group
  sd_control <- sd(data[[covariate]][data[[treatment]] == 0], na.rm = TRUE)
  
  # Compute the SMD using the means and pooled standard deviations
  smd <- (mean_treatment - mean_control) / sqrt((sd_treatment^2 + sd_control^2) / 2)
  
  return(smd)  # Return the calculated SMD
}

# List of covariates to analyze for balance
covariates <- c("feelings_w1_sc", "sex_w1", "age_w1_sc", "white_w1", 
                "black_w1", "parentalEdu_w1_sc", "familyStruct_w1", 
                "healthInsur_w1", "selfEst_w1_sc")

# Calculate SMD for each covariate before applying weights
smd_before <- sapply(covariates, function(cov) calculate_smd(data, "sportPartic_w1", cov))

# Convert the SMD results into a data frame for easier manipulation
smd_before_df <- data.frame(covariate = covariates, SMD = smd_before)
smd_before_df$type <- "Unweighted"  # Label the type as 'Unweighted'

# Function to calculate weighted SMD, allowing for weights as an argument
calculate_weighted_smd <- function(data, treatment, covariate, weights_col) {
  # Subset the data into treatment and control groups
  treatment_group <- data[data[[treatment]] == 1, ]
  control_group <- data[data[[treatment]] == 0, ]
  
  # Check if either group is empty; return NA if so
  if (nrow(treatment_group) == 0 || nrow(control_group) == 0) {
    return(NA)  # Return NA if either group is empty
  }
  
  # Remove rows with missing values in the covariate or weights
  treatment_group <- treatment_group[!is.na(treatment_group[[covariate]]) & !is.na(treatment_group[[weights_col]]), ]
  control_group <- control_group[!is.na(control_group[[covariate]]) & !is.na(control_group[[weights_col]]), ]
  
  # Check lengths after removing NAs; return NA if either group is empty
  if (nrow(treatment_group) == 0 || nrow(control_group) == 0) {
    return(NA)  # Return NA if either group is empty after filtering
  }
  
  # Calculate weighted means for both groups
  mean_treatment <- weighted.mean(treatment_group[[covariate]], treatment_group[[weights_col]], na.rm = TRUE)
  mean_control <- weighted.mean(control_group[[covariate]], control_group[[weights_col]], na.rm = TRUE)
  
  # Calculate weighted standard deviations for both groups
  sd_treatment <- sqrt(weighted.var(treatment_group[[covariate]], treatment_group[[weights_col]], na.rm = TRUE))
  sd_control <- sqrt(weighted.var(control_group[[covariate]], control_group[[weights_col]], na.rm = TRUE))
  
  # Compute the SMD using the weighted means and pooled standard deviations
  smd <- (mean_treatment - mean_control) / sqrt((sd_treatment^2 + sd_control^2) / 2)
  
  return(smd)  # Return the calculated weighted SMD
}

# Calculate SMD for each covariate after applying different weighting methods
smd_sl_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_sl"))
smd_fe_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_fe"))
smd_re_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_re"))

# Convert the results into data frames for easier handling
smd_sl_after_df <- data.frame(covariate = covariates, SMD = smd_sl_after, type = "Single-Level")
smd_fe_after_df <- data.frame(covariate = covariates, SMD = smd_fe_after, type = "Fixed-Effect")
smd_re_after_df <- data.frame(covariate = covariates, SMD = smd_re_after, type = "Random-Effects")

# Combine all SMD data frames into one for comprehensive analysis
smd_combined <- rbind(smd_before_df, smd_sl_after_df, smd_fe_after_df, smd_re_after_df)

# Calculate the Absolute Standardized Mean Difference (ASMD)
smd_combined$ASMD <- abs(smd_combined$SMD)

# Define a custom order for the y-axis in the plot
custom_order <- c("black_w1", "white_w1", "familyStruct_w1", 
                  "healthInsur_w1", "age_w1_sc", "sex_w1", 
                  "selfEst_w1_sc", "feelings_w1_sc", "parentalEdu_w1_sc")

# Define new labels for the y-axis to enhance readability
new_labels <- c("Race: Black", "Race: White", "Family Structure", 
                "Health Insurance \n Coverage Gap", "Age", "Sex", 
                "Self-Esteem Score", "Feelings Scale Score", "Parental Education")


# Create a Love Plot to visualize the Absolute SMD, using custom ordering
## Save visual 
# svg(filename = "Application/Output/Covariate-Balance.svg")
pdf("Application/Output/Covariate-Balance_Love-Plot.pdf")
## Visual 
ggplot(smd_combined, aes(x = ASMD, y = factor(covariate, levels = custom_order), color = type, shape = type)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +  # Reference line for SMD threshold
  geom_vline(xintercept = 0, color = "black") +  # Line at zero for clarity
  geom_point(size = 3, stroke = 1.5) +  # Points with increased size and stroke for emphasis
  labs(title = "Love Plot",
       subtitle = "Covariate Balance of Individual-Level Covariates",
       x = "Absolute Standardized Mean Difference (ASMD)",
       y = "") +
  theme_minimal() +  # Apply a minimal theme for clean aesthetics
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10),  # Adjust y-axis text for clarity
        axis.title = element_text(size = 14),  # Increase axis title size for visibility
        plot.title = element_text(size = 16, face = "bold"),  # Bold title for emphasis
        plot.subtitle = element_text(size = 14),  # Subtitle size for readability
        legend.position = "top") +  # Position legend at the top
  scale_color_manual(values = c("Unweighted" = "#1f77b4",  # Blue for unweighted
                                "Single-Level" = "#2ca02c",  # Green for single-level weighting
                                "Fixed-Effect" = "#ff7f0e",  # Orange for fixed-effect weighting
                                "Random-Effects" = "#9467bd"),  # Purple for random-effects weighting
                     name = NULL) +
  scale_shape_manual(values = c("Unweighted" = 16,  # Circle shape for unweighted
                                "Single-Level" = 17,  # Triangle shape for single-level
                                "Fixed-Effect" = 15,  # Square shape for fixed-effect
                                "Random-Effects" = 18),  # Diamond shape for random-effects
                     name = NULL) +
  scale_y_discrete(labels = new_labels)  # Use new labels for the y-axis

dev.off()
# Save visual 
ggsave(filename = "Application/Output/Covariate-Balance_Love-Plot.png", plot = last_plot())



# Create a Love Plot for paper 
library(extrafont)
loadfonts()
## Save visual 
pdf("Application/Output/Covariate-Balance_QP-Doc.pdf")
## Visual 
ggplot(smd_combined, aes(x = ASMD, y = factor(covariate, levels = custom_order), color = type, shape = type)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +  # Reference line for SMD threshold
  geom_vline(xintercept = 0, color = "black") +  # Line at zero for clarity
  geom_point(size = 3, stroke = 1.5) +  # Points with increased size and stroke for emphasis
  labs(#title = "Love Plot",
       #subtitle = "Covariate Balance of Individual-Level Covariates",
       x = "\n Absolute Standardized Mean Difference (ASMD)",
       y = "") +
  theme_minimal() +  # Apply a minimal theme for clean aesthetics
  theme(text = element_text(family = "Times New Roman"),
        axis.text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 10),  # Adjust y-axis text for clarity
        axis.title = element_text(size = 14),  # Increase axis title size for visibility
        plot.title = element_text(size = 16, face = "bold"),  # Bold title for emphasis
        plot.subtitle = element_text(size = 14),  # Subtitle size for readability
        legend.position = "top") +  # Position legend at the top
  scale_color_manual(values = c("Unweighted" = "#1f77b4",  # Blue for unweighted
                                "Single-Level" = "#2ca02c",  # Green for single-level weighting
                                "Fixed-Effect" = "#ff7f0e",  # Orange for fixed-effect weighting
                                "Random-Effects" = "#9467bd"),  # Purple for random-effects weighting
                     name = NULL) +
  scale_shape_manual(values = c("Unweighted" = 16,  # Circle shape for unweighted
                                "Single-Level" = 17,  # Triangle shape for single-level
                                "Fixed-Effect" = 15,  # Square shape for fixed-effect
                                "Random-Effects" = 18),  # Diamond shape for random-effects
                     name = NULL) +
  # geom_text()
  scale_y_discrete(labels = new_labels)  # Use new labels for the y-axis

dev.off()
# Save visual 
ggsave(filename = "Application/Output/Covariate-Balance_QP-Doc.png", plot = last_plot())






# 
# # ADD VISUALS TO SHOW HISTOGRAM OR DENSITY OVERLAP 
# # ADD VISUAL SHOWING DOT PLOT ACROSS MODELS ON EACH COVARIATE OF SAMPLE BALANCE (SMD LIKE LOVEPLOT)
# 
# colnames(data)
# ## THIS IS BEFORE TRUNCATION 
# # SL 
# data %>% 
#   mutate(sportPartic_w1 = as.factor(sportPartic_w1)) %>% 
#   ggplot(aes(x = ps_sl, 
#              group = sportPartic_w1, 
#              fill = sportPartic_w1)) +
#   geom_histogram(position = "identity",
#                  alpha = 0.5,
#                  binwidth = 0.01) +
#   geom_density(aes(y = ..count.. * 0.01, 
#                    color = sportPartic_w1),  
#                size = 0.5,        
#                alpha = 0, 
#                trim = TRUE, 
#                show.legend = FALSE) +      # Set alpha to 1 (no transparency)
#   scale_fill_manual(values = c("blue", "darkorange")) +
#   scale_color_manual(values = c("blue", "darkorange")) +
#   theme_minimal() +
#   labs(# title = "Histogram of ps_sl_logit by sportPartic_w1",
#     y = "count",
#     fill = "trt (Sport Participation)") +
#   theme(legend.position = "bottom")
# 
# # FE 
# data %>% 
#   mutate(sportPartic_w1 = as.factor(sportPartic_w1)) %>% 
#   ggplot(aes(x = ps_fe, 
#              group = sportPartic_w1, 
#              fill = sportPartic_w1)) +
#   geom_histogram(position = "identity",
#                  alpha = 0.5,
#                  binwidth = 0.01) +
#   geom_density(aes(y = ..count.. * 0.01, 
#                    color = sportPartic_w1),  
#                size = 0.5,        
#                alpha = 0, 
#                trim = TRUE, 
#                show.legend = FALSE) +      # Set alpha to 1 (no transparency)
#   scale_fill_manual(values = c("blue", "darkorange")) +
#   scale_color_manual(values = c("blue", "darkorange")) +
#   theme_minimal() +
#   labs(# title = "Histogram of ps_sl_logit by sportPartic_w1",
#     y = "count",
#     fill = "trt (Sport Participation)") +
#   theme(legend.position = "bottom")
# 
# # RE 
# data %>% 
#   mutate(sportPartic_w1 = as.factor(sportPartic_w1)) %>% 
#   ggplot(aes(x = ps_re, 
#              group = sportPartic_w1, 
#              fill = sportPartic_w1)) +
#   geom_histogram(position = "identity",
#                  alpha = 0.5,
#                  binwidth = 0.01) +
#   geom_density(aes(y = ..count.. * 0.01, 
#                    color = sportPartic_w1),  
#                size = 0.5,        
#                alpha = 0, 
#                trim = TRUE, 
#                show.legend = FALSE) +      # Set alpha to 1 (no transparency)
#   scale_fill_manual(values = c("blue", "darkorange")) +
#   scale_color_manual(values = c("blue", "darkorange")) +
#   theme_minimal() +
#   labs(# title = "Histogram of ps_sl_logit by sportPartic_w1",
#     y = "count",
#     fill = "trt (Sport Participation)") +
#   theme(legend.position = "bottom")
#   





# Estimate Effects --------------------------------------------------------

## Mediator models ---------------------------------------------------------

#### Single-Level (SL) Models ------------------------------------------------

# Model: SL PS & SL - Mediation/Outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_sl
  )

# Model: FE PS & SL - Mediation/Outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_fe
  )

# Model: RE PS & SL - Mediation/Outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_re
  )


#### Fixed-Effect (FE) Models ------------------------------------------------

# Model: SL PS & FE - Mediation/Outcome
med_slfe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_sl
)

# Model: FE PS & FE - Mediation/Outcome 
med_fefe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_fe
)

# Model: RE PS & FE - Mediation/Outcome
med_refe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_re
)


#### Random-Effect (RE) Models ------------------------------------------------

### Add a column of ones for level-2 weights
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Model: SL PS & RE - Mediation/Outcome
med_slre <-
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )

# Model: FE PS & RE - Mediation/Outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )

# Model: RE PS & RE - Mediation/Outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )



## Outcome Models (TNDE & PNIE) --------------------------------------------

#### Single-Level (SL) Models ------------------------------------------------

# Model: SL PS & SL Mediation/Outcome
out_slsl <-
  glm(
    formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_sl
  )

# Model: FE PS & SL Mediation/Outcome
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_fe
  )

# Model: RE PS & SL Mediation/Outcome
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1",
    data = data,
    weights = iptw_re
  )

#### Fixed-Effect (FE) Models ------------------------------------------------

# Model: SL PS & FE Mediation/Outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)

# Model: FE PS & FE Mediation/Outcome
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Model: RE PS & FE Mediation/Outcome
out_refe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
      data = data,
      weights = iptw_re
)

#### Random-Effect (RE) Models ------------------------------------------------

# Prepare data for Random-Effect Models: Add a column for level-2 weights
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Model: SL PS & RE Mediation/Outcome
out_slre <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )

# Model: FE PS & RE Mediation/Outcome
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )

# Model: RE PS & RE Mediation/Outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )



## Outcome Models (PNDE & TNIE) --------------------------------------------

#### Single-Level (SL) Models ------------------------------------------------

# Model: SL PS & SL Mediation/Outcome
out_slsl_interac <-
  glm(
    formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
    selfEst_w3_sc:sportPartic_w1",
    data = data,
    weights = iptw_sl
  )

# Model: FE PS & SL Mediation/Outcome
out_fesl_interac <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
    selfEst_w3_sc:sportPartic_w1",
    data = data,
    weights = iptw_fe
  )

# Model: RE PS & SL Mediation/Outcome
out_resl_interac <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
    selfEst_w3_sc:sportPartic_w1",
    data = data,
    weights = iptw_re
  )

#### Fixed-Effect (FE) Models ------------------------------------------------

# Model: SL PS & FE Mediation/Outcome
out_slfe_interac <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)

# Model: FE PS & FE Mediation/Outcome
out_fefe_interac <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Model: RE PS & FE Mediation/Outcome
out_refe_interac <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + as.factor(CLUSTER2)", 
      data = data,
      weights = iptw_re
)

#### Random-Effect (RE) Models ------------------------------------------------

# Prepare data for Random-Effect Models: Add a column for level-2 weights
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Model: SL PS & RE Mediation/Outcome
out_slre_interac <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )

# Model: FE PS & RE Mediation/Outcome
out_fere_interac <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )

# Model: RE PS & RE Mediation/Outcome
out_rere_interac <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + 
      selfEst_w3_sc:sportPartic_w1 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )



# Display Estimates -------------------------------------------------------

# Define conditions
conditions <- c("slsl", "fesl", "resl", "slfe", "fefe", "refe", "slre", "fere", "rere")

# Extract TNDE estimates
TNDE <- sapply(conditions, function(cond) {
  model_name <- paste0("out_", cond)
  summary(get(model_name))$coef["sportPartic_w1", "Estimate"]
})

# Extract PNIE estimates
PNIE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond)
  summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"] * 
    summary(get(out_model_name))$coef["selfEst_w3_sc", "Estimate"]
})

# Extract PNDE estimates
PNDE <- sapply(conditions, function(cond) {
  model_name <- paste0("out_", cond, "_interac")
  summary(get(model_name))$coef["sportPartic_w1", "Estimate"]
})

# Extract TNIE estimates
TNIE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond, "_interac")
  summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"] * 
    summary(get(out_model_name))$coef["selfEst_w3_sc", "Estimate"] + 
    summary(get(out_model_name))$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"]
})




# Create results DataFrame
results_DF <- data.frame(
  cond = conditions,
  TNDE = TNDE,
  PNDE = PNDE,
  PNIE = PNIE,
  TNIE = TNIE
)

# Display results
rownames(results_DF) <- NULL
results_DF
#   cond       TNDE       PNDE        PNIE          TNIE
# 1 slsl -0.2764533 -0.2764801 -0.03027051  0.0180582390
# 2 fesl -0.1973009 -0.1976276 -0.09752048 -0.0367624900
# 3 resl -0.2400019 -0.2401109 -0.07276088 -0.0045414414
# 4 slfe -0.1804319 -0.1808842 -0.09704172 -0.0363361536
# 5 fefe -0.1889825 -0.1894340 -0.09570423 -0.0002112213
# 6 refe -0.1958658 -0.1963458 -0.09240229  0.0030725382
# 7 slre -0.2208416 -0.2210517 -0.06717527 -0.0117677236
# 8 fere -0.1946794 -0.1951001 -0.09536315 -0.0126787393
# 9 rere -0.2171489 -0.2174313 -0.08069368  0.0039460448



# Conduct Bootstrap Confidence Intervals (CI) ---------------------------------

## TNDE & PNIE -------------------------------------------------------------

#### Single-Level (SL) Med/Out Models  ---------------------------------------
# SL PS Model 
slsl_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                iptw = iptw_sl,
                                data = data,
                                model = "SL",
                                cores = 6,
                                core_seeds = c(4561:4566), 
                                effect_type = "PNIE")

# FE PS Model 
fesl_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                              iptw = iptw_fe,
                              data = data,
                              model = "SL",
                              cores = 6,
                              core_seeds = c(4561:4566), 
                              effect_type = "PNIE")

# RE PS Model 
resl_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                              iptw = iptw_re,
                              data = data,
                              model = "SL",
                              cores = 6,
                              core_seeds = c(4561:4566),
                              effect_type = "PNIE")

#### Fixed-Effect (FE) Med/Out Models ----------------------------------------
# SL PS Model 
slfe_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                              iptw = iptw_sl,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566), 
                              effect_type = "PNIE")

# FE PS Model 
fefe_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                              iptw = iptw_fe,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566), 
                              effect_type = "PNIE")

# RE PS Model 
refe_ci_PNIE <- bootstrap_ci_paral_2(iterations = 1000,
                              iptw = iptw_re,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566), 
                              effect_type = "PNIE")

#### Random-Effect (RE) Med/Out Models ---------------------------------------
# SL PS Model 
execution_time_slre <- system.time({ # Track computation time 
  slre_ci_PNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                   iptw = iptw_sl, 
                                   data = data, 
                                   cores = 6, 
                                   core_seeds = c(4561:4566), 
                                   effect_type = "PNIE")
})
# Print the execution time
print(execution_time_slre)
#    user   system  elapsed 
# 8656.601 1124.463 3393.450 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_slre["elapsed"], "seconds\n")
# Elapsed time: 3393.45 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", slre_ci_PNIE$mediator_converged_count, 
       " (", (slre_ci_PNIE$mediator_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", slre_ci_PNIE$outcome_converged_count, 
       " (", (slre_ci_PNIE$outcome_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", slre_ci_PNIE$both_converged_count, 
       " (", (slre_ci_PNIE$both_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1146 (76.4%)"

# FE PS Model
execution_time_fere <- system.time({ # Track computation time 
  fere_ci_PNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                   iptw = iptw_fe, 
                                   data = data, 
                                   cores = 6, 
                                   core_seeds = c(4561:4566), 
                                   effect_type = "PNIE")
})
# Print the execution time
print(execution_time_fere)
#     user    system   elapsed 
# 12087.304  1923.361  6444.097 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_fere["elapsed"], "seconds\n")
# Elapsed time: 6444.097 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", fere_ci_PNIE$mediator_converged_count, 
       " (", (fere_ci_PNIE$mediator_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", fere_ci_PNIE$outcome_converged_count, 
       " (", (fere_ci_PNIE$outcome_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", fere_ci_PNIE$both_converged_count, 
       " (", (fere_ci_PNIE$both_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1146 (76.4%)"

# RE PS Model
execution_time_rere <- system.time({ # Track computation time 
  rere_ci_PNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                   iptw = iptw_re, 
                                   data = data, 
                                   cores = 6, 
                                   core_seeds = c(4561:4566), 
                                   effect_type = "PNIE")
})
# Print the execution time
print(execution_time_rere)
#     user    system   elapsed 
# 10658.977  1202.641  2922.896 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_rere["elapsed"], "seconds\n")
# Elapsed time: 2922.896 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", rere_ci_PNIE$mediator_converged_count, 
       " (", (rere_ci_PNIE$mediator_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", rere_ci_PNIE$outcome_converged_count, 
       " (", (rere_ci_PNIE$outcome_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", rere_ci_PNIE$both_converged_count, 
       " (", (rere_ci_PNIE$both_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1146 (76.4%)"


## PNDE & TNIE -------------------------------------------------------------

#### Single-Level (SL) Med/Out Models  ---------------------------------------
# SL PS Model 
slsl_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_sl,
                                     data = data,
                                     model = "SL",
                                     cores = 6,
                                     core_seeds = c(4561:4566), 
                                     effect_type = "TNIE")

# FE PS Model 
fesl_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_fe,
                                     data = data,
                                     model = "SL",
                                     cores = 6,
                                     core_seeds = c(4561:4566), 
                                     effect_type = "TNIE")

# RE PS Model 
resl_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_re,
                                     data = data,
                                     model = "SL",
                                     cores = 6,
                                     core_seeds = c(4561:4566),
                                     effect_type = "TNIE")

#### Fixed-Effect (FE) Med/Out Models ----------------------------------------
# SL PS Model 
slfe_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_sl,
                                     data = data,
                                     model = "FE",
                                     cores = 6,
                                     core_seeds = c(4561:4566), 
                                     effect_type = "TNIE")

# FE PS Model 
fefe_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_fe,
                                     data = data,
                                     model = "FE",
                                     cores = 6,
                                     core_seeds = c(4561:4566), 
                                     effect_type = "TNIE")

# RE PS Model 
refe_ci_TNIE <- bootstrap_ci_paral_2(iterations = 1000,
                                     iptw = iptw_re,
                                     data = data,
                                     model = "FE",
                                     cores = 6,
                                     core_seeds = c(4561:4566), 
                                     effect_type = "TNIE")

#### Random-Effect (RE) Med/Out Models ---------------------------------------
# SL PS Model 
execution_time_slre <- system.time({ # Track computation time 
  slre_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                          iptw = iptw_sl, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_slre)
#    user   system  elapsed 
# 6963.170 1059.563 2095.475 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_slre["elapsed"], "seconds\n")
# Elapsed time: 2095.475 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", slre_ci_TNIE$mediator_converged_count, 
       " (", (slre_ci_TNIE$mediator_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", slre_ci_TNIE$outcome_converged_count, 
       " (", (slre_ci_TNIE$outcome_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", slre_ci_TNIE$both_converged_count, 
       " (", (slre_ci_TNIE$both_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1144 (76.2666666666667%)"

# FE PS Model
execution_time_fere <- system.time({ # Track computation time 
  fere_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                          iptw = iptw_fe, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_fere)
#   user   system  elapsed 
# 6804.245  544.192 3079.205 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_fere["elapsed"], "seconds\n")
# Elapsed time: 3079.205 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", fere_ci_TNIE$mediator_converged_count, 
       " (", (fere_ci_TNIE$mediator_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", fere_ci_TNIE$outcome_converged_count, 
       " (", (fere_ci_TNIE$outcome_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", fere_ci_TNIE$both_converged_count, 
       " (", (fere_ci_TNIE$both_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1144 (76.2666666666667%)"

# RE PS Model
execution_time_rere <- system.time({ # Track computation time 
  rere_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 1500, 
                                          iptw = iptw_re, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_rere)
#    user   system  elapsed 
# 6713.972  607.371 1913.623 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_rere["elapsed"], "seconds\n")
# Elapsed time: 1913.623 seconds

# Print convergence statistics
paste0("Number of converged mediator models: ", rere_ci_TNIE$mediator_converged_count, 
       " (", (rere_ci_TNIE$mediator_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", rere_ci_TNIE$outcome_converged_count, 
       " (", (rere_ci_TNIE$outcome_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", rere_ci_TNIE$both_converged_count, 
       " (", (rere_ci_TNIE$both_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)")
# [1] "Number of converged mediator models: 1154 (76.9333333333333%)"
# [1] "Number of converged outcome models: 1488 (99.2%)"
# [1] "Number of iterations with both models converged: 1144 (76.2666666666667%)"



# Store & Join Results ----------------------------------------------------

#### RE Mediator/Outcome Models CI -----------------------------------------

###### Obtain 1,000 completed iterations -----------------------------------

# Function to get non-NA pairs (i.e., first 1,000 completed iterations)
get_non_na_pairs <- function(direct, indirect, n = 1000) {
  combined <- data.frame(direct = direct, indirect = indirect)  # Combine vectors into a dataframe
  combined <- na.omit(combined)  # Remove rows with any NA values
  return(head(combined, n))  # Return the first n rows (or all if less than n)
}

# Apply the function to our data
slre_ci_PNIE_DF <- get_non_na_pairs(slre_ci_PNIE$direct_effects, slre_ci_PNIE$indirect_effects) 
fere_ci_PNIE_DF <- get_non_na_pairs(fere_ci_PNIE$direct_effects, fere_ci_PNIE$indirect_effects)
rere_ci_PNIE_DF <- get_non_na_pairs(rere_ci_PNIE$direct_effects, rere_ci_PNIE$indirect_effects)

slre_ci_TNIE_DF <- get_non_na_pairs(slre_ci_TNIE$direct_effects, slre_ci_TNIE$indirect_effects) 
fere_ci_TNIE_DF <- get_non_na_pairs(fere_ci_TNIE$direct_effects, fere_ci_TNIE$indirect_effects)
rere_ci_TNIE_DF <- get_non_na_pairs(rere_ci_TNIE$direct_effects, rere_ci_TNIE$indirect_effects)


###### Store RE Med/Outcome Model CIs --------------------------------------

# Create the results dataframe
results_DF_RE <- data.frame(
  cond = c("slre", "fere", "rere"),
  PNIE_LL = numeric(3),
  PNIE_UL = numeric(3),
  TNDE_LL = numeric(3),
  TNDE_UL = numeric(3),
  TNIE_LL = numeric(3),
  TNIE_UL = numeric(3),
  PNDE_LL = numeric(3),
  PNDE_UL = numeric(3),
  stringsAsFactors = FALSE
)

# List of dataframes to process
df_list_PNIE <- list(slre_ci_PNIE_DF, fere_ci_PNIE_DF, rere_ci_PNIE_DF)
df_list_TNIE <- list(slre_ci_TNIE_DF, fere_ci_TNIE_DF, rere_ci_TNIE_DF)

# Calculate CIs and fill the dataframe
for (i in 1:3) {
  results_DF_RE[i, c("PNIE_LL", "PNIE_UL")] <- quantile(df_list_PNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("TNDE_LL", "TNDE_UL")] <- quantile(df_list_PNIE[[i]]$direct, probs = c(0.025, 0.975))
  
  results_DF_RE[i, c("TNIE_LL", "TNIE_UL")] <- quantile(df_list_TNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("PNDE_LL", "PNDE_UL")] <- quantile(df_list_TNIE[[i]]$indirect, probs = c(0.025, 0.975))
}

# results_DF_RE


###### Join CI & point estimates for RE med/outcome ------------------------

# Merge results with existing data
results_DF_RE <- merge(results_DF[results_DF$cond %in% c("slre", "fere", "rere"), ], 
                       results_DF_RE)

# Clean up environment (drop functions & objects from this section that are no longer needed)
rm(get_non_na_pairs, slre_ci_PNIE_DF, slre_ci_TNIE_DF, fere_ci_PNIE_DF, fere_ci_TNIE_DF, rere_ci_PNIE_DF, rere_ci_TNIE_DF)


#### SL & FE Mediator/Outcome Models CI ------------------------------------

###### Join CI & point estimates for SL & FE med/outcome -------------------

# List of your lists
list_names <- c("slsl_ci_TNIE", "slsl_ci_PNIE", 
                "fesl_ci_TNIE", "fesl_ci_PNIE",
                "resl_ci_TNIE", "resl_ci_PNIE", 
                "slfe_ci_TNIE", "slfe_ci_PNIE", 
                "fefe_ci_TNIE", "fefe_ci_PNIE",
                "refe_ci_TNIE", "refe_ci_PNIE")
lists <- list(slsl_ci_TNIE, slsl_ci_PNIE, 
              fesl_ci_TNIE, fesl_ci_PNIE, 
              resl_ci_TNIE, resl_ci_PNIE, 
              slfe_ci_TNIE, slfe_ci_PNIE, 
              fefe_ci_TNIE, fefe_ci_PNIE, 
              refe_ci_TNIE, refe_ci_PNIE)

# Extracting 'indirect_ci'
indirect_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$indirect_ci
  data.frame(
    list_name = list_names[i],
    indirect_ci_LL = ci_values[1],  # First value
    indirect_ci_UL = ci_values[2],  # Second value
    stringsAsFactors = FALSE
  )
}))

# Extracting 'direct_ci'
direct_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$direct_ci
  data.frame(
    list_name = list_names[i],
    direct_ci_LL = ci_values[1],  # First value
    direct_ci_UL = ci_values[2],  # Second value
    stringsAsFactors = FALSE
  )
}))

# Combine both into a single data frame
combined_df <- merge(indirect_df, direct_df, by = "list_name")

# Create the 'effect' column & modify 'list_name' to keep only the first 4 characters
combined_df <- combined_df %>%
  mutate(effect = substr(list_name, nchar(list_name) - 3, nchar(list_name))) %>% 
  mutate(list_name = substr(list_name, 1, 4))

# Pivot the dataframe to wide format
combined_df_wide <- combined_df %>%
  pivot_longer(cols = c(indirect_ci_LL, indirect_ci_UL, direct_ci_LL, direct_ci_UL), 
               names_to = "variable", values_to = "value") %>%
  unite("variable", effect, variable, sep = "_") %>%
  pivot_wider(names_from = variable, values_from = value) %>% 
  as.data.frame()

# Change column names 
colnames(combined_df_wide) <- c("cond", 
                                "PNIE_LL", "PNIE_UL", "TNDE_LL", "TNDE_UL", 
                                "TNIE_LL", "TNIE_UL", "PNDE_LL", "PNDE_UL")

# Merge with existing results
results_DF_noRE <- merge(results_DF[1:6, ], combined_df_wide, by = "cond")

# Clean up environment (drop objects from this section that are no longer needed)
rm(lists, list_names, direct_df, indirect_df, combined_df, combined_df_wide)


#### Create dataframe of final estimates -----------------------------------

# Combine RE and non-RE results
results_DF <- rbind(results_DF_noRE, results_DF_RE)

# View the results
print(results_DF)
#   cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL     PNDE_UL
# 1 fefe -0.1889825 -0.1894340 -0.09570423 -0.0002112213 -0.2966453 0.10177927 -0.4418273  0.059211894 -0.3019212 0.10584214 -0.4397014  0.06128061
# 2 fesl -0.1973009 -0.1976276 -0.09752048 -0.0367624900 -0.2846837 0.09221394 -0.4424248  0.046017498 -0.2909517 0.09654751 -0.4436347  0.04782185
# 3 refe -0.1958658 -0.1963458 -0.09240229  0.0030725382 -0.2852170 0.09831936 -0.4478837  0.046145304 -0.2902472 0.10347047 -0.4487281  0.04482294
# 4 resl -0.2400019 -0.2401109 -0.07276088 -0.0045414414 -0.2526444 0.10399138 -0.4826027 -0.004350140 -0.2584021 0.10995042 -0.4826292 -0.00266291
# 5 slfe -0.1804319 -0.1808842 -0.09704172 -0.0363361536 -0.2882729 0.09692833 -0.4341171  0.062637519 -0.2928385 0.10242150 -0.4350416  0.06349128
# 6 slsl -0.2764533 -0.2764801 -0.03027051  0.0180582390 -0.2159972 0.15090354 -0.5105975 -0.029469066 -0.2171044 0.15043562 -0.5132218 -0.03095103
# 7 fere -0.1946794 -0.1951001 -0.09536315 -0.0126787393 -0.3103216 0.09548134 -0.4548479  0.045418542 -0.3232379 0.09920591 -0.3232379  0.09920591
# 8 rere -0.2171489 -0.2174313 -0.08069368  0.0039460448 -0.2753214 0.10296388 -0.4714517  0.012829715 -0.2980896 0.10930377 -0.2980896  0.10930377
# 9 slre -0.2208416 -0.2210517 -0.06717527 -0.0117677236 -0.2650267 0.11912664 -0.4662056 -0.001344533 -0.2763510 0.12138392 -0.2763510  0.12138392

# Add PS model & Mediator/Outcome model labels 
results_DF <- results_DF %>%
  mutate(
    PS = case_when(
      startsWith(cond, "fe") ~ "Fixed-Effect",
      startsWith(cond, "sl") ~ "Single-Level",
      startsWith(cond, "re") ~ "Random-Effect",
      TRUE ~ NA_character_  # Default case
    ),
    Model = case_when(
      endsWith(cond, "fe") ~ "Fixed-Effect",
      endsWith(cond, "sl") ~ "Single-Level",
      endsWith(cond, "re") ~ "Random-Effect",
      TRUE ~ NA_character_  # Default case
    )
  )

# Print the modified dataframe
print(results_DF)
#   cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL     PNDE_UL            PS         Model
# 1 fefe -0.1889825 -0.1894340 -0.09570423 -0.0002112213 -0.2966453 0.10177927 -0.4418273  0.059211894 -0.3019212 0.10584214 -0.4397014  0.06128061  Fixed-Effect  Fixed-Effect
# 2 fesl -0.1973009 -0.1976276 -0.09752048 -0.0367624900 -0.2846837 0.09221394 -0.4424248  0.046017498 -0.2909517 0.09654751 -0.4436347  0.04782185  Fixed-Effect  Single-Level
# 3 refe -0.1958658 -0.1963458 -0.09240229  0.0030725382 -0.2852170 0.09831936 -0.4478837  0.046145304 -0.2902472 0.10347047 -0.4487281  0.04482294 Random-Effect  Fixed-Effect
# 4 resl -0.2400019 -0.2401109 -0.07276088 -0.0045414414 -0.2526444 0.10399138 -0.4826027 -0.004350140 -0.2584021 0.10995042 -0.4826292 -0.00266291 Random-Effect  Single-Level
# 5 slfe -0.1804319 -0.1808842 -0.09704172 -0.0363361536 -0.2882729 0.09692833 -0.4341171  0.062637519 -0.2928385 0.10242150 -0.4350416  0.06349128  Single-Level  Fixed-Effect
# 6 slsl -0.2764533 -0.2764801 -0.03027051  0.0180582390 -0.2159972 0.15090354 -0.5105975 -0.029469066 -0.2171044 0.15043562 -0.5132218 -0.03095103  Single-Level  Single-Level
# 7 fere -0.1946794 -0.1951001 -0.09536315 -0.0126787393 -0.3103216 0.09548134 -0.4548479  0.045418542 -0.3232379 0.09920591 -0.3232379  0.09920591  Fixed-Effect Random-Effect
# 8 rere -0.2171489 -0.2174313 -0.08069368  0.0039460448 -0.2753214 0.10296388 -0.4714517  0.012829715 -0.2980896 0.10930377 -0.2980896  0.10930377 Random-Effect Random-Effect
# 9 slre -0.2208416 -0.2210517 -0.06717527 -0.0117677236 -0.2650267 0.11912664 -0.4662056 -0.001344533 -0.2763510 0.12138392 -0.2763510  0.12138392  Single-Level Random-Effect

# Save dataframe of results 
write_rds(results_DF, file = "Application/Output/Estimates.rds")


# Result visuals ----------------------------------------------------------

# TNDE 
## Save visual 
pdf("Application/Output/TNDE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(TNDE_LL > 0 | TNDE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = TNDE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Total Natural Direct Effect (TNDE) with 95% Confidence Intervals",
       x = "Total Natural Direct Effect (TNDE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/TNDE-Estimates.png", plot = last_plot())


# PNDE 
## Save visual 
pdf("Application/Output/PNDE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(PNDE_LL > 0 | PNDE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = PNDE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Pure Natural Direct Effect (PNDE) with 95% Confidence Intervals",
       x = "Pure Natural Direct Effect (PNDE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/PNDE-Estimates.png", plot = last_plot())




# TNIE 
## Save visual 
pdf("Application/Output/TNIE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(TNIE_LL > 0 | TNIE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = TNIE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Total Natural Indirect Effect (TNIE) with 95% Confidence Intervals",
       x = "Total Natural Indirect Effect (TNIE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/TNIE-Estimates.png", plot = last_plot())


# PNIE 
## Save visual 
pdf("Application/Output/PNIE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(PNIE_LL > 0 | PNIE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = PNIE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Pure Natural Indirect Effect (PNIE) with 95% Confidence Intervals",
       x = "Pure Natural Indirect Effect (PNIE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/PNIE-Estimates.png", plot = last_plot())




# Result Table ------------------------------------------------------------

results_DF[, c("PS", "Model", 
               "TNDE", "TNDE_LL", "TNDE_UL", 
               "PNDE", "PNDE_LL", "PNDE_UL", 
               "TNIE", "TNIE_LL", "TNIE_UL", 
               "PNIE", "PNIE_LL", "PNIE_UL")]


