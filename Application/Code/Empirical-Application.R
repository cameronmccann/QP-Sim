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
# Last Updated: 07/24/2024 
#
#
# Notes:
#   To-Do:
#     + Refer to "NEXT" area in code
#   Then: 
#     + Check & clean up comments for the data cleaning sections (Wave I, III, IV)
#     + Also double check, delete & clean the following sections: Drop small clusters, Missing Data
#     + Add visuals of covariate balance & overlap 
# 
# 

# 
#     + MAYBE go back & create viz to check overlap
#     + Figure out how to get covariate balance for all models on 1 visual (like Chang et al. article) 
# 
#     + increase number of clusters to all of those greater or equal to 5 (Liu email)
#       - write RE mediator & outcome code 
#       - create bootstrap function
#       - address model convergence issues: with PS (& maybe mediator & outcome) 
#     + I had model convergence issues with 15 clusters (cs: 41-81 students; N = 742). So I will attempt to increase the number of clusters 
# 
#     + Obtain variables names from codebook
#       - School ID
# 
#   Next: 
#     + 
#     + Figure out how to get covariate balance for all models on 1 visual (like Chang et al. article) 
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
source("Application/Functions/bootstrap_ci_paral.R")
source("Application/Functions/bootstrap_ci_re_paral.R")




# Import & Clean Data -----------------------------------------------------

## Wave I ------------------------------------------------------------------

# Import student data
w1 <- readr::read_tsv(file = "Application/Data/ICPSR_21600/DS0001/21600-0001-Data.tsv")

# Drop variables in wave I student data
w1 <- w1 %>% 
  select("AID":"SCH_YR", "S1":"S7", "S12", "S18", "S11", "S17", "PC22",
         "H1FS1":"H1FS19", # Feelings Scale (contains mostly CES-D items)
         "H1NF12B", "S44A14", "S44A18":"S44A29")

# In School Questionnaire Code Book Public Use Sample (pg. 697 of Questionnaire.pdf)
# age (S1)
# sex (S2)
# race/ethnicity (S4 - S7)
data <- w1 %>% 
  rename(age = S1, 
         sex = S2, 
         ethnicity = S4, # Create ethnicity variable (indicates Hispanic or Spanish origin)
         white = S6A, 
         black = S6B, 
         asian = S6C, 
         nativeAmerican = S6D, 
         raceOther = S6E)

# will those indicating an 8 (for S4 or ethnicity variable) mess things up? (same with sex)
summary(data$sex) #summary(data$S2)
summary(data$ethnicity) #summary(data$S4)
# Change 8 to NA
data$ethnicity <- ifelse(data$ethnicity == 8, NA, data$ethnicity)
# sum(I(data$ethnicity == 8), na.rm = TRUE)

# parental education 
## mother - How far in school did she go? (S12)
## father - How far in school did he go? (S18)
### [recode this to --> parental education (high school or less, some college, and college graduate or more)]

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

data$parentalEdu <- apply(data[, stringr::str_detect(colnames(data), pattern = "Edu$")], 1, function(x) max(x, na.rm = FALSE))

# data %>% 
#   select("dadEdu", "momEdu", "parentalEdu") %>% 
#   head()


# family structure 
## mother - Do you live with your biological mother, stepmother, foster mother, or adoptive mother? (S11)
## father - Do you live with your biological father, stepfather, foster father, or adoptive father? (S17)
### [recode this to --> family structure (2 parents, single mom, and single dad or other)]


data$familyStruct <- data$S11 + data$S17

table(data$S11); sum(is.na(data$S11))
table(data$S17); sum(is.na(data$S17))
table(data$S11, data$S17)
data$familyStruct <- ifelse(data$familyStruct == 10, NA, data$familyStruct) # change "multiple response" to missing 

table(data$familyStruct)

# gap in health insurance (wave 1) 
## In the past 12 months, has there been a time when {NAME} had no health insurance? (PC22)
data <-
  data %>%
  rename(healthInsur = PC22) %>%
  mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know; 99=not applicable


# Sports participation 
## Have you played a sport? (H1NF12B)
## Here is a list of clubs, organizations, and teams found at many schools. Darken the oval next to any of them that you are participating in this year, or that you plan to participate in later in the school year.
## (S44A14 & S44A18 - S44A29)

# sum number of sports 
data$sport <-
  rowSums(data[, colnames(data)[stringr::str_starts(colnames(data), pattern = "^S44A")]])
# Dichotomize 
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
         "healthInsur", "parentalEdu", "familyStruct", "sport", "sportPartic", "feelings")

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

# gap in health insurance (wave 3) 
# 1. Over the past 12 months, how many months did you have health insurance? (H3HS1)
w3 <- w3 %>%
  rename(healthInsur = H3HS1) %>%
  mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know; 99=not applicable

# Add wave III student data to data set
data <- merge(data, w3[, c("AID", "selfEst", "healthInsur")], by = "AID")

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
         "H4MH18":"H4MH27", # CES-D-10 scale [NOTE: DOUBLE CHECK THESE ARE CORRECT ITEMS]
         "H4HS3") # Health insurance item 


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

# check 
hist(w4$depress)
table(I(w4$depress > 10)); sum(I(w4$depress > 10), na.rm = T) / sum(!is.na(w4$depress)) # 16% would be classified with depression in Easterlin study


# gap in health insurance (wave 4) 
# 3. Over the past 12 months, how many months did you have health insurance? (H4HS3)
w4 <- w4 %>%
  rename(healthInsur = H4HS3) %>%
  mutate(healthInsur = ifelse(healthInsur >= 98, NA, healthInsur)) # 98=don't know

# Add wave IV student data to data set
data <- merge(data, w4[, c("AID", "healthInsur", "depress")], by = "AID")

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
data$healthInsur_w3 <- ifelse(data$healthInsur_w3 > 0, 1, 0)
data$healthInsur_w4 <- ifelse(data$healthInsur_w4 > 0, 1, 0)
data$healthInsur_gap <- ifelse(data$healthInsur_w1 == 0 | data$healthInsur_w3 == 0 | data$healthInsur_w4 == 0,
                               1, 0)

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
  "healthInsur_w3",
  "healthInsur_w4",
  "depress_w4"
)])

# PS model includes: "sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 + parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2)"

# Need to scale: 
# Scale, rename, & drop non-scaled variable
data <- data %>%
  mutate(
    age_w1_sc = as.vector(scale(age_w1)),
    parentalEdu_w1_sc = as.vector(scale(parentalEdu_w1)),
    feelings_w1_sc = as.vector(scale(feelings_w1))
  ) %>%
  select(-age_w1, -parentalEdu_w1, -feelings_w1)  # Remove original columns if desired

# Select variables in PS model 
t <- data %>% 
  select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_w3, healthInsur_w4, healthInsur_w1)) %>% 
  # select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_gap))
  # healthInsur_w3, healthInsur_w4, 
  # healthInsur_w1))
  select(c("sex_w1", "white_w1", "black_w1", 
           "sportPartic_w1", "selfEst_w3", "depress_w4", 
           "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_gap", 
           "age_w1_sc", "feelings_w1_sc"))

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



# PS & IPTW Calculation ---------------------------------------------------------------
# This script calculates propensity scores and Inverse Probability of Treatment Weights (IPTW)
# for three models: SL (Standard Logistic), FE (Fixed Effects), and RE (Random Effects).

### SL Propensity Score Model -----------------------------------------------
# Standard logistic regression model for propensity score estimation
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
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
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap +
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
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap +
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



# [INSERT VISUALS] --------------------------------------------------------

# ADD VISUALS TO SHOW HISTOGRAM OR DENSITY OVERLAP 
# ADD VISUAL SHOWING DOT PLOT ACROSS MODELS ON EACH COVARIATE OF SAMPLE BALANCE (SMD LIKE LOVEPLOT)



# Estimate Effects --------------------------------------------------------

## Mediator models ---------------------------------------------------------

#### Single-Level (SL) Models ------------------------------------------------

# Model: SL PS & SL - Mediation/Outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_sl
  )

# Model: FE PS & SL - Mediation/Outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_fe
  )

# Model: RE PS & SL - Mediation/Outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_re
  )


#### Fixed-Effect (FE) Models ------------------------------------------------

# Model: SL PS & FE - Mediation/Outcome
med_slfe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_sl
)

# Model: FE PS & FE - Mediation/Outcome 
med_fefe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_fe
)

# Model: RE PS & FE - Mediation/Outcome
med_refe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
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
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )

# Model: FE PS & RE - Mediation/Outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )

# Model: RE PS & RE - Mediation/Outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )



## Outcome Models ----------------------------------------------------------

#### Single-Level (SL) Models ------------------------------------------------

# Model: SL PS & SL Mediation/Outcome
out_slsl <-
  glm(
    formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_sl
  )

# Model: FE PS & SL Mediation/Outcome
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_fe
  )

# Model: RE PS & SL Mediation/Outcome
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_re
  )

#### Fixed-Effect (FE) Models ------------------------------------------------

# Model: SL PS & FE Mediation/Outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)

# Model: FE PS & FE Mediation/Outcome
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Model: RE PS & FE Mediation/Outcome
out_refe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
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
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )

# Model: FE PS & RE Mediation/Outcome
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )

# Model: RE PS & RE Mediation/Outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


# Display Estimates -------------------------------------------------------

# Define conditions
conditions <- c("slsl", "fesl", "resl", "slfe", "fefe", "refe", "slre", "fere", "rere")

# Extract NDE estimates
NDE <- sapply(conditions, function(cond) {
  model_name <- paste0("out_", cond)
  summary(get(model_name))$coef["sportPartic_w1", "Estimate"]
})

# Extract NIE estimates
NIE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond)
  summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"] * 
    summary(get(out_model_name))$coef["selfEst_w3_sc", "Estimate"]
})

# Create results DataFrame
results_DF <- data.frame(
  cond = conditions,
  NDE = NDE,
  NIE = NIE
)

# Display results
rownames(results_DF) <- NULL
results_DF

# cond        NDE        NIE
# 1 slsl -0.3279340 -0.1349289
# 2 fesl -0.2337515 -0.2072494
# 3 resl -0.2772344 -0.1727740
# 4 slfe -0.2024146 -0.1935592
# 5 fefe -0.2212818 -0.2102262
# 6 refe -0.2245736 -0.1958965
# 7 slre -0.2510726 -0.1651808
# 8 fere -0.2279858 -0.2074429
# 9 rere -0.2482880 -0.1817158






# Conduct Bootstrap CI ------------------------------------------------------------

# SL Mediator/Outcome
# SL PS Model 
slsl_ci <- bootstrap_ci_paral(iterations = 1000,
                   iptw = iptw_sl,
                   data = data,
                   model = "SL",
                   cores = 6,
                   core_seeds = c(4561:4566))
# FE PS Model 
fesl_ci <- bootstrap_ci_paral(iterations = 1000,
                              iptw = iptw_fe,
                              data = data,
                              model = "SL",
                              cores = 6,
                              core_seeds = c(4561:4566))
# RE PS Model 
resl_ci <- bootstrap_ci_paral(iterations = 1000,
                              iptw = iptw_re,
                              data = data,
                              model = "SL",
                              cores = 6,
                              core_seeds = c(4561:4566))

# FE Mediator/Outcome 
# SL PS Model 
slfe_ci <- bootstrap_ci_paral(iterations = 1000,
                              iptw = iptw_sl,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566))
# FE PS Model 
fefe_ci <- bootstrap_ci_paral(iterations = 1000,
                              iptw = iptw_fe,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566))
# RE PS Model 
refe_ci <- bootstrap_ci_paral(iterations = 1000,
                              iptw = iptw_re,
                              data = data,
                              model = "FE",
                              cores = 6,
                              core_seeds = c(4561:4566))

# RE Mediator/Outcome 
# SL PS Model 

execution_time <- system.time({ # track computation time 
slre_ci <- bootstrap_ci_re_paral(iterations = 1500, 
                                 iptw = iptw_sl, 
                                 data = data, 
                                 cores = 6, 
                                 core_seeds = c(4561:4566))
})
# Print the execution time
print(execution_time)

# If you want to see the elapsed time specifically
cat("Elapsed time:", execution_time["elapsed"], "seconds\n")
# user   system  elapsed 
# 4062.031   70.951  872.969 
# Elapsed time: 872.969 seconds

# 
paste0("Number of converged mediator models: ", slre_ci$mediator_converged_count, " (", (slre_ci$mediator_converged_count/length(slre_ci$direct_effects))*100, "%)")
paste0("Number of converged outcome models: ", slre_ci$outcome_converged_count, " (", (slre_ci$outcome_converged_count/length(slre_ci$direct_effects))*100, "%)")
paste0("Number of iterations with both models converged: ", slre_ci$both_converged_count, " (", (slre_ci$both_converged_count/length(slre_ci$direct_effects))*100, "%)")
# [1] "Number of converged mediator models: 765 (76.5%)"
# [1] "Number of converged outcome models: 998 (99.8%)"
# [1] "Number of iterations with both models converged: 765 (76.5%)"





##################
# NEXT
#   + Check time took to compute & how many iterations converged for slre_ci, then up the iterations (maybe something like 1,500)
#   
#   + proceed to repeat with FE & RE PS models 
#   + Check that you have 1,000 good/converged iterations
#   + grab first 1,000 good/converged (in)direct effect estimates & compute percentile bootstrap CI (with quantile())
#       + Store this in a dataframe with columns & column names like results_DF_noRE (cond, NIE_LL, NIE_UL, NDE_LL, NDE_UL)
#       + Then merge this with their respective point estimates 

# 
# 


# FE PS Model 

# RE PS Model 







# Join results ------------------------------------------------------------



## Join CI and point estimates for SL & FE med/outcome  --------------------

# List of your lists
list_names <- c("slsl_ci", "fesl_ci", "resl_ci", "slfe_ci", "fefe_ci", "refe_ci")
lists <- list(slsl_ci, fesl_ci, resl_ci, slfe_ci, fefe_ci, refe_ci)

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

# combine both into a single data frame
combined_df <- merge(indirect_df, direct_df, by = "list_name")
# print(combined_df)

# 
combined_df$list_name <- sub("_ci$", "", combined_df$list_name)

# Rename columns starting with "indirect_ci" to "NIE"
colnames(combined_df) <- gsub("^indirect_ci", "NIE", colnames(combined_df))

# Rename columns starting with "direct_ci" to "NDE"
colnames(combined_df) <- gsub("^direct_ci", "NDE", colnames(combined_df))

results_DF_noRE <- merge(results_DF[1:6, ], combined_df, by.x = "cond", by.y = "list_name")




## Join CI and point estimates for RE med/outcome --------------------------










