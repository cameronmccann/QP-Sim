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
# Last Updated: 07/17/2024 
#
#
# Notes:
#   To-Do:
#     + Drop nonoverlaping cases 
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
  WeMix
  # ggdag, 
  # dagitty, 
  # huxtable
)




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
data$familyStruct <- ifelse(data$familyStruct == 9, NA, data$familyStruct) # change "multiple response" to missing 

# table(data$familyStruct)

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

# Add cluster (school) size variable 
data <- data %>% 
  group_by(CLUSTER2) %>% 
  mutate(n = n()) %>% 
  ungroup()



# Missing Data ------------------------------------------------------------

# Drop those missing sport participation (treatment)
data <- data %>% 
  filter(is.na(sportPartic_w1) != TRUE)

# Select variables in PS model 
t <- data %>% 
  select(!c(AID, CLUSTER2, sport_w1))

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



# Descriptives & choose cluster size  -------------------------------------

# Calculate cluster sizes again 
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


# Drop cluster sizes < 41 
# data <- data %>%
#   filter(n > 40) # To avoid convergence issues, we restricted our analysis to the 15 largest schools
# NOTE: I had convergence issues with 15 clusters, so I am trying more

# Drop cluster sizes < 5
data <- data %>%
  filter(n >= 5) # To avoid convergence issues, we restricted our analysis to the schools with 5 or more students (121 schools)


# cluster sizes
data %>% 
  group_by(CLUSTER2) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n = 132) # 15 schools (cluster sizes ranged from 41 to 81 students)
                 # 121 schools (cluster sizes ranged from 5 to 81 students)

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
med_icc # 0.02 for med icc
        # with 121 schools (5+ in size) med icc = 0.007111455

# ICC for outcome 
out_unconditional <- lme4::lmer(depress_w4 ~ (1 | CLUSTER2), data = data)
summary(out_unconditional)

out_var <- data.frame(lme4::VarCorr(out_unconditional))[1, 4] # variance due to level-2
out_res <- data.frame(lme4::VarCorr(out_unconditional))[2, 4] # residual variance 
out_icc <- out_var / (out_var + out_res)
out_icc # 0.016 for outcome icc
        # with 121 schools (5+ in size) outcome icc = 0.01866436



# PS ----------------------------------------------------------------------

## SL PS model -------------------------------------------------------------
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + 
                ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                parentalEdu_w1 + familyStruct_w1 + healthInsur_w1",
                family = "binomial", 
                data = data)
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))


## FE PS model -------------------------------------------------------------
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + 
                ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))


## RE PS model -------------------------------------------------------------
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + 
                        ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                        parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data, 
                        control = lme4::glmerControl(optimizer = "bobyqa")) # Changing optimizer
# control = glmerControl(optCtrl = list(maxfun = 100000))) # Increase max iterations to 1000
# control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Increase max iterations to 1000
data$ps_re <- predict(psmod_re, type = "response")
data$ps_re_logit <- predict(psmod_re, type = "link")
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))



# Viz for balancing PS variables ------------------------------------------

# Display density plot for covariate balance before and after PS weighting 
# SL
cobalt::bal.plot(
  WeightIt::weightit(
    formula = sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 +
      ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w1,
    data = data,
    estimand = "ATE",
    method = "ps"
  ),
  which = "both",
  type = "density"
)
# FE PS 
cobalt::bal.plot(
  WeightIt::weightit(
    formula = sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 +
      ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 +
      as.factor(CLUSTER2),
    data = data,
    estimand = "ATE",
    method = "ps"
  ),
  which = "both",
  type = "density"
)
# RE PS 
# cobalt::bal.plot(
#   WeightIt::weightit(
#     formula = sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 +
#       ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
#       parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2),
#     data = data,
#     estimand = "ATE",
#     method = "ps"
#   ),
#   which = "both",
#   type = "density"
# )

# Dot plot of Absolute Std Mean Diff
# SL PS 
cobalt::bal.tab(sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 +
                  ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                  parentalEdu_w1 + familyStruct_w1 + healthInsur_w1,
                data = data, 
                weights = data$iptw_sl, 
                un = TRUE) %>% 
  cobalt::love.plot(stars = "raw", 
                    grid = TRUE, 
                    abs = TRUE, 
                    var.order = "unadjusted", 
                    alpha = 0.8, 
                    thresholds = 0.1, 
                    line = FALSE)
# FE PS 
cobalt::bal.tab(sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + 
                  ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                  parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 +
                  as.factor(CLUSTER2),
                data = data, 
                weights = data$iptw_fe, 
                un = TRUE) %>% 
  cobalt::love.plot(stars = "raw", 
                    grid = TRUE, 
                    abs = TRUE, 
                    var.order = "unadjusted", 
                    alpha = 0.8, 
                    thresholds = 0.1, 
                    line = FALSE)
# RE PS 
cobalt::bal.tab(sportPartic_w1 ~ feelings_w1 + age_w1 + sex_w1 + 
                  ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                  parentalEdu_w1 + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2),
                data = data, 
                weights = data$iptw_re, 
                un = TRUE) %>% 
  cobalt::love.plot(stars = "raw", 
                    grid = TRUE, 
                    abs = TRUE, 
                    var.order = "unadjusted", 
                    alpha = 0.8, 
                    thresholds = 0.1, 
                    line = FALSE)




# Drop Nonoverlap ---------------------------------------------------------

## RE ----------------------------------------------------------------------
### fixed threshold ---------------------------------------------------------
# Identify nonoverlapping cases with a caliper of 0.05
caliper <- 0.05
data_re <- data
data_re <- data_re %>% 
  mutate(
    nonoverlap = case_when(
      sportPartic_w1 = 1 & (ps_re < caliper | ps_re > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_re < caliper | ps_re > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_left = case_when(
      sportPartic_w1 = 1 & (ps_re < caliper) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_re < caliper) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_right = case_when(
      sportPartic_w1 = 1 & (ps_re > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_re > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

paste0(sum(data_re$nonoverlap), " non-overlapping cases (left: ", 
       sum(data_re$nonoverlap_left), "; right: ", sum(data_re$nonoverlap_right), ")")


### dynamic threshold -------------------------------------------------------
caliper <- 0.05
data_re <- data
data_re <- data_re %>% 
  mutate(
    overlap_left = max(tapply(ps_re_logit, sportPartic_w1, min)) - sd(ps_re_logit) * caliper,
    overlap_right = min(tapply(ps_re_logit, sportPartic_w1, max)) + sd(ps_re_logit) * caliper, 
    # ind_nonoverlap = case_when(
    #   ps_re_logit > overlap_right ~ TRUE,
    #   ps_re_logit < overlap_left ~ TRUE,
    #   TRUE ~ FALSE
    # ), 
    ind_nonoverlap_left = case_when(
      ps_re_logit < overlap_left ~ TRUE,
      TRUE ~ FALSE
    ), 
    ind_nonoverlap_right = case_when(
      ps_re_logit > overlap_right ~TRUE, 
      TRUE ~ FALSE
    ), 
    ind_nonoverlap = case_when(
      (ind_nonoverlap_left == TRUE | ind_nonoverlap_right == TRUE) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

paste0("Number of non-overlapping cases: ", sum(data_re$ind_nonoverlap), " (left: ", sum(data_re$ind_nonoverlap_left), "; right: ", sum(data_re$ind_nonoverlap_right), ")")



# # drop nonoverlap cases 
# data_re <- data_re[data_re$ind_nonoverlap == FALSE, ]


# SL 
# Identify nonoverlapping cases with a caliper of 0.05
caliper <- 0.05
data_sl <- data
data_sl <- data_sl %>% 
  mutate(
    nonoverlap = case_when(
      sportPartic_w1 = 1 & (ps_sl < caliper | ps_sl > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_sl < caliper | ps_sl > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_left = case_when(
      sportPartic_w1 = 1 & (ps_sl < caliper) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_sl < caliper) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_right = case_when(
      sportPartic_w1 = 1 & (ps_sl > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_sl > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

sum(data_sl$nonoverlap)
sum(data_sl$nonoverlap_left)
sum(data_sl$nonoverlap_right)

# FE 
# Identify nonoverlapping cases with a caliper of 0.05
caliper <- 0.05
data_fe <- data
data_fe <- data_fe %>% 
  mutate(
    nonoverlap = case_when(
      sportPartic_w1 = 1 & (ps_fe < caliper | ps_fe > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_fe < caliper | ps_fe > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_left = case_when(
      sportPartic_w1 = 1 & (ps_fe < caliper) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_fe < caliper) ~ TRUE, 
      TRUE ~ FALSE
    ), 
    nonoverlap_right = case_when(
      sportPartic_w1 = 1 & (ps_fe > (1 - caliper)) ~ TRUE, 
      sportPartic_w1 = 0 & (ps_fe > (1 - caliper)) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

sum(data_fe$nonoverlap)
sum(data_fe$nonoverlap_left)
sum(data_fe$nonoverlap_right)









# Estimate Effects --------------------------------------------------------


## Mediator models ---------------------------------------------------------

### Mediator models (SL) ----------------------------------------------------

# Single-Level PS & SL med/outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect PS & SL med/outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_fe
  )
# Random-Effect PS & SL med/outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_re
  )


### Mediator models (FE) --------------------------------------------------

# Single-Level PS & FE med/outcome
med_slfe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
med_fefe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Random-Effect PS & FE med/outcome
med_refe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_re)


### Mediator models (RE) ----------------------------------------------------

### Add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Single-Level PS & RE med/outcome
med_slre <-
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )



## Outcome models ----------------------------------------------------------

### Outcome models (SL) -----------------------------------------------------

# Single-Level
out_slsl <-
  glm(
    formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_fe
  )
# Random-Effect
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_re
  )


### Outcome models (FE) -----------------------------------------------------

# Single-Level PS & FE med/outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)
# Random-Effect PS & FE med/outcome
out_refe <-glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)",
      data = data,
      weights = iptw_re
)


### Outcome models (RE) -----------------------------------------------------

# ### If not done so with mediator model, add column of just 1s for level-2 weight
# data <- cbind(data, L2weight = rep(1, nrow(data)))

# Single-Level PS & RE med/outcome
out_slre <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


## Display Estimates -------------------------------------------------------


# Add NDE & NIE estimates 
results_DF <- data.frame(
  cond = c("slsl", "fesl", "resl",
           "slfe", "fefe", "refe", 
           "slre", "fere", "rere"),
  # est = rep("NDE", 3),
  NDE = c(
    summary(out_slsl)$coef["sportPartic_w1", "Estimate"],
    summary(out_fesl)$coef["sportPartic_w1", "Estimate"],
    summary(out_resl)$coef["sportPartic_w1", "Estimate"],
    
    summary(out_slfe)$coef["sportPartic_w1", "Estimate"],
    summary(out_fefe)$coef["sportPartic_w1", "Estimate"],
    summary(out_refe)$coef["sportPartic_w1", "Estimate"],
    
    summary(out_slre)$coef["sportPartic_w1", "Estimate"], 
    summary(out_fere)$coef["sportPartic_w1", "Estimate"], 
    summary(out_rere)$coef["sportPartic_w1", "Estimate"]
  ),
  
  NIE = c(
    summary(med_slsl)$coef["sportPartic_w1", "Estimate"] * summary(out_slsl)$coef["selfEst_w3", "Estimate"],
    summary(med_fesl)$coef["sportPartic_w1", "Estimate"] * summary(out_fesl)$coef["selfEst_w3", "Estimate"],
    summary(med_resl)$coef["sportPartic_w1", "Estimate"] * summary(out_resl)$coef["selfEst_w3", "Estimate"],
    
    summary(med_slfe)$coef["sportPartic_w1", "Estimate"] * summary(out_slfe)$coef["selfEst_w3", "Estimate"],
    summary(med_fefe)$coef["sportPartic_w1", "Estimate"] * summary(out_fefe)$coef["selfEst_w3", "Estimate"],
    summary(med_refe)$coef["sportPartic_w1", "Estimate"] * summary(out_refe)$coef["selfEst_w3", "Estimate"], 
    
    summary(med_slre)$coef["sportPartic_w1", "Estimate"] * summary(out_slre)$coef["selfEst_w3", "Estimate"], 
    summary(med_fere)$coef["sportPartic_w1", "Estimate"] * summary(out_fere)$coef["selfEst_w3", "Estimate"], 
    summary(med_rere)$coef["sportPartic_w1", "Estimate"] * summary(out_rere)$coef["selfEst_w3", "Estimate"]
  )
)

# results_DF
#   cond       NDE         NIE
# 1 slsl 0.4034979 -0.01127713  # With 15 largest schools 
# 2 fesl 0.5285937  0.07458896
# 3 resl 0.5139536  0.03039552
# 4 slfe 0.5049625 -0.05323353
# 5 fefe 0.5015928  0.05567506
# 6 refe 0.5086497  0.02351801
# 7 slre 0.4578713 -0.03441880
# 8 fere 0.5097119  0.06037420
# 9 rere 0.5099829  0.02680802

#   cond        NDE         NIE
# 1 slsl -0.2368474 -0.02117536
# 2 fesl -0.1184377 -0.03255229
# 3 resl -0.1793095 -0.02837372
# 4 slfe -0.1407412 -0.04280035
# 5 fefe -0.1439630 -0.03339718
# 6 refe -0.1445497 -0.03370138
# 7 slre -0.1783744 -0.03331416
# 8 fere -0.1399102 -0.03423868
# 9 rere -0.1611487 -0.03028568








