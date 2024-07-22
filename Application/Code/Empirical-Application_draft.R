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
# Last Updated: 07/22/2024 
#
#
# Notes:
#   To-Do:
#     + Standardize variables, impute missing, data, & check model collinearity & proportion of students in each demographic group (likely add race categories: other, asian, native american)
#     + Drop nonoverlaping cases 
#     + Med & Outcome model formulas 
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

# Drop cluster sizes < 5
data <- data %>%
  filter(n >= 5) # To avoid convergence issues, we restricted our analysis to the schools with 5 or more students (121 schools)

# Drop cluster sizes < 41 
# data <- data %>%
#   filter(n > 40) # To avoid convergence issues, we restricted our analysis to the 15 largest schools
# NOTE: I had convergence issues with 15 clusters, so I am trying more


# Drop cluster sizes < 25
# data <- data %>%
#   filter(n >= 40) # To avoid convergence issues, we restricted our analysis to the schools with 25 or more students (62 schools)


# cluster sizes
data %>% 
  group_by(CLUSTER2) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n = 132) # 15 schools (cluster sizes ranged from 41 to 81 students)
                # 121 schools (cluster sizes ranged from 5 to 81 students)
                # 62 schools (cluster sizes ranged from 25 to 81 students)
                # 40 schools (cluster sizes ranged from 40 to 81 students)



# Missing Data ------------------------------------------------------------

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
  select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_w3, healthInsur_w4, healthInsur_w1))
  # select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_gap))
            # healthInsur_w3, healthInsur_w4, 
            # healthInsur_w1))

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




# Descriptives & choose cluster size  -------------------------------------

# # Calculate cluster sizes again 
# data <- data %>% 
#   group_by(CLUSTER2) %>% 
#   mutate(n = n()) %>% 
#   ungroup()
# 
# # Check clusters in order of size 
# data %>% 
#   group_by(CLUSTER2) %>% 
#   summarize(n = max(n)) %>% 
#   arrange(desc(n)) %>% 
#   print(n = 132) # 123 schools (cluster sizes ranged from 1 to 81 students)
# 
# 
# # Drop cluster sizes < 41 
# # data <- data %>%
# #   filter(n > 40) # To avoid convergence issues, we restricted our analysis to the 15 largest schools
# # NOTE: I had convergence issues with 15 clusters, so I am trying more
# 
# # Drop cluster sizes < 5
# # data <- data %>%
# #   filter(n >= 5) # To avoid convergence issues, we restricted our analysis to the schools with 5 or more students (121 schools)
# 
# # Drop cluster sizes < 25
# data <- data %>%
#   filter(n >= 40) # To avoid convergence issues, we restricted our analysis to the schools with 25 or more students (62 schools)
# 
# 
# # cluster sizes
# data %>% 
#   group_by(CLUSTER2) %>% 
#   summarize(n = n()) %>% 
#   arrange(desc(n)) %>% 
#   print(n = 132) # 15 schools (cluster sizes ranged from 41 to 81 students)
#                  # 121 schools (cluster sizes ranged from 5 to 81 students)
#                  # 62 schools (cluster sizes ranged from 25 to 81 students)
#                  # 40 schools (cluster sizes ranged from 40 to 81 students)

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
        # with 121 schools (5+ in size) med icc = 0.007176327
        # with 62 schools (25+ in size) med icc = 0.01020648
        # with 17 schools (40+ in size) med icc = 0.02445009

# ICC for outcome 
out_unconditional <- lme4::lmer(depress_w4 ~ (1 | CLUSTER2), data = data)
summary(out_unconditional)

out_var <- data.frame(lme4::VarCorr(out_unconditional))[1, 4] # variance due to level-2
out_res <- data.frame(lme4::VarCorr(out_unconditional))[2, 4] # residual variance 
out_icc <- out_var / (out_var + out_res)
out_icc # 0.016 for outcome icc
        # with 121 schools (5+ in size) outcome icc = 0.0185748
        # with 62 schools (25+ in size) outcome icc = 0.01981785
        # with 17 schools (40+ in size) outcome icc = 0.01729325


# PS ----------------------------------------------------------------------

## SL PS model -------------------------------------------------------------
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1",
                family = "binomial", 
                data = data)
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")
# data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))


## FE PS model -------------------------------------------------------------
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
# data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))


## RE PS model -------------------------------------------------------------
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                        ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                        parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) 
                        # control = lme4::glmerControl(optimizer = "bobyqa")) # Changing optimizer
# control = glmerControl(optCtrl = list(maxfun = 100000))) # Increase max iterations to 1000
# control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Increase max iterations to 1000
data$ps_re <- predict(psmod_re, type = "response")
data$ps_re_logit <- predict(psmod_re, type = "link")
# data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))



# Potential collinearity issues ------------------------------------------------------------

# Check VIF
performance::check_collinearity(psmod_sl)
performance::check_collinearity(psmod_fe)
performance::check_collinearity(psmod_re)

# Racial categories frequency & corr 
paste0("Number of white: ", prettyNum(sum(data$white_w1), big.mark = ","), " (", round((sum(data$white_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of black: ", prettyNum(sum(data$black_w1), big.mark = ","), " (", round((sum(data$black_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of asian: ", prettyNum(sum(data$asian_w1), big.mark = ","), " (", round((sum(data$asian_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of nativeAm: ", prettyNum(sum(data$nativeAmerican_w1), big.mark = ","), " (", round((sum(data$nativeAmerican_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of other: ", prettyNum(sum(data$raceOther_w1), big.mark = ","), " (", round((sum(data$raceOther_w1) / nrow(data)) * 100, 2), "%)")
 # "Number of white: 1,965 (62.36%)"
 # "Number of black: 804 (25.52%)"
 # "Number of asian: 147 (4.67%)"
 # "Number of nativeAm: 183 (5.81%)"
 # "Number of other: 215 (6.82%)"
round(cor(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", 
             "asian_w1", "nativeAmerican_w1", "raceOther_w1", 
             "sportPartic_w1")]), 
      digits = 3)
#                   sex_w1 ethnicity_w1 white_w1 black_w1 asian_w1 nativeAmerican_w1 raceOther_w1 sportPartic_w1
# sex_w1             1.000       -0.022   -0.025    0.049   -0.016             0.005        0.009         -0.097
# ethnicity_w1      -0.022        1.000   -0.322   -0.088    0.014             0.048        0.382         -0.011
# white_w1          -0.025       -0.322    1.000   -0.686   -0.220            -0.062       -0.177          0.013
# black_w1           0.049       -0.088   -0.686    1.000   -0.040            -0.005       -0.101         -0.003
# asian_w1          -0.016        0.014   -0.220   -0.040    1.000            -0.003        0.030         -0.021
# nativeAmerican_w1  0.005        0.048   -0.062   -0.005   -0.003             1.000        0.078         -0.022
# raceOther_w1       0.009        0.382   -0.177   -0.101    0.030             0.078        1.000          0.004
# sportPartic_w1    -0.097       -0.011    0.013   -0.003   -0.021            -0.022        0.004          1.000

corrplot::corrplot(cor(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", "asian_w1", "nativeAmerican_w1", 
                                "raceOther_w1", "sportPartic_w1")]),
                   method = c("ellipse"))

PerformanceAnalytics::chart.Correlation(data[, c("ethnicity_w1", "white_w1", "black_w1", 
                                                 "asian_w1", "nativeAmerican_w1", 
                                                 "raceOther_w1", "sportPartic_w1")])

# Generate frequency matrix for racial categories 
dat <- as.data.frame(data[, c("white_w1", "black_w1", "asian_w1", "nativeAmerican_w1", "raceOther_w1")])
vars <- names(dat)
total_samples <- nrow(dat)
# Create an empty matrix to store results, with an extra row and column for totals
result_matrix <- matrix(NA, nrow = length(vars) + 1, ncol = length(vars) + 1, 
                        dimnames = list(c(vars, "Total"), c(vars, "Total")))
# Fill the matrix
for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    if (i == j) {
      # Diagonal: count of 1s in each variable
      count <- sum(dat[[vars[i]]] == 1)
    } else {
      # Off-diagonal: intersection count
      count <- sum(dat[[vars[i]]] == 1 & dat[[vars[j]]] == 1)
    }
    percentage <- round(count / total_samples * 100, 2)
    result_matrix[i, j] <- paste0(count, " (", percentage, "%)")
  }
}
# Calculate row totals
for (i in seq_along(vars)) {
  count <- sum(dat[[vars[i]]] == 1)
  percentage <- round(count / total_samples * 100, 2)
  result_matrix[i, length(vars) + 1] <- paste0(count, " (", percentage, "%)")
}
# Calculate column totals
for (j in seq_along(vars)) {
  count <- sum(dat[[vars[j]]] == 1)
  percentage <- round(count / total_samples * 100, 2)
  result_matrix[length(vars) + 1, j] <- paste0(count, " (", percentage, "%)")
}
# Set the overall total
overall_count <- sum(dat == 1)
overall_percentage <- round(overall_count / (total_samples * length(vars)) * 100, 2)
result_matrix[length(vars) + 1, length(vars) + 1] <- paste0(overall_count, " (", overall_percentage, "%)")
# Convert to data frame
result_df <- as.data.frame(result_matrix)
# Print the result
print(result_df)
#                       white_w1     black_w1    asian_w1 nativeAmerican_w1 raceOther_w1         Total
# white_w1          1965 (62.36%)   45 (1.43%)  21 (0.67%)        92 (2.92%)   66 (2.09%) 1965 (62.36%)
# black_w1             45 (1.43%) 804 (25.52%)  26 (0.83%)        45 (1.43%)   20 (0.63%)  804 (25.52%)
# asian_w1             21 (0.67%)   26 (0.83%) 147 (4.67%)         8 (0.25%)   15 (0.48%)   147 (4.67%)
# nativeAmerican_w1    92 (2.92%)   45 (1.43%)   8 (0.25%)       183 (5.81%)   27 (0.86%)   183 (5.81%)
# raceOther_w1         66 (2.09%)   20 (0.63%)  15 (0.48%)        27 (0.86%)  215 (6.82%)   215 (6.82%)
# Total             1965 (62.36%) 804 (25.52%) 147 (4.67%)       183 (5.81%)  215 (6.82%) 3314 (21.03%)


## RE PS model stepwise from intercept only model --------------------------

# Run RE PS model & add predictors 
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                        ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) # Model fails to converge when black_w1 &/or raceOther_w1 is added 
# Note: white_w1 corr with black_w1 (r = -0.69) & raceOther_w1 (r = -0.18); and ethnicity_w1 corr with raceOther_w1 (r = 0.38)



## Check black_w1 & raceOther_w1 relation with mediator & outcome ----------

round(cor(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", 
                   "asian_w1", "nativeAmerican_w1", "raceOther_w1", 
                   "sportPartic_w1", "selfEst_w3_sc", "depress_w4")]), 
      digits = 3)
#                   sex_w1 ethnicity_w1 white_w1 black_w1 asian_w1 nativeAmerican_w1 raceOther_w1 sportPartic_w1 selfEst_w3_sc depress_w4
# sex_w1             1.000       -0.022   -0.025    0.049   -0.016             0.005        0.009         -0.097        -0.079      0.125
# ethnicity_w1      -0.022        1.000   -0.322   -0.088    0.014             0.048        0.382         -0.011        -0.027      0.046
# white_w1          -0.025       -0.322    1.000   -0.686   -0.220            -0.062       -0.177          0.013        -0.031     -0.091
# black_w1           0.049       -0.088   -0.686    1.000   -0.040            -0.005       -0.101         -0.003         0.058      0.063
# asian_w1          -0.016        0.014   -0.220   -0.040    1.000            -0.003        0.030         -0.021        -0.029      0.042
# nativeAmerican_w1  0.005        0.048   -0.062   -0.005   -0.003             1.000        0.078         -0.022        -0.014      0.048
# raceOther_w1       0.009        0.382   -0.177   -0.101    0.030             0.078        1.000          0.004        -0.029     -0.008
# sportPartic_w1    -0.097       -0.011    0.013   -0.003   -0.021            -0.022        0.004          1.000         0.053     -0.094
# selfEst_w3_sc     -0.079       -0.027   -0.031    0.058   -0.029            -0.014       -0.029          0.053         1.000     -0.311
# depress_w4         0.125        0.046   -0.091    0.063    0.042             0.048       -0.008         -0.094        -0.311      1.000

corrplot::corrplot(cor(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", "asian_w1", "nativeAmerican_w1", 
                                "raceOther_w1", "sportPartic_w1",  "selfEst_w3_sc", "depress_w4")]),
                   method = c("ellipse"))

PerformanceAnalytics::chart.Correlation(data[, c("ethnicity_w1", "white_w1", "black_w1", 
                                                 "asian_w1", "nativeAmerican_w1", 
                                                 "raceOther_w1", "sportPartic_w1",  "selfEst_w3_sc", "depress_w4")])

# Using FE models to check black_w1 & raceOther_w1 possible relations to mediator & outcome 

## PS 
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                ethnicity_w1 + white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))

## Mediator 
med_fe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + black_w1 + raceOther_w1 + as.factor(CLUSTER2)", 
              data = data, 
              weights = iptw_fe)
## Outcome 
out_fe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + black_w1 + raceOther_w1 + as.factor(CLUSTER2)",
              data = data,
              weights = iptw_fe)

summary(med_fe)
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            18.65502    0.70251  26.555  < 2e-16 ***
#   sportPartic_w1          0.20465    0.10301   1.987   0.0471 *  
#   age_w1_sc              -0.01065    0.06394  -0.167   0.8677    
# sex_w1                 -0.49201    0.09340  -5.268 1.48e-07 ***
#   parentalEdu_w1_sc       0.05964    0.05123   1.164   0.2445    
# familyStruct_w1         0.39100    0.09676   4.041 5.46e-05 ***
#   healthInsur_w3          0.78904    0.13287   5.938 3.21e-09 ***
#   ethnicity_w1           -0.11074    0.18729  -0.591   0.5544    
# white_w1               -0.14313    0.18259  -0.784   0.4332    
# asian_w1               -0.49341    0.26422  -1.867   0.0619 .  
# nativeAmerican_w1      -0.12144    0.20794  -0.584   0.5592    
# black_w1               -0.03602    0.19924  -0.181   0.8566    
# raceOther_w1           -0.18510    0.21092  -0.878   0.3802    
# AIC: 15005

#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            17.80757    0.72026  24.724  < 2e-16 ***
#   sportPartic_w1          0.07551    0.09175   0.823  0.41055    
# age_w1_sc               0.04970    0.06267   0.793  0.42780    
# sex_w1                 -0.42352    0.07816  -5.419 6.48e-08 ***
#   parentalEdu_w1_sc       0.06534    0.05075   1.287  0.19804    
# familyStruct_w1         0.44539    0.09850   4.522 6.37e-06 ***
#   healthInsur_w3          0.83042    0.12978   6.399 1.81e-10 ***
#   ethnicity_w1           -0.25120    0.18717  -1.342  0.17968    
# white_w1               -0.10971    0.18064  -0.607  0.54366    
# asian_w1               -0.52110    0.26997  -1.930  0.05368 .  
# nativeAmerican_w1      -0.03361    0.20111  -0.167  0.86727    
# black_w1                0.11970    0.19980   0.599  0.54914    
# raceOther_w1           -0.10991    0.20542  -0.535  0.59263    
# AIC: 15262
summary(out_fe)
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             7.75468    1.15400   6.720 2.17e-11 ***
#   selfEst_w3_sc          -1.35032    0.07826 -17.254  < 2e-16 ***
#   sportPartic_w1         -0.32940    0.16865  -1.953  0.05089 .  
# age_w1_sc               0.02045    0.10452   0.196  0.84488    
# sex_w1                  0.93863    0.15404   6.094 1.24e-09 ***
#   parentalEdu_w1_sc      -0.24334    0.08371  -2.907  0.00368 ** 
#   familyStruct_w1        -0.67175    0.15798  -4.252 2.18e-05 ***
#   healthInsur_w4         -0.98189    0.22402  -4.383 1.21e-05 ***
#   ethnicity_w1            0.12796    0.30674   0.417  0.67658    
# white_w1               -0.53639    0.29832  -1.798  0.07227 .  
# asian_w1                0.21890    0.43249   0.506  0.61281    
# nativeAmerican_w1       0.70689    0.33919   2.084  0.03724 *  
#   black_w1                0.12222    0.32544   0.376  0.70729    
# raceOther_w1           -0.67769    0.34608  -1.958  0.05030 .  
# AIC: 18159

#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             5.77607    1.17845   4.901 1.00e-06 ***
#   selfEst_w3_sc          -1.38810    0.07825 -17.739  < 2e-16 ***
#   sportPartic_w1         -0.14384    0.15054  -0.956 0.339377    
# age_w1_sc               0.07599    0.10249   0.741 0.458472    
# sex_w1                  1.20881    0.12919   9.357  < 2e-16 ***
#   parentalEdu_w1_sc      -0.27807    0.08298  -3.351 0.000815 ***
#   familyStruct_w1        -0.69908    0.16059  -4.353 1.39e-05 ***
#   healthInsur_w4         -1.02305    0.21809  -4.691 2.84e-06 ***
#   ethnicity_w1            0.14202    0.30749   0.462 0.644203    
# white_w1               -0.51976    0.29555  -1.759 0.078748 .  
# asian_w1                0.34828    0.44229   0.787 0.431082    
# nativeAmerican_w1       0.86326    0.32722   2.638 0.008379 ** 
#   black_w1                0.29481    0.32659   0.903 0.366758    
# raceOther_w1           -0.97026    0.33726  -2.877 0.004045 ** 
#   AIC: 18419


# Perform bootstrap resampling 
med_black_boot <- numeric(1000)
med_other_boot <- numeric(1000)
out_black_boot <- numeric(1000)
out_other_boot <- numeric(1000)

set.seed(456)
for (i in 1:1000) {
  # Resample with replacement at the cluster level
  cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
  data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
  # Fit the FE models for the bootstrap sample
  mediator_fe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + black_w1 + raceOther_w1 + as.factor(CLUSTER2)", 
      data = data_boot, #)
      weights = iptw_fe)
  outcome_fe <-glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + black_w1 + raceOther_w1 + as.factor(CLUSTER2)",
      data = data_boot, #)
      weights = iptw_fe)
  # Calculate the effects for the bootstrap sample
  med_black_boot[i] <- summary(mediator_fe)$coef["black_w1", "Estimate"] 
  med_other_boot[i] <- summary(mediator_fe)$coef["raceOther_w1", "Estimate"] 
  out_black_boot[i] <- summary(outcome_fe)$coef["black_w1", "Estimate"] 
  out_other_boot[i] <- summary(outcome_fe)$coef["raceOther_w1", "Estimate"] 

  # Progress message 
  if (i %% 100 == 0) cat("Completed iteration", i, "\n")
}

# Calculate the percentile bootstrap CI
(boot_med_black <- quantile(med_black_boot, probs = c(0.025, 0.975)))
(boot_med_other <- quantile(med_other_boot, probs = c(0.025, 0.975)))
(boot_out_black <- quantile(out_black_boot, probs = c(0.025, 0.975)))
(boot_out_other <- quantile(out_other_boot, probs = c(0.025, 0.975)))

summary(med_fe)$coef["black_w1", "Estimate"]; (boot_med_black <- quantile(med_black_boot, probs = c(0.025, 0.975)))
summary(out_fe)$coef["black_w1", "Estimate"]; (boot_out_black <- quantile(out_black_boot, probs = c(0.025, 0.975)))

summary(med_fe)$coef["raceOther_w1", "Estimate"]; (boot_med_other <- quantile(med_other_boot, probs = c(0.025, 0.975)))
summary(out_fe)$coef["raceOther_w1", "Estimate"]; (boot_out_other <- quantile(out_other_boot, probs = c(0.025, 0.975))) 
# Only the other racial category appears to have a significant relation with the outcome, -0.678, 95% CI[-1.209, -0.171] (when IPTW is used: -0.970, 95% CI[-1.549, -0.422] )


## Group native american & other and drop black dummy variables ------------

# Group asian, native american, & other 
data$Other <- rowSums(data[, c("nativeAmerican_w1", "raceOther_w1")])

paste0("Number of white: ", prettyNum(sum(data$white_w1), big.mark = ","), " (", round((sum(data$white_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of black: ", prettyNum(sum(data$black_w1), big.mark = ","), " (", round((sum(data$black_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of asian: ", prettyNum(sum(data$asian_w1), big.mark = ","), " (", round((sum(data$asian_w1) / nrow(data)) * 100, 2), "%)")
paste0("Number of other: ", prettyNum(sum(data$Other), big.mark = ","), " (", round((sum(data$Other) / nrow(data)) * 100, 2), "%)")
# [1] "Number of white: 1,965 (62.36%)"
# [1] "Number of black: 804 (25.52%)"
# [1] "Number of asian: 147 (4.67%)"
# [1] "Number of other: 398 (12.63%)"
round(cor(data[, c("sex_w1",
                   "ethnicity_w1",
                   "white_w1",
                   "black_w1",
                   "asian_w1", 
                   "Other",
                   "sportPartic_w1")]),
      digits = 3)
#                 sex_w1 ethnicity_w1 white_w1 black_w1 asian_w1  Other sportPartic_w1
# sex_w1          1.000       -0.022   -0.025    0.049   -0.016  0.010         -0.097
# ethnicity_w1   -0.022        1.000   -0.322   -0.088    0.014  0.301         -0.011
# white_w1       -0.025       -0.322    1.000   -0.686   -0.220 -0.165          0.013
# black_w1        0.049       -0.088   -0.686    1.000   -0.040 -0.075         -0.003
# asian_w1       -0.016        0.014   -0.220   -0.040    1.000  0.019         -0.021
# Other           0.010        0.301   -0.165   -0.075    0.019  1.000         -0.012
# sportPartic_w1 -0.097       -0.011    0.013   -0.003   -0.021 -0.012          1.000

corrplot::corrplot(cor(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", 
                                "asian_w1", "Other", "sportPartic_w1")]),
                   method = c("ellipse"))

PerformanceAnalytics::chart.Correlation(data[, c("sex_w1", "ethnicity_w1", "white_w1", "black_w1", 
                                                     "Other", "sportPartic_w1")])

# Run RE PS model & add predictors (with new Other variable)
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + 
                        ethnicity_w1 + white_w1 + asian_w1 + Other +  
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) # Model converges when Other (other racial groups + native american) is added, but fails to converge when black_w1 is added 
# Note: white_w1 corr with Other (r = -0.17); and ethnicity_w1 corr with Other (r = 0.30)






# PS & IPTW (Dropped black_w1 & otherRace_w1) ------------------------------------

## SL PS model -------------------------------------------------------------
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
                parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1",
                family = "binomial", 
                data = data)
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))


## FE PS model -------------------------------------------------------------
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
                parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))


## RE PS model -------------------------------------------------------------
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
                        ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
                        parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) 
# control = lme4::glmerControl(optimizer = "bobyqa")) # Changing optimizer
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
    formula = sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
      ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1, 
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
    formula = sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
      ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + as.factor(CLUSTER2),
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
#     formula = sportPartic_w1 ~ feelings_w1_sc + age_w1_sc + sex_w1 + 
#       ethnicity_w1 + white_w1 + asian_w1 + nativeAmerican_w1 + 
#       parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + (1 | CLUSTER2),
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

# SL 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_sl < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_sl > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.01)), 
  "): ", sum(I(data$iptw_sl < quantile(data$iptw_sl, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.99)), 
  "): ", sum(I(data$iptw_sl > quantile(data$iptw_sl, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_sl, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < percentiles[1], percentiles[1], 
                          ifelse(iptw_sl > percentiles[2], percentiles[2], 
                                 iptw_sl)))

# FE 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_fe < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_fe > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.01)), 
  "): ", sum(I(data$iptw_fe < quantile(data$iptw_fe, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.99)), 
  "): ", sum(I(data$iptw_fe > quantile(data$iptw_fe, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_fe, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < percentiles[1], percentiles[1], 
                          ifelse(iptw_fe > percentiles[2], percentiles[2], 
                                 iptw_fe)))

# RE
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_re < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_re > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_re, probs = c(0.01)), 
  "): ", sum(I(data$iptw_re < quantile(data$iptw_re, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_re, probs = c(0.99)), 
  "): ", sum(I(data$iptw_re > quantile(data$iptw_re, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_re, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_re = ifelse(iptw_re < percentiles[1], percentiles[1], 
                          ifelse(iptw_re > percentiles[2], percentiles[2], 
                                 iptw_re)))






### SL ----------------------------------------------------------------------

##### fixed 

data$ps_sl
data %>% 
  summarize()
sum(I(data$ps_sl < 0.0111))
sum(I(data$ps_sl > 0.999))

quantile(data$iptw_sl, probs = c(0.01, 0.1, 0.90, 0.99))
sum(I(data$iptw_sl < 1.259809))
sum(I(data$iptw_sl > 4.439679))

hist(data$ps_sl)
hist(data$iptw_sl)

##### dynamic threshold -------------------------------------------------------
# Identify nonoverlapping cases with a caliper of 0.05
caliper <- 0.05
data_sl <- data
data_sl <- data_sl %>% 
  mutate(
    overlap_left = max(tapply(ps_sl_logit, sportPartic_w1, min)) - sd(ps_sl_logit) * caliper,
    overlap_right = min(tapply(ps_sl_logit, sportPartic_w1, max)) + sd(ps_sl_logit) * caliper, 
    # ind_nonoverlap = case_when(
    #   ps_sl_logit > overlap_right ~ TRUE,
    #   ps_sl_logit < overlap_left ~ TRUE,
    #   TRUE ~ FALSE
    # ), 
    ind_nonoverlap_left = case_when(
      ps_sl_logit < overlap_left ~ TRUE,
      TRUE ~ FALSE
    ), 
    ind_nonoverlap_right = case_when(
      ps_sl_logit > overlap_right ~TRUE, 
      TRUE ~ FALSE
    ), 
    ind_nonoverlap = case_when(
      (ind_nonoverlap_left == TRUE | ind_nonoverlap_right == TRUE) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

paste0(
  "Number of non-overlapping cases: ",
  sum(data_sl$ind_nonoverlap),
  " (left: ",
  sum(data_sl$ind_nonoverlap_left),
  "; right: ",
  sum(data_sl$ind_nonoverlap_right),
  ")"
)

# drop nonoverlap cases
data_sl <- data_sl[data_sl$ind_nonoverlap == FALSE, ]


### FE ----------------------------------------------------------------------
##### dynamic threshold -------------------------------------------------------
caliper <- 0.05
data_fe <- data
data_fe <- data_fe %>% 
  mutate(
    overlap_left = max(tapply(ps_fe_logit, sportPartic_w1, min)) - sd(ps_fe_logit) * caliper,
    overlap_right = min(tapply(ps_fe_logit, sportPartic_w1, max)) + sd(ps_fe_logit) * caliper, 
    # ind_nonoverlap = case_when(
    #   ps_fe_logit > overlap_right ~ TRUE,
    #   ps_fe_logit < overlap_left ~ TRUE,
    #   TRUE ~ FALSE
    # ), 
    ind_nonoverlap_left = case_when(
      ps_fe_logit < overlap_left ~ TRUE,
      TRUE ~ FALSE
    ), 
    ind_nonoverlap_right = case_when(
      ps_fe_logit > overlap_right ~TRUE, 
      TRUE ~ FALSE
    ), 
    ind_nonoverlap = case_when(
      (ind_nonoverlap_left == TRUE | ind_nonoverlap_right == TRUE) ~ TRUE, 
      TRUE ~ FALSE
    )
  )

paste0(
  "Number of non-overlapping cases: ",
  sum(data_fe$ind_nonoverlap),
  " (left: ",
  sum(data_fe$ind_nonoverlap_left),
  "; right: ",
  sum(data_fe$ind_nonoverlap_right),
  ")"
)

# drop nonoverlap cases 
data_fe <- data_fe[data_fe$ind_nonoverlap == FALSE, ]


### RE ----------------------------------------------------------------------
##### fixed threshold ---------------------------------------------------------
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


##### dynamic threshold -------------------------------------------------------
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

paste0(
  "Number of non-overlapping cases: ",
  sum(data_re$ind_nonoverlap),
  " (left: ",
  sum(data_re$ind_nonoverlap_left),
  "; right: ",
  sum(data_re$ind_nonoverlap_right),
  ")"
)

# drop nonoverlap cases
data_re <- data_re[data_re$ind_nonoverlap == FALSE, ]


# Check Overlap -----------------------------------------------------------

### SL ----------------------------------------------------------------------

sl_overlap <- list()

# Add stacked overlap visual to compare before & after dropping 
(sl_overlap$histogram <- 
    rbind(
      data.frame(ps_sl_logit = data$ps_sl_logit, 
                 sportPartic_w1 = as.factor(data$sportPartic_w1),
                 source = "pre"),
      data.frame(ps_sl_logit = data_sl$ps_sl_logit, 
                 sportPartic_w1 = as.factor(data_sl$sportPartic_w1),
                 source = "post")
    ) %>% 
    ggplot(aes(x = ps_sl_logit, 
               group = sportPartic_w1, 
               fill = sportPartic_w1)) +
      geom_histogram(position = "identity",
                     alpha = 0.5,
                     binwidth = 0.1) +
      geom_density(aes(y = ..count.. * 0.1, 
                       color = sportPartic_w1),  
                   size = 0.5,        
                   alpha = 0, 
                   trim = TRUE, 
                   show.legend = FALSE) +      # Set alpha to 1 (no transparency)
    facet_wrap(~ source, ncol = 1) +
    scale_fill_manual(values = c("blue", "darkorange")) +
    scale_color_manual(values = c("blue", "darkorange")) +
    theme_minimal() +
    labs(# title = "Histogram of ps_sl_logit by sportPartic_w1",
      y = "count",
      fill = "trt (Sport Participation)") +
    theme(legend.position = "bottom"))

# Add summary stats for PS logit before & after dropping cases 
sl_overlap$pre_summary <- data %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_sl_logit), 
            mean = mean(ps_sl_logit), 
            max = max(ps_sl_logit))

sl_overlap$post_summary <- data_sl %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_sl_logit), 
            mean = mean(ps_sl_logit), 
            max = max(ps_sl_logit))


# Display 
sl_overlap$histogram
sl_overlap$pre_summary
sl_overlap$post_summary

# > sl_overlap$pre_summary            # With 17 largest schools 
#   sportPartic_w1    min  mean   max
# 1              0 -1.17  0.370  2.32
# 2              1 -0.816 0.681  1.79
# > sl_overlap$post_summary
#   sportPartic_w1    min  mean   max
# 1              0 -0.837 0.383  1.76
# 2              1 -0.816 0.681  1.79

# > sl_overlap$pre_summary            # With 62 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.24 0.355  1.83
# 2              1 -2.64 0.560  1.56
# > sl_overlap$post_summary
#   sportPartic_w1    min  mean   max
# 1              0 -1.24  0.351  1.46
# 2              1 -0.764 0.562  1.56

# > sl_overlap$pre_summary            # With 121 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.14 0.361  1.54
# 2              1 -2.67 0.546  1.46
# > sl_overlap$post_summary
#   sportPartic_w1    min  mean   max
# 1              0 -1.14  0.359  1.46
# 2              1 -0.700 0.547  1.46


### FE ----------------------------------------------------------------------

fe_overlap <- list()

# Add stacked overlap visual to compare before & after dropping 
(fe_overlap$histogram <- 
    rbind(
      data.frame(ps_fe_logit = data$ps_fe_logit, 
                 sportPartic_w1 = as.factor(data$sportPartic_w1),
                 source = "pre"),
      data.frame(ps_fe_logit = data_fe$ps_fe_logit, 
                 sportPartic_w1 = as.factor(data_fe$sportPartic_w1),
                 source = "post")
    ) %>% 
    ggplot(aes(x = ps_fe_logit, 
               group = sportPartic_w1, 
               fill = sportPartic_w1)) +
    geom_histogram(position = "identity",
                   alpha = 0.7,
                   binwidth = 0.1) +
    geom_density(aes(y = ..count.. * 0.1, 
                     color = sportPartic_w1),  
                 size = 0.5,        
                 alpha = 0, 
                 trim = TRUE, 
                 show.legend = FALSE) +      # Set alpha to 1 (no transparency)
    facet_wrap(~ source, ncol = 1) +
    scale_fill_manual(values = c("blue", "darkorange")) +
    scale_color_manual(values = c("blue", "darkorange")) +
    theme_minimal() +
    labs(# title = "Histogram of ps_fe_logit by sportPartic_w1",
      y = "count",
      fill = "trt (Sport Participation)") +
    theme(legend.position = "bottom"))

# Add summary stats for PS logit before & after dropping cases 
fe_overlap$pre_summary <- data %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_fe_logit), 
            mean = mean(ps_fe_logit), 
            max = max(ps_fe_logit))

fe_overlap$post_summary <- data_fe %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_fe_logit), 
            mean = mean(ps_fe_logit), 
            max = max(ps_fe_logit))


# Display 
fe_overlap$histogram
fe_overlap$pre_summary
fe_overlap$post_summary

# > fe_overlap$pre_summary            # With 17 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -2.14 0.146  3.69
# 2              1 -1.27 1.02   4.47
# > fe_overlap$post_summary
#   sportPartic_w1   min  mean   max
# 1              0 -1.30 0.176  3.69
# 2              1 -1.27 0.965  3.71

# > fe_overlap$pre_summary            # With 62 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.90 0.118  3.70
# 2              1 -2.77 0.869  4.30
# > fe_overlap$post_summary
#   sportPartic_w1   min  mean   max
# 1              0 -1.90 0.118  3.70
# 2              1 -1.41 0.832  3.75

# > fe_overlap$pre_summary            # With 121 largest schools 
# sportPartic_w1   min  mean   max
# 1              0 -2.00 0.111  3.73
# 2              1 -3.14 0.928 15.2 
# > fe_overlap$post_summary
# sportPartic_w1   min  mean   max
# 1              0 -2.00 0.111  3.73
# 2              1 -1.82 0.827  3.79


### RE ----------------------------------------------------------------------

re_overlap <- list()

# Add stacked overlap visual to compare before & after dropping 
(re_overlap$histogram <- 
  rbind(
    data.frame(ps_re_logit = data$ps_re_logit, 
               sportPartic_w1 = as.factor(data$sportPartic_w1),
               source = "pre"),
    data.frame(ps_re_logit = data_re$ps_re_logit, 
               sportPartic_w1 = as.factor(data_re$sportPartic_w1),
               source = "post")
  ) %>% 
  ggplot(aes(x = ps_re_logit, 
             group = sportPartic_w1, 
             fill = sportPartic_w1)) +
  geom_histogram(position = "identity",
                 alpha = 0.7,
                 binwidth = 0.1) +
    geom_density(aes(y = ..count.. * 0.1, 
                     color = sportPartic_w1),  
                 size = 0.5,        
                 alpha = 0, 
                 trim = TRUE, 
                 show.legend = FALSE) +      # Set alpha to 1 (no transparency)
    facet_wrap(~ source, ncol = 1) +
    scale_fill_manual(values = c("blue", "darkorange")) +
    scale_color_manual(values = c("blue", "darkorange")) +
  theme_minimal() +
  labs(# title = "Histogram of ps_re_logit by sportPartic_w1",
    y = "count",
    fill = "trt (Sport Participation)") +
  theme(legend.position = "bottom"))

# Add summary stats for PS logit before & after dropping cases 
re_overlap$pre_summary <- data %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_re_logit), 
            mean = mean(ps_re_logit), 
            max = max(ps_re_logit))

re_overlap$post_summary <- data_re %>% 
  group_by(sportPartic_w1) %>% 
  summarize(min = min(ps_re_logit), 
            mean = mean(ps_re_logit), 
            max = max(ps_re_logit))


# Display 
re_overlap$histogram
re_overlap$pre_summary
re_overlap$post_summary

# > re_overlap$pre_summary            # With 17 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.78 0.186  2.83
# 2              1 -1.11 0.903  3.54
# > re_overlap$post_summary
#   sportPartic_w1   min  mean   max
# 1              0 -1.07 0.218  2.83
# 2              1 -1.11 0.843  2.85

# > re_overlap$pre_summary            # With 62 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.46 0.178  2.83
# 2              1 -2.72 0.746  3.02
# > re_overlap$post_summary
#   sportPartic_w1   min  mean   max
# 1              0 -1.46 0.178  2.83
# 2              1 -1.10 0.743  2.86

# > re_overlap$pre_summary            # With 121 largest schools 
#   sportPartic_w1   min  mean   max
# 1              0 -1.48 0.193  2.72
# 2              1 -2.91 0.715  2.96
# > re_overlap$post_summary
#   sportPartic_w1   min  mean   max
# 1              0 -1.48 0.193  2.72
# 2              1 -1.09 0.710  2.74












# Estimate Effects --------------------------------------------------------

## Mediator models ---------------------------------------------------------

### Mediator models (SL) ----------------------------------------------------

# Single-Level PS & SL med/outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect PS & SL med/outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_fe
  )
# Random-Effect PS & SL med/outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_re
  )


### Mediator models (FE) --------------------------------------------------

# Single-Level PS & FE med/outcome
med_slfe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
med_fefe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Random-Effect PS & FE med/outcome
med_refe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_re)


### Mediator models (RE) ----------------------------------------------------

### Add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Single-Level PS & RE med/outcome
med_slre <-
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
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
    data = data_sl,
    weights = iptw_sl
  )
# Fixed-Effect
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4",
    data = data_fe,
    weights = iptw_fe
  )
# Random-Effect
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4",
    data = data_re,
    weights = iptw_re
  )


### Outcome models (FE) -----------------------------------------------------

# Single-Level PS & FE med/outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data_sl, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data_fe, 
      weights = iptw_fe)
# Random-Effect PS & FE med/outcome
out_refe <-glm(formula = "depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)",
      data = data_re,
      weights = iptw_re
)


### Outcome models (RE) -----------------------------------------------------

# ### If not done so with mediator model, add column of just 1s for level-2 weight
# data <- cbind(data, L2weight = rep(1, nrow(data)))
# data_sl <- cbind(data_sl, L2weight = rep(1, nrow(data_sl)))
# data_fe <- cbind(data_fe, L2weight = rep(1, nrow(data_fe)))
# data_re <- cbind(data_re, L2weight = rep(1, nrow(data_re)))

# Single-Level PS & RE med/outcome
out_slre <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data_sl,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data_fe,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
      white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
      parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data_re,
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

results_DF
#   cond        NDE          NIE  # With 62 largest schools & dropped nonoverlapping cases (SL: 3, FE: 18, RE: 4) 
# 1 slsl -0.3440008  0.001201179
# 2 fesl -0.2160818 -0.009420466
# 3 resl -0.2486070 -0.016260947
# 4 slfe -0.2526180 -0.045080793
# 5 fefe -0.2335502 -0.018090723
# 6 refe -0.2283702 -0.031976047
# 7 slre -0.2902622 -0.021014822
# 8 fere -0.2302072 -0.016367209
# 9 rere -0.2372861 -0.023401019

#   cond        NDE         NIE   # With 121 largest schools & dropped nonoverlapping cases (SL: 3, FE: 27, RE: 7) 
# 1 slsl -0.2438572 -0.01447974
# 2 fesl -0.1110024 -0.02847153
# 3 resl -0.1706718 -0.02770936
# 4 slfe -0.1397979 -0.04151474
# 5 fefe -0.1224656 -0.03418641
# 6 refe -0.1258049 -0.03783174
# 7 slre -0.1837684 -0.02902276
# 8 fere -0.1216793 -0.03356498
# 9 rere -0.1471450 -0.03248475

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

#   cond        NDE         NIE # With 121 largest schools 
# 1 slsl -0.2368474 -0.02117536
# 2 fesl -0.1184377 -0.03255229
# 3 resl -0.1793095 -0.02837372
# 4 slfe -0.1407412 -0.04280035
# 5 fefe -0.1439630 -0.03339718
# 6 refe -0.1445497 -0.03370138
# 7 slre -0.1783744 -0.03331416
# 8 fere -0.1399102 -0.03423868
# 9 rere -0.1611487 -0.03028568






# Bootstrap CI ------------------------------------------------------------

###### Investigating convergence issues when bootstrapping with RERE ----------------

# only 50 iterations 

# Initialize vectors to store the bootstrap indirect effects and convergence statuses
indirect_effects_boot_rere <- numeric(50)
mediator_converged <- logical(50)
outcome_converged <- logical(50)

# Create the progress bar
pb <- txtProgressBar(min = 0, max = length(outcome_converged), style = 3, width = 50, char = "=")

# Perform bootstrap resampling
set.seed(456)
for (i in 1:50) {
  # Resample with replacement at the cluster level
  cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
  data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
  
  # Fit the RE models with RE ps for the bootstrap sample
  ### Add column of just 1s for level-2 weight
  data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
  
  mediator_rere <- tryCatch({
    WeMix::mix(formula = selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
                 white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                 parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2), 
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  outcome_rere <- tryCatch({
    WeMix::mix(formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
                 white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                 parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2), 
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  # Check convergence and calculate indirect effect if both models converged
  mediator_converged[i] <- !is.null(mediator_rere)
  outcome_converged[i] <- !is.null(outcome_rere)
  
  if (mediator_converged[i] && outcome_converged[i]) {
    indirect_effects_boot_rere[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
      summary(outcome_rere)$coef["selfEst_w3", "Estimate"]
  } else {
    indirect_effects_boot_rere[i] <- NA
  }
  
  # Update the progress bar
  setTxtProgressBar(pb, i)
}

# Close the progress bar & print completion message
close(pb)
cat("\nBootstrap resampling completed!\n")

# Calculate the percentile bootstrap CI
boot_ci_rere <- quantile(indirect_effects_boot_rere, probs = c(0.025, 0.975), na.rm = TRUE)

# Print results
print(boot_ci_rere)
paste0("Number of converged mediator models: ", sum(mediator_converged), " (", (sum(mediator_converged)/length(mediator_converged))*100, "%)")
paste0("Number of converged outcome models: ", sum(outcome_converged), " (", (sum(outcome_converged)/length(outcome_converged))*100, "%)")
paste0("Number of iterations with both models converged: ", sum(mediator_converged & outcome_converged), " (", (sum(mediator_converged & outcome_converged)/length(mediator_converged & outcome_converged))*100, "%)")
# 2.5%      97.5% 
# -0.1764424  0.1432817 
# [1] "Number of converged mediator models: 28 (56%)"
# [1] "Number of converged outcome models: 28 (56%)"
# [1] "Number of iterations with both models converged: 17 (34%)"   # With 15 largest schools  

# 2.5%       97.5% 
# -0.08573095  0.04281563 
# [1] "Number of converged mediator models: 40 (80%)"
# [1] "Number of converged outcome models: 50 (100%)"
# [1] "Number of iterations with both models converged: 40 (80%)"   # With 121 largest schools 

# Print all indirect effects
cat("All indirect effects:\n")
print(indirect_effects_boot_rere)



###### After dropping nonoverlap - Investigating convergence issues when bootstrapping with RERE ----------------

# only 50 iterations 

# Initialize vectors to store the bootstrap indirect effects and convergence statuses
indirect_effects_boot_rere <- numeric(50)
mediator_converged <- logical(50)
outcome_converged <- logical(50)

# Create the progress bar
pb <- txtProgressBar(min = 0, max = length(outcome_converged), style = 3, width = 50, char = "=")

# Perform bootstrap resampling
set.seed(456)
for (i in 1:50) {
  # Resample with replacement at the cluster level
  cluster_boot <- sample(unique(data_re$CLUSTER2), replace = TRUE)
  data_boot <- data_re[data_re$CLUSTER2 %in% cluster_boot, ]
  
  # Fit the RE models with RE ps for the bootstrap sample
  ### Add column of just 1s for level-2 weight
  data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
  
  mediator_rere <- tryCatch({
    WeMix::mix(formula = selfEst_w3 ~ sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
                 white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                 parentalEdu_w1 + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2), 
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  outcome_rere <- tryCatch({
    WeMix::mix(formula = depress_w4 ~ selfEst_w3 + sportPartic_w1 + age_w1 + sex_w1 + ethnicity_w1 +
                 white_w1 + black_w1 + asian_w1 + nativeAmerican_w1 + raceOther_w1 +
                 parentalEdu_w1 + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2), 
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  # Check convergence and calculate indirect effect if both models converged
  mediator_converged[i] <- !is.null(mediator_rere)
  outcome_converged[i] <- !is.null(outcome_rere)
  
  if (mediator_converged[i] && outcome_converged[i]) {
    indirect_effects_boot_rere[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
      summary(outcome_rere)$coef["selfEst_w3", "Estimate"]
  } else {
    indirect_effects_boot_rere[i] <- NA
  }
  
  # Update the progress bar
  setTxtProgressBar(pb, i)
}

# Close the progress bar & print completion message
close(pb)
cat("\nBootstrap resampling completed!\n")

# Calculate the percentile bootstrap CI
boot_ci_rere <- quantile(indirect_effects_boot_rere, probs = c(0.025, 0.975), na.rm = TRUE)

# Print results
print(boot_ci_rere)
paste0("Number of converged mediator models: ", sum(mediator_converged), " (", (sum(mediator_converged)/length(mediator_converged))*100, "%)")
paste0("Number of converged outcome models: ", sum(outcome_converged), " (", (sum(outcome_converged)/length(outcome_converged))*100, "%)")
paste0("Number of iterations with both models converged: ", sum(mediator_converged & outcome_converged), " (", (sum(mediator_converged & outcome_converged)/length(mediator_converged & outcome_converged))*100, "%)")

# 2.5%       97.5% 
# -0.29022769  0.02240406 
# [1] "Number of converged mediator models: 34 (68%)"
# [1] "Number of converged outcome models: 20 (40%)"
# [1] "Number of iterations with both models converged: 12 (24%)"   # With 17 largest schools & dropped nonoverlapping cases (SL: 5, FE: 15, RE: 20) 

# 2.5%       97.5% 
# -0.10019603  0.06098333 
# [1] "Number of converged mediator models: 43 (86%)"
# [1] "Number of converged outcome models: 50 (100%)"
# [1] "Number of iterations with both models converged: 43 (86%)"   # With 62 largest schools & dropped nonoverlapping cases (SL: 3, FE: 18, RE: 4) 

# 2.5%       97.5% 
# -0.07437584  0.02966024 
# [1] "Number of converged mediator models: 36 (72%)"
# [1] "Number of converged outcome models: 50 (100%)"
# [1] "Number of iterations with both models converged: 36 (72%)"   # With 121 largest schools & dropped nonoverlapping cases (SL: 3, FE: 27, RE: 7) 


# 2.5%      97.5% 
# -0.1764424  0.1432817 
# [1] "Number of converged mediator models: 28 (56%)"
# [1] "Number of converged outcome models: 28 (56%)"
# [1] "Number of iterations with both models converged: 17 (34%)"   # With 15 largest schools  

# 2.5%       97.5% 
# -0.08573095  0.04281563 
# [1] "Number of converged mediator models: 40 (80%)"
# [1] "Number of converged outcome models: 50 (100%)"
# [1] "Number of iterations with both models converged: 40 (80%)"   # With 121 largest schools 

# Print all indirect effects
cat("All indirect effects:\n")
print(indirect_effects_boot_rere)





