################################################################################
################### QP Empirical Application - Data Cleaning ###################
################################################################################

############################ Script Description ################################
#
# Author: Cameron
# 
# Date Created: 04/16/24
#
#
# Script Description:
# This script imports and cleans data from the Add Health study, focusing on 
# multiple waves (I, III, and IV) to prepare it for analysis. It includes 
# importing raw data, selecting and renaming relevant variables, recoding 
# demographic and scale variables, scoring self-esteem and depression scales, 
# handling missing data using multiple imputation (MICE), and scaling selected 
# variables. The script also merges data across waves, filters clusters, and 
# ultimately exports a cleaned dataset ("/Data/Cleaned/Empirical-Application-Data.rds").
#
#
# Last Updated: 08/15/2024 
#
#
# Notes:
# 
# 
################################################################################


# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  tidyverse, 
  readr, 
  ggplot2, 
  dplyr, 
  stringr, 
  mice 
)



# Wave I ------------------------------------------------------------------

## Import Data -------------------------------------------------------------
# Import student data from Wave I of the Add Health study
# The Add Health study is a longitudinal study of adolescents in the United States
w1 <- readr::read_tsv(file = "Application/Data/Raw/ICPSR_21600/DS0001/21600-0001-Data.tsv")

## Select Relevant Variables -----------------------------------------------
# Select relevant variables from Wave I student data
# This includes demographic information, feelings scale (CES-D items), and self-esteem scale
w1 <- w1 %>% 
  select("AID":"SCH_YR", "S1":"S7", "S12", "S18", "S11", "S17", "PC22",
         "H1FS1":"H1FS19", # Feelings Scale (contains mostly CES-D items for depression assessment)
         "H1NF12B", "S44A14", "S44A18":"S44A29", 
         "H1PF30", "H1PF32":"H1PF34", "H1PF36") # Self-esteem scale items

## Rename Demographic Variables --------------------------------------------
# Rename key demographic variables for clarity
data <- w1 %>% 
  rename(age = S1, 
         sex = S2, 
         white = S6A, 
         black = S6B) 

## Clean and Score Self-Esteem Scale ---------------------------------------
# Replace invalid responses (6=refused; 8=don't know) with NA for self-esteem items
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]] <- 
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]], 2,
        function(x) ifelse(x >= 6, NA, x))

# Reverse code all 5-point scale items for self-esteem
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]] <- 
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1PF\\d{2}")]], 2, function(x) 6 - x) 

# Calculate total self-esteem score by summing individual item scores
data$selfEst <- rowSums(data[, grep(pattern = "^H1PF\\d{2}", colnames(data))])


## Recode Parental Education -----------------------------------------------
# Recode parental education
data$momEdu <- ifelse(data$S12 %in% c(1, 2, 3, 10), 1, 
                      ifelse(data$S12 >= 4 & data$S12 <= 6, 2, 
                             ifelse(data$S12 %in% c(7, 8), 3, NA)))

data$dadEdu <- ifelse(data$S18 %in% c(1, 2, 3, 10), 1, 
                      ifelse(data$S18 >= 4 & data$S18 <= 6, 2, 
                             ifelse(data$S18 %in% c(7, 8), 3, NA)))

# Create a parental education variable (highest level of either parent)
data$parentalEdu <- apply(data[, c("momEdu", "dadEdu")], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)  # Assign NA if both momEdu and dadEdu are NA
  } else {
    return(max(x, na.rm = TRUE))  # Take the max value, ignoring NA
  }
})

## Recode Family Structure -------------------------------------------------
# Recode family structure
data$familyStruct <- data$S11 + data$S17
data$familyStruct <- ifelse(data$familyStruct == 10, NA, data$familyStruct) # Change "multiple response" to missing

## Recode Sex  -------------------------------------------------------------
# Recode Sex 
data <- data %>% 
  mutate(sex = ifelse(sex >= 3, NA, sex))

## Calculate Sports Participation ------------------------------------------
# Calculate sports participation (H1NF12B)
data$sport <- rowSums(data[, colnames(data)[stringr::str_starts(colnames(data), pattern = "^S44A")]])
data$sportPartic <- ifelse(data$sport >= 1, 1, 
                           ifelse(data$sport == 0, 0, NA))

## Score Feelings Scale ----------------------------------------------------
# Change values for Feelings Scale items to NA
data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]] <-
  apply(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]], 2, function(x) ifelse(x >= 6, NA, x)) # 6=refused; 8=don't know

# Reverse code items 4, 8, 11, & 15 on the Feelings Scale
data[, c("H1FS4", "H1FS8", "H1FS11", "H1FS15")] <-
  apply(data[, c("H1FS4", "H1FS8", "H1FS11", "H1FS15")], 2, function(x) 3 - x)

# Score Feelings Scale 
data$feelings <- rowSums(data[, colnames(data)[str_detect(colnames(data), pattern = "^H1FS")]])

## Import Cluster Data -----------------------------------------------------
# Import wave I data to obtain cluster indicators (CLUSTER2)
# Note: This dataset may also contain sampling weights, but our focus here is on clusters
w1.w <- readr::read_tsv(file = "Application/Data/Raw/ICPSR_21600/DS0004/21600-0004-Data.tsv")

# Add cluster indicator (CLUSTER2) to the main dataset
data <- merge(data, w1.w[, c("AID", "CLUSTER2")], by = "AID")

## Reorder Variables -------------------------------------------------------
# Reorder variables & add wave number to variable names
data <- data %>%
  select(AID, CLUSTER2, everything()) %>%
  select("AID", "CLUSTER2", "age", "sex",
         "white", "black", 
         "parentalEdu", "familyStruct", "sport", "sportPartic", "feelings", "selfEst")

colnames(data)[-c(1:2)] <- paste0(colnames(data)[-c(1:2)], "_w1")

## Clean Environment -------------------------------------------------------
# Remove all objects from the environment except for the 'data' object
rm(list = setdiff(ls(), "data"))



# Wave III ----------------------------------------------------------------
# Wave III (21600-0008-Codebook) pg. 395

## Import Data -------------------------------------------------
# Import student data
w3 <- readr::read_tsv(file = "Application/Data/Raw/ICPSR_21600/DS0008/21600-0008-Data.tsv")

# Drop variables in wave III student data
w3 <- w3 %>% 
  select("AID":"BIO_SEX3", 
         "H3SP7", "H3SP19":"H3SP22") # self-esteem scale [NOTE: STILL MISSING AN ITEM]

## Self-Esteem Scale ------------------------------------------------

## “You have a lot of good qualities,” (H3SP19)
## “you have a lot to be proud of,” (H3SP20)
## “you like yourself just the way you are,” (H3SP21)
# used "Do you agree or disagree that you feel you are doing things just about right?" (H3SP22) 
### with responses on a 5-point scale ranging from 1 (strongly agree) to 5 (strongly disagree) (reverse coded)
## and “How often was the following been true over the past week? You felt that you were just as good as other people,” (H3SP7)
### with responses ranging from 0 (never or rarely) to 3 (most of the time or all of the time). Higher scores indicate greater self-esteem.28

# Change values for self-esteem items to NA 
w3$H3SP7 <- ifelse(w3$H3SP7 >= 6, NA, w3$H3SP7) # 6=refused; 8=don't know; 9=not applicable 

w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]] <-
  apply(w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]], 2,
        function(x) ifelse(x >= 96, NA, x)) # 96=refused; 98=don't know; 99=not applicable 

# Reverse code the four 5-point scale items
w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]] <-
  apply(w3[, colnames(w3)[str_detect(colnames(w3), pattern = "^H3SP\\d{2}")]], 2, function(x) 6 - x) 

# Score self-esteem scale 
w3$selfEst <- rowSums(w3[, grep(pattern = "^H3SP", colnames(w3))])


## Merge Wave III Data with Main Data Set -------------------------------
# Merge wave III student data to data set
data <- merge(data, w3[, c("AID", "selfEst")], by = "AID")

## Update Variable Names with Wave Number --------------------------------
# Add wave number to variable names
colnames(data)[!grepl(pattern = paste(c("_w1$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))] <-
  paste0(colnames(data)[!grepl(pattern = paste(c("_w1$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))], "_w3")



# Wave IV -----------------------------------------------------------------
# Wave IV (21600-0022-Codebook) pg. 184 

## Import Data -------------------------------------------------
# Import student data
w4 <- readr::read_tsv(file = "Application/Data/Raw/ICPSR_21600/DS0022/21600-0022-Data.tsv")

# Drop variables in Wave IV student data
w4 <- w4 %>% 
  select("AID":"BIO_SEX4", 
         "H4MH18":"H4MH27") # CES-D-10 scale 

## CES-D-10 Scale ----------------------------------------------------

# 1. You were bothered by things that usually don't bother you. (H4MH18)
# 2. (During the past seven days:) You could not shake off the blues, even with help from your family and your friends. (H4MH19)
# 3. (During the past seven days:) You felt you were just as good as other people. (H4MH20)
# 4. (During the past seven days:) You had trouble keeping your mind on what you were doing. (H4MH21)
# 5. (During the past seven days:) You felt depressed. (H4MH22)
# 6. (During the past seven days:) You felt that you were too tired to do things. (H4MH23)
# 7. (During the past seven days:) You felt happy. (H4MH24)
# 8. (During the past seven days:) You enjoyed life. (H4MH25)
# 9. (During the past seven days:) You felt sad. (H4MH26)
# 10. (During the past seven days:) You felt that people disliked you, during the past seven days. (H4MH27)

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

## Merge Data --------------------------------------------------------
# Add Wave IV depression score to the main dataset
data <- merge(data, w4[, c("AID", "depress")], by = "AID")

# Add wave number to variable names
colnames(data)[!grepl(pattern = paste(c("_w\\d{1}$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))] <- 
  paste0(colnames(data)[!grepl(pattern = paste(c("_w\\d{1}$", "AID", "CLUSTER2"), collapse = "|"), colnames(data))], "_w4")

## Clean Environment -------------------------------------------------------
# Remove all objects from the environment except for the 'data' object
rm(list = setdiff(ls(), "data"))



# Drop small clusters -----------------------------------------------------

## Calculate and Filter Cluster Sizes --------------------------------
# Add cluster (school) size variable
data <- data %>% 
  group_by(CLUSTER2) %>% 
  mutate(n = n()) %>% 
  ungroup()

# Drop observations with missing sport participation (treatment)
data <- data %>% 
  filter(!is.na(sportPartic_w1)) # 4208 - 3155 = 1,053 dropped 

# Recalculate cluster sizes after filtering
data <- data %>% 
  group_by(CLUSTER2) %>% 
  mutate(n = n()) %>% 
  ungroup()

# Drop clusters with fewer than 5 students
data <- data %>% 
  filter(n >= 5) # 4 dropped 

## Examine Final Cluster Sizes ---------------------------------------
# Check final cluster sizes in order of size
data %>% 
  group_by(CLUSTER2) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n = 132)  # 121 schools (cluster sizes ranged from 5 to 81 students)



# Missing Data ------------------------------------------------------------

## Data Summary and Scaling -------------------------------------------
# Check proportions of categorical variables & summary of scale variables
summary(data[, c(
  "sex_w1", "white_w1", "black_w1", 
  "parentalEdu_w1", "familyStruct_w1", 
  "sport_w1", "sportPartic_w1", "feelings_w1", "depress_w4"
)])   # Consider collapsing "asian", "nativeAmerican", and "raceOther" due to low proportions

# Scale selected variables, rename, & drop original non-scaled variables
data <- data %>%
  mutate(
    age_w1_sc = as.vector(scale(age_w1)),
    parentalEdu_w1_sc = as.vector(scale(parentalEdu_w1)),
    feelings_w1_sc = as.vector(scale(feelings_w1)), 
    selfEst_w1_sc = as.vector(scale(selfEst_w1))
  ) %>%
  select(-age_w1, -parentalEdu_w1, -feelings_w1)  # Remove original columns if desired

## Missing Data Handling ----------------------------------------------
# Select variables needed for propensity score (PS) model
t <- data %>% 
  select(c("sex_w1", "white_w1", "black_w1", 
           "sportPartic_w1", "selfEst_w3", "depress_w4", 
           "familyStruct_w1", "parentalEdu_w1_sc", 
           "age_w1_sc", "feelings_w1_sc", "selfEst_w1_sc"))

# Examine missing data pattern
md.pattern(t, rotate.names = TRUE)

# Impute missing data using MICE
imp <- mice(t, m = 1, seed = 87542) 

# Compare original and imputed data
# Example with stripplot
stripplot(imp, pch = 20, cex = 1.2)
# Example with xyplot
xyplot(imp, selfEst_w3 ~ depress_w4 | .imp, pch = 20, cex = 1.2)

# Replace missing values in the original dataset
t_imp <- complete(imp, "long", include = TRUE)
var <- colnames(t)
for (i in 1:length(var)) {
  data[, var[i]] <- t_imp[t_imp$.imp == 1, var[i]]
}

# Add a scaled version of the mediator for use in the outcome model 
data <- data %>%
  mutate(selfEst_w3_sc = as.vector(scale(selfEst_w3)))

## Final Data Check --------------------------------------------------
# Check proportions of categorical variables & summary after imputing values
summary(data[, c(
  "sex_w1", "white_w1", "black_w1", 
  "parentalEdu_w1_sc", "familyStruct_w1", 
  "sport_w1", "sportPartic_w1", "feelings_w1_sc", 
  "selfEst_w3_sc", "depress_w4"
)])

## Clean Environment -------------------------------------------------------
# Remove all objects from the environment except for the 'data' object
rm(list = setdiff(ls(), "data"))



# Export Clean Data -------------------------------------------------------
# Save cleaned data 
saveRDS(data, file = "Application/Data/Cleaned/Empirical-Application-Data.rds")


