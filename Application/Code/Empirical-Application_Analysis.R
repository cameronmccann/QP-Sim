################################################################################
##################### QP Empirical Application - Analysis ######################
################################################################################

############################ Script Description ################################
#
# Author: 
# 
# Date Created: 2024-04-16
#
#
# Script Description:
#   This R script performs the mediation analysis for the empirical application. 
#   First, it calculate the Intraclass Correlation Coefficients (ICC) for the 
#   mediator and outcome variables. Next, it computes propensity scores (PSs) and 
#   Inverse Probability of Treatment Weight (IPTW) using the Single-Level, 
#   Fixed-Effect, and Random-Effect PS models, with percentile bootstrap confidence 
#   intervals. Covariate balance is visualized. Finally, it visualizes the 
#   estimated effects to facilitate interpretation of the results.
# 
# Last Updated: 2025-03-19 
#
#
# Notes:
# 
################################################################################


# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  tidyverse, 
  ggplot2, 
  extrafont, 
  stringr, 
  mice, 
  cobalt, 
  WeightIt, 
  boot, 
  utils, 
  lme4, 
  WeMix, 
  parallel
)

# Load Functions 
# source("Application/Functions/monteCarloCI.R")
source("Application/Functions/monteCarloCIb.R")



# Import Data -----------------------------------------------------
# Load clean dataset 
data <- read_rds(file = "Application/Data/Cleaned/Empirical-Application-Data.rds")



# Intraclass Correlation Coefficients (ICC) --------------------------

## ICC for Mediator ---------------------------------------------------
# Calculate the ICC for the mediator (self-esteem) to assess the proportion of variance
# that is attributable to differences between clusters (schools).
med_unconditional <- lme4::lmer(selfEst_w3 ~ (1 | CLUSTER2), data = data) # Unconditional model
# summary(med_unconditional)

med_var <- data.frame(lme4::VarCorr(med_unconditional))[1, 4] # Level-2 (between-cluster) variance
med_res <- data.frame(lme4::VarCorr(med_unconditional))[2, 4] # Residual (within-cluster) variance
med_icc <- med_var / (med_var + med_res) # Calculate ICC for mediator
med_icc # Display the ICC (with 121 schools, ICC ~ 0.006)

## ICC for Outcome ----------------------------------------------------
# Calculate the ICC for the outcome (depression) to assess the proportion of variance
# that is attributable to differences between clusters (schools).
out_unconditional <- lme4::lmer(depress_w4 ~ (1 | CLUSTER2), data = data) # Unconditional model
# summary(out_unconditional)

out_var <- data.frame(lme4::VarCorr(out_unconditional))[1, 4] # Level-2 (between-cluster) variance
out_res <- data.frame(lme4::VarCorr(out_unconditional))[2, 4] # Residual (within-cluster) variance
out_icc <- out_var / (out_var + out_res) # Calculate ICC for outcome
out_icc # Display the ICC (with 121 schools, ICC ~ 0.018)

## Optional: ICC for Variables in Propensity Score (PS) Models -------------
# Uncomment the following lines to calculate ICCs for other variables used in PS models.

# ## Parent Education
# parent_uncond <- lme4::lmer(parentalEdu_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(parent_uncond))[1, 4] / 
#   (data.frame(lme4::VarCorr(parent_uncond))[1, 4] + data.frame(lme4::VarCorr(parent_uncond))[2, 4])

# ## Feelings Scale
# feeling_uncond <- lme4::lmer(feelings_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(feeling_uncond))[1, 4] / 
#   (data.frame(lme4::VarCorr(feeling_uncond))[1, 4] + data.frame(lme4::VarCorr(feeling_uncond))[2, 4])

# ## Age
# age_uncond <- lme4::lmer(age_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(age_uncond))[1, 4] / 
#   (data.frame(lme4::VarCorr(age_uncond))[1, 4] + data.frame(lme4::VarCorr(age_uncond))[2, 4])

# ## Self-Esteem (Wave 1)
# selfEst_uncond <- lme4::lmer(selfEst_w1_sc ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(selfEst_uncond))[1, 4] / 
#   (data.frame(lme4::VarCorr(selfEst_uncond))[1, 4] + data.frame(lme4::VarCorr(selfEst_uncond))[2, 4])

# ## White Ethnicity
# white_uncond <- lme4::lmer(white_w1 ~ (1 | CLUSTER2), data = data)
# data.frame(lme4::VarCorr(white_uncond))[1, 4] / 
#   (data.frame(lme4::VarCorr(white_uncond))[1, 4] + data.frame(lme4::VarCorr(white_uncond))[2, 4])



# PS & IPTW Calculation ---------------------------------------------------------------
# This script calculates propensity scores and Inverse Probability of Treatment Weights (IPTW)
# for three models: SL (Single-Level), FE (Fixed-Effect), and RE (Random-Effect).

### SL Propensity Score Model -----------------------------------------------
# Single-level logistic regression model for propensity score estimation
psmod_sl <- glm(
  formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + selfEst_w1_sc",
  family = "binomial", 
  data = data
)

# Predict propensity scores and log-odds for the SL model
data$ps_sl <- predict(psmod_sl, type = "response") # Propensity scores
data$ps_sl_logit <- predict(psmod_sl, type = "link") # Log-odds

# Calculate IPTW for the SL model
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))

### FE Propensity Score Model ----------------------------------------
# Fixed-effect logistic regression model for propensity score estimation.
psmod_fe <- glm(
  formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + selfEst_w1_sc +
                as.factor(CLUSTER2)",
  family = "binomial", 
  data = data
)

# Predict propensity scores and log-odds for the FE model
data$ps_fe <- predict(psmod_fe, type = "response") # Propensity scores
data$ps_fe_logit <- predict(psmod_fe, type = "link") # Log-odds

# Calculate IPTW for the FE model
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))

### RE Propensity Score Model ----------------------------------------
# Random-effect logistic regression model for propensity score estimation.
psmod_re <- lme4::glmer(
  formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + selfEst_w1_sc + 
                        (1 | CLUSTER2)",
  family = "binomial", 
  data = data
)

# Predict propensity scores and log-odds for the RE model
data$ps_re <- predict(psmod_re, type = "response") # Propensity scores
data$ps_re_logit <- predict(psmod_re, type = "link") # Log-odds

# Calculate IPTW for the RE model
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))

### RE-Mean Propensity Score Model -----------------------------------------------

# # Calculate the cluster mean of covariates 
# data <- data |> 
#   group_by(CLUSTER2) |> 
#   mutate(
#     cluster_mean_feelings_w1_sc = mean(feelings_w1_sc, na.rm = TRUE), 
#     cluster_mean_sex_w1 = mean(sex_w1, na.rm = TRUE), 
#     cluster_mean_age_w1_sc = mean(age_w1_sc, na.rm = TRUE),
#     cluster_mean_white_w1 = mean(white_w1, na.rm = TRUE),
#     cluster_mean_black_w1 = mean(black_w1, na.rm = TRUE),
#     cluster_mean_parentalEdu_w1_sc = mean(parentalEdu_w1_sc, na.rm = TRUE),
#     cluster_mean_familyStruct_w1 = mean(familyStruct_w1, na.rm = TRUE),
#     cluster_mean_selfEst_w1_sc = mean(selfEst_w1_sc, na.rm = TRUE)
#   )
# 
# # RE with cluster means logistic regression model for propensity score estimation
# psmod_remean <- lme4::glmer(
#   formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
#   white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + selfEst_w1_sc +
#   cluster_mean_feelings_w1_sc + cluster_mean_sex_w1 + cluster_mean_age_w1_sc + 
#   cluster_mean_white_w1 + cluster_mean_black_w1 + cluster_mean_parentalEdu_w1_sc + 
#   cluster_mean_familyStruct_w1 + cluster_mean_selfEst_w1_sc + (1 | CLUSTER2)", 
#   family = "binomial", 
#   data = data
# )
# 
# # Predict propensity scores and log-odds for the RE-Mean model
# data$ps_remean <- predict(psmod_remean, type = "response") # Propensity scores
# data$ps_remean_logit <- predict(psmod_remean, type = "link") # Log-odds
# 
# # Calculate IPTW for the RE model
# data <- cbind(data, iptw_remean = with(data, (sportPartic_w1 / ps_remean) + (1 - sportPartic_w1) / (1 - ps_remean)))







# Truncate Non-overlap Cases ---------------------------------------------------------

### SL (Single-Level) ------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_sl < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_sl > 0.99)))
# [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"

# Determine the 1st and 99th percentiles of the IPTW distribution for the SL model.
# Cases below the 1st percentile or above the 99th percentile are considered outliers.
first_percentile <- quantile(data$iptw_sl, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_sl, probs = 0.99)

# Count the number of cases below the 1st percentile and above the 99th percentile of IPTW.
paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_sl < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_sl > ninety_ninth_percentile))
)
# [1] "Number of cases < 1st percentile of IPTW (1.25675888749596): 32; Number of cases > 99th percentile of IPTW (4.31652336904921): 32"

# Adjust IPTW values to the 1st and 99th percentile thresholds
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < first_percentile, first_percentile, 
                          ifelse(iptw_sl > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_sl)))

### FE (Fixed-Effect) -------------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_fe < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_fe > 0.99)))
# [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 11"

# Determine the 1st and 99th percentiles of the IPTW distribution for the FE model.
first_percentile <- quantile(data$iptw_fe, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_fe, probs = 0.99)

# Count the number of cases below the 1st percentile and above the 99th percentile of IPTW.
paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_fe < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_fe > ninety_ninth_percentile))
)
# [1] "Number of cases < 1st percentile of IPTW (1.02483454975511): 32; Number of cases > 99th percentile of IPTW (6.61712855938224): 32"

# Adjust IPTW values to the 1st and 99th percentile thresholds 
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < first_percentile, first_percentile, 
                          ifelse(iptw_fe > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_fe)))

### RE (Random-Effect) -------------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_re < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_re > 0.99)))
# [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"

# Determine the 1st and 99th percentiles of the IPTW distribution for the RE model.
first_percentile <- quantile(data$iptw_re, probs = 0.01)
ninety_ninth_percentile <- quantile(data$iptw_re, probs = 0.99)

# Count the number of cases below the 1st percentile and above the 99th percentile of IPTW.
paste0(
  "Number of cases < 1st percentile of IPTW (", first_percentile, 
  "): ", sum(I(data$iptw_re < first_percentile)), 
  "; Number of cases > 99th percentile of IPTW (", ninety_ninth_percentile, 
  "): ", sum(I(data$iptw_re > ninety_ninth_percentile))
)
# [1] "Number of cases < 1st percentile of IPTW (1.10193371886911): 32; Number of cases > 99th percentile of IPTW (4.95620683646533): 32"

# Adjust IPTW values to the 1st and 99th percentile thresholds 
data <- data %>% 
  mutate(iptw_re = ifelse(iptw_re < first_percentile, first_percentile, 
                          ifelse(iptw_re > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_re)))



# Covariate Balance Visualization -----------------------------------

## Create Functions --------------------------------------------------------
# Function to calculate weighted variance, which accounts for different sample weights
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    # Remove NA values from the variables and weights
    valid <- !is.na(x) & !is.na(w)
    x <- x[valid]
    w <- w[valid]
  }
  
  # Calculate the weighted mean
  wm <- weighted.mean(x, w)
  
  # Calculate the weighted variance using the weighted mean
  variance <- sum(w * (x - wm)^2) / sum(w)
  
  return(variance)  # Return the calculated weighted variance
}

# Function to calculate the Standardized Mean Difference (SMD) between two groups
calculate_smd <- function(data, treatment, covariate) {
  # Calculate the mean of the covariate for the treatment group
  mean_treatment <- mean(data[[covariate]][data[[treatment]] == 1], na.rm = TRUE)
  
  # Calculate the mean of the covariate for the control group
  mean_control <- mean(data[[covariate]][data[[treatment]] == 0], na.rm = TRUE)
  
  # Calculate the standard deviation for the treatment group
  sd_treatment <- sd(data[[covariate]][data[[treatment]] == 1], na.rm = TRUE)
  
  # Calculate the standard deviation for the control group
  sd_control <- sd(data[[covariate]][data[[treatment]] == 0], na.rm = TRUE)
  
  # Compute the SMD using the pooled standard deviations
  smd <- (mean_treatment - mean_control) / sqrt((sd_treatment^2 + sd_control^2) / 2)
  
  return(smd)  # Return the calculated SMD
}

# Function to calculate weighted SMD, which includes the use of weights for covariate balance
calculate_weighted_smd <- function(data, treatment, covariate, weights_col) {
  # Subset the data into treatment and control groups
  treatment_group <- data[data[[treatment]] == 1, ]
  control_group <- data[data[[treatment]] == 0, ]
  
  # Return NA if either group is empty after subsetting
  if (nrow(treatment_group) == 0 || nrow(control_group) == 0) {
    return(NA)
  }
  
  # Remove rows with missing values in the covariate or weights
  treatment_group <- treatment_group[!is.na(treatment_group[[covariate]]) & !is.na(treatment_group[[weights_col]]), ]
  control_group <- control_group[!is.na(control_group[[covariate]]) & !is.na(control_group[[weights_col]]), ]
  
  # Return NA if either group is empty after removing NAs
  if (nrow(treatment_group) == 0 || nrow(control_group) == 0) {
    return(NA)
  }
  
  # Calculate weighted means for both treatment and control groups
  mean_treatment <- weighted.mean(treatment_group[[covariate]], treatment_group[[weights_col]], na.rm = TRUE)
  mean_control <- weighted.mean(control_group[[covariate]], control_group[[weights_col]], na.rm = TRUE)
  
  # Calculate weighted standard deviations for both groups
  sd_treatment <- sqrt(weighted.var(treatment_group[[covariate]], treatment_group[[weights_col]], na.rm = TRUE))
  sd_control <- sqrt(weighted.var(control_group[[covariate]], control_group[[weights_col]], na.rm = TRUE))
  
  # Compute the weighted SMD using the weighted means and pooled standard deviations
  smd <- (mean_treatment - mean_control) / sqrt((sd_treatment^2 + sd_control^2) / 2)
  
  return(smd)  # Return the calculated weighted SMD
}

## Covariate Balance Calculations -----------------------------------

# Define the covariates to be assessed for balance
covariates <- c("feelings_w1_sc", "sex_w1", "age_w1_sc", "white_w1", 
                "black_w1", "parentalEdu_w1_sc", "familyStruct_w1", 
                "selfEst_w1_sc")

# Calculate the SMD for each covariate before applying weights (unweighted)
smd_before <- sapply(covariates, function(cov) calculate_smd(data, "sportPartic_w1", cov))

# Convert the SMD results into a data frame for easier handling and labeling
smd_before_df <- data.frame(covariate = covariates, SMD = smd_before)
smd_before_df$type <- "Unweighted"  # Label the SMDs as 'Unweighted'

# Calculate weighted SMD for each covariate using different weighting methods
smd_sl_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_sl"))
smd_fe_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_fe"))
smd_re_after <- sapply(covariates, function(cov) calculate_weighted_smd(data, "sportPartic_w1", cov, "iptw_re"))

# Convert the weighted SMD results into data frames
smd_sl_after_df <- data.frame(covariate = covariates, SMD = smd_sl_after, type = "Single-Level")
smd_fe_after_df <- data.frame(covariate = covariates, SMD = smd_fe_after, type = "Fixed-Effect")
smd_re_after_df <- data.frame(covariate = covariates, SMD = smd_re_after, type = "Random-Effect")

# Combine all SMD data frames into one for comprehensive analysis
smd_combined <- rbind(smd_before_df, smd_sl_after_df, smd_fe_after_df, smd_re_after_df)

# Calculate the Absolute Standardized Mean Difference (ASMD) for all covariates
smd_combined$ASMD <- abs(smd_combined$SMD)

## Visualization: Love Plot ------------------------------------------

# Define a custom order for the covariates on the y-axis in the plot
custom_order <- c("black_w1", "white_w1", "familyStruct_w1", 
                  "age_w1_sc", "sex_w1", 
                  "selfEst_w1_sc", "feelings_w1_sc", "parentalEdu_w1_sc")

# Define new labels for the y-axis to improve readability
new_labels <- c("Race: Black", "Race: White", "Family Structure", 
                "Age", "Sex", 
                "Self-Esteem Score", "Feelings Scale Score", "Parental Education")

# Load Times New Roman font for publication-quality visuals
# loadfonts()

# Save the Love Plot visualization for the paper as a PDF
pdf("Application/Output/Visuals/Covariate-Balance_QP-Doc.pdf")

# Create and customize the Love Plot for the paper, using Times New Roman font
ggplot(smd_combined, aes(x = ASMD, y = factor(covariate, levels = custom_order), color = type, shape = type)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +  # Reference line for SMD threshold (0.1)
  geom_vline(xintercept = 0, color = "black") +  # Line at zero to indicate no difference
  geom_point(size = 3, stroke = 1.5) +  # Plot points with increased size and stroke 
  labs(#title = "Love Plot",  
    #subtitle = "Covariate Balance of Individual-Level Covariates",  
    x = "\n Absolute Standardized Mean Difference (ASMD)",
    y = "") +
  theme_minimal() +  
  theme(text = element_text(family = "Times New Roman"),
        axis.text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 10),  # Adjust y-axis text 
        # axis.title = element_text(size = 14),  # Increase axis title size for better visibility
        # plot.title = element_text(size = 16, face = "bold"),  # Bold title for emphasis
        # plot.subtitle = element_text(size = 14),  # Increase subtitle size for readability
        legend.position = "top") +  # Position legend at the top for easy reference
  scale_color_manual(values = c("Unweighted" = "#00A9B7",  # Teal for unweighted
                                "Single-Level" = "#333F48",  # Gray for single-level weighting
                                "Fixed-Effect" = "#BF5700",  # Orange for fixed-effect weighting
                                "Random-Effect" = "#A6CD57"),  # Green for random-effect weighting
                     name = NULL) +
  # scale_color_manual(values = c("Unweighted" = "#1f77b4",  # Blue for unweighted
  #                               "Single-Level" = "#2ca02c",  # Green for single-level weighting
  #                               "Fixed-Effect" = "#ff7f0e",  # Orange for fixed-effect weighting
  #                               "Random-Effect" = "#9467bd"),  # Purple for random-effect weighting
  #                    name = NULL) +
  scale_shape_manual(values = c("Unweighted" = 16,  # Circle shape for unweighted
                                "Single-Level" = 17,  # Triangle shape for single-level
                                "Fixed-Effect" = 15,  # Square shape for fixed-effect
                                "Random-Effect" = 18),  # Diamond shape for random-effect
                     name = NULL) +
  scale_y_discrete(labels = new_labels)  # Apply new labels to the y-axis

# Close the PDF device
dev.off()

# Save the Love Plot visualization for the paper as a PNG
ggsave(filename = "Application/Output/Visuals/Covariate-Balance_QP-Doc.pdf", 
       plot = last_plot(), 
       width = 6, 
       height = 7, 
       units = "in", 
       dpi = 300)

ggsave(filename = "Application/Output/Visuals/Covariate-Balance_QP-Doc.png", plot = last_plot())



# Estimate Effects --------------------------------------------------------

## Mediator models ---------------------------------------------------

### Single-Level (SL) ------------------------------------------
# Single-level mediator model with PS weights from single-level PS model.
med_slsl <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_sl
)

# Single-level mediator model with PS weights from fixed-effect PS model.
med_fesl <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_fe
)

# Single-level mediator model with PS weights from random-effect PS model.
med_resl <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_re
)

### Fixed-Effect (FE) ------------------------------------------
# Fixed-effect mediator model with PS weights from single-level PS model.
med_slfe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_sl
)

# Fixed-effect mediator model with PS weights from fixed-effect PS model.
med_fefe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_fe
)

# Fixed-effect mediator model with PS weights from random-effect PS model.
med_refe <- glm(
  formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_re
)

### Random-Effect (RE) -----------------------------------------
# Add a column of ones for level-2 weights, required for the WeMix package.
data <- cbind(data, L2weight = rep(1, nrow(data)))

# Random-effect mediator model with PS weights from single-level PS model.
med_slre <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effect mediator model with PS weights from fixed-effect PS model.
med_fere <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effect mediator model with PS weights from random-effect PS model.
med_rere <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)

### Random-Effect with Cluster Means (RE-Mean) ---------------------------------
# Calculate the cluster mean of the treatment (sportPartic_w1) and mediator (selfEst_w3_sc)
data <- data %>%
  group_by(CLUSTER2) %>%
  mutate(
    cluster_mean_sportPartic_w1 = mean(sportPartic_w1, na.rm = TRUE),
    cluster_mean_selfEst_w3_sc = mean(selfEst_w3_sc, na.rm = TRUE)
  ) %>%
  ungroup()

# Random-effect mediator model with cluster means and PS weights from single-level PS model.
med_slre_cm <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + cluster_mean_sportPartic_w1 + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effect mediator model with cluster means and PS weights from fixed-effect PS model.
med_fere_cm <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + cluster_mean_sportPartic_w1 + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effect mediator model with cluster means and PS weights from random-effect PS model.
med_rere_cm <- WeMix::mix(
  formula = selfEst_w3 ~ sportPartic_w1 + cluster_mean_sportPartic_w1 + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)


## Outcome Models (TNDE & PNIE) --------------------------------------

### Single-Level (SL) ------------------------------------------
# Single-level outcome model with PS weights from single-level PS model.
out_slsl <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_sl
)

# Single-level outcome model with PS weights from fixed-effect PS model.
out_fesl <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_fe
)

# Single-level outcome model with PS weights from random-effect PS model.
out_resl <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_re
)

### Fixed-Effect (FE) ------------------------------------------
# Fixed-effect outcome model with PS weights from single-level PS model.
out_slfe <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_sl
)

# Fixed-effect outcome model with PS weights from fixed-effect PS model.
out_fefe <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_fe
)

# Fixed-effect outcome model with PS weights from random-effect PS model.
out_refe <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data,
  weights = iptw_re
)

### Random-Effect (RE) -----------------------------------------
# Random-effects outcome model with PS weights from single-level PS model.
out_slre <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effects outcome model with PS weights from fixed-effect PS model.
out_fere <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effects outcome model with PS weights from random-effect PS model.
out_rere <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)

### Random-Effect with Cluster Means (RE-Mean) ---------------------------------
# Random-effect outcome model with cluster means and PS weights from single-level PS model.
out_slre_cm <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + 
    parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effect outcome model with cluster means and PS weights from fixed-effect PS model.
out_fere_cm <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + 
    parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effect outcome model with cluster means and PS weights from random-effect PS model.
out_rere_cm <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    age_w1_sc + sex_w1 + white_w1 + black_w1 + 
    parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)


## Outcome Models (PNDE & TNIE) --------------------------------------

### Single-Level (SL) ------------------------------------------
# Single-level outcome model with interaction term and PS weights from single-level PS model.
out_slsl_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_sl
)

# Single-level outcome model with interaction term and PS weights from fixed-effect PS model.
out_fesl_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_fe
)

# Single-level outcome model with interaction term and PS weights from random-effect PS model.
out_resl_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1",
  data = data,
  weights = iptw_re
)

### Fixed-Effect (FE) ------------------------------------------
# Fixed-effect outcome model with interaction term and PS weights from single-level PS model.
out_slfe_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_sl
)

# Fixed-effect outcome model with interaction term and PS weights from fixed-effect PS model.
out_fefe_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data, 
  weights = iptw_fe
)

# Fixed-effect outcome model with interaction term and PS weights from random-effect PS model.
out_refe_interac <- glm(
  formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)", 
  data = data,
  weights = iptw_re
)

### Random-Effect (RE) -----------------------------------------
# Random-effects outcome model with interaction term and PS weights from single-level PS model.
out_slre_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effects outcome model with interaction term and PS weights from fixed-effect PS model.
out_fere_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effects outcome model with interaction term and PS weights from random-effect PS model.
out_rere_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)

### Random-Effect with Cluster Means (RE-Mean) ---------------------------------
# Random-effect outcome model with cluster means, interaction term, and PS weights from single-level PS model.
out_slre_cm_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_sl", "L2weight")
)

# Random-effect outcome model with cluster means, interaction term, and PS weights from fixed-effect PS model.
out_fere_cm_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_fe", "L2weight")
)

# Random-effect outcome model with cluster means, interaction term, and PS weights from random-effect PS model.
out_rere_cm_interac <- WeMix::mix(
  formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
    cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc + 
    selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
    white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
  data = data,
  weights = c("iptw_re", "L2weight")
)



# Display Estimates -------------------------------------------------------

# Define conditions
conditions <- c("slsl", "fesl", "resl", 
                "slfe", "fefe", "refe", 
                "slre", "fere", "rere", 
                "slre_cm", "fere_cm", "rere_cm")

# Extract TNDE estimates
TNDE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond, "_interac")
  # model_name <- paste0("out_", cond)
  # summary(get(out_model_name))$coef["sportPartic_w1", "Estimate"]
  summary(get(out_model_name))$coef["sportPartic_w1", "Estimate"] +
    summary(get(out_model_name))$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"] *
    (summary(get(med_model_name))$coef["(Intercept)", "Estimate"] +
       summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"])
})

# Extract PNIE estimates
PNIE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond, "_interac")
  # out_model_name <- paste0("out_", cond)
  summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"] * 
    summary(get(out_model_name))$coef["selfEst_w3_sc", "Estimate"]
})

# Extract PNDE estimates
PNDE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond, "_interac")
  summary(get(out_model_name))$coef["sportPartic_w1", "Estimate"] + 
    summary(get(out_model_name))$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"] * 
    summary(get(med_model_name))$coef["(Intercept)", "Estimate"] 
})

# Extract TNIE estimates
TNIE <- sapply(conditions, function(cond) {
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond, "_interac")
  summary(get(med_model_name))$coef["sportPartic_w1", "Estimate"] * 
    (summary(get(out_model_name))$coef["selfEst_w3_sc", "Estimate"] + 
    summary(get(out_model_name))$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
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
#       cond      TNDE      PNDE        PNIE        TNIE
# 1     slsl 0.2321494 0.2311803 -0.04519194 -0.04422288
# 2     fesl 1.0025387 0.9966576 -0.12429266 -0.11841160
# 3     resl 0.9065915 0.9023210 -0.09379753 -0.08952703
# 4     slfe 0.4792298 0.4759510 -0.11979506 -0.11651626
# 5     fefe 1.6211455 1.6123121 -0.12406800 -0.11523462
# 6     refe 1.4425614 1.4348636 -0.11983126 -0.11213347
# 7     slre 0.3762131 0.3740996 -0.08527099 -0.08315750
# 8     fere 1.4563319 1.4485534 -0.12320456 -0.11542608
# 9     rere 1.2580573 1.2520506 -0.10488250 -0.09887577
# 10 slre_cm 0.3489953 0.3465779 -0.11219647 -0.10977906
# 11 fere_cm 1.3915754 1.3839682 -0.12813238 -0.12052516
# 12 rere_cm 1.2172823 1.2109209 -0.11767835 -0.11131691

# Save the final dataframe of results
write_rds(results_DF, file = "Application/Output/Estimates/Effect-Estimates_noCIs.rds")

# Save mediator and outcome models
sapply(conditions, function(cond) {
  # Collect object names
  med_model_name <- paste0("med_", cond)
  out_model_name <- paste0("out_", cond)
  outWithInterac_model_name <- paste0("out_", cond, "_interac")
  # Save objects
  write_rds(get(med_model_name), file = paste0("Application/Output/mediator-and-outcome-models/", med_model_name, ".rds"))
  write_rds(get(out_model_name), file = paste0("Application/Output/mediator-and-outcome-models/", out_model_name, ".rds"))
  write_rds(get(outWithInterac_model_name), file = paste0("Application/Output/mediator-and-outcome-models/", outWithInterac_model_name, ".rds"))
})

# object_names <- ls(pattern = "^(out_|med_)")
# for (obj_name in object_names) {
#   file_path <- file.path("Application/Output/mediator-and-outcome-models", paste0(obj_name, ".rds"))
#   saveRDS(get(obj_name), file = file_path)
# }


## Clean Environment -------------------------------------------------------
# Remove all objects from the environment except for 'data', 'results_DF', and functions
# rm(list = setdiff(ls(), c("data", "results_DF", lsf.str())))







# Monte Carlo Confidence Intervals (CI) -----------------------------------
# This section computes Monte Carlo confidence intervals for each effect 
# using each PS model (SL, FE, & RE) and mediator/outcome model (SL, FE, RE, & RE-Mean). 

# Import all mediator and outcome models

# List all .rds files in the folder (with full paths)
rds_files <- list.files(path = "Application/Output/mediator-and-outcome-models", pattern = "\\.rds$", full.names = TRUE)

# Loop through each file and assign its content to an object in the global environment
for (file in rds_files) {
  # Get the object name from the file name by removing the folder path and ".rds" extension
  obj_name <- tools::file_path_sans_ext(basename(file))
  # Read the RDS file and assign it to the object name in the global environment
  assign(obj_name, readRDS(file), envir = .GlobalEnv)
}


## TNDE & PNIE -------------------------------------------------------------
# This subsection focuses on TNDE (Total Natural Direct Effect) and PNIE (Pure Natural Indirect Effect) across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_slsl, 
  outcome_fit = out_slsl_interac, # out_slsl, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fesl_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_fesl, 
  outcome_fit = out_fesl_interac, # out_fesl, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# RE PS Model
resl_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_resl, 
  outcome_fit = out_resl_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_slfe, 
  outcome_fit = out_slfe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fefe_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_fefe, 
  outcome_fit = out_fefe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# RE PS Model
refe_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_refe, 
  outcome_fit = out_refe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
slre_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_slre, 
  outcome_fit = out_slre_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "RE", 
  outcome_model_type = "RE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fere_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_fere, 
  outcome_fit = out_fere_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "RE", 
  outcome_model_type = "RE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# RE PS Model
rere_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_rere, 
  outcome_fit = out_rere_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "RE", 
  outcome_model_type = "RE", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
slre_cm_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_slre_cm, 
  outcome_fit = out_slre_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "RE-Mean", 
  outcome_model_type = "RE-Mean", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fere_cm_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_fere_cm, 
  outcome_fit = out_fere_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "RE-Mean", 
  outcome_model_type = "RE-Mean", 
  effect_type = "PNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# RE PS Model
rere_cm_ci_PNIE <- monteCarloCIb(
  mediator_fit = med_rere_cm, 
  outcome_fit = out_rere_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "RE-Mean", 
  outcome_model_type = "RE-Mean", 
  effect_type = "PNIE", 
  output_type = "dataframe",
  use_joint = TRUE
)

## PNDE & TNIE -------------------------------------------------------------
# This subsection focuses on PNDE (Pure Natural Direct Effect) and TNIE (Total Natural Indirect Effect)
# effects across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_slsl,
  outcome_fit = out_slsl_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "TNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fesl_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_fesl, 
  outcome_fit = out_fesl_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "TNIE", 
  output_type = "dataframe",
  use_joint = TRUE
)

# RE PS Model
resl_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_resl, 
  outcome_fit = out_resl_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "SL", 
  outcome_model_type = "SL", 
  effect_type = "TNIE", 
  output_type = "dataframe",
  use_joint = TRUE
)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_slfe, 
  outcome_fit = out_slfe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

# FE PS Model
fefe_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_fefe, 
  outcome_fit = out_fefe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

# RE PS Model
refe_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_refe, 
  outcome_fit = out_refe_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "FE", 
  outcome_model_type = "FE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
slre_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_slre, 
  outcome_fit = out_slre_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "RE",
  outcome_model_type = "RE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

# FE PS Model
fere_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_fere, 
  outcome_fit = out_fere_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "RE",
  outcome_model_type = "RE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

# RE PS Model
rere_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_rere, 
  outcome_fit = out_rere_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "RE",
  outcome_model_type = "RE", 
  effect_type = "TNIE", 
  output_type = "dataframe"
)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
slre_cm_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_slre_cm, 
  outcome_fit = out_slre_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "SL", 
  mediator_model_type = "RE-Mean",
  outcome_model_type = "RE-Mean", 
  effect_type = "TNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# FE PS Model
fere_cm_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_fere_cm, 
  outcome_fit = out_fere_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "FE", 
  mediator_model_type = "RE-Mean",
  outcome_model_type = "RE-Mean", 
  effect_type = "TNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)

# RE PS Model
rere_cm_ci_TNIE <- monteCarloCIb(
  mediator_fit = med_rere_cm, 
  outcome_fit = out_rere_cm_interac, 
  n_MC = 1000, 
  seed_MC = 4561, 
  PS_model_type = "RE", 
  mediator_model_type = "RE-Mean",
  outcome_model_type = "RE-Mean", 
  effect_type = "TNIE", 
  output_type = "dataframe", 
  use_joint = TRUE
)


## Display & Save Monte Carlo CI Estimates -------------------------------------------------------

# # Create PNIE dataframe
# pnie_df <- rbind(
#   slsl_ci_PNIE,
#   fesl_ci_PNIE,
#   resl_ci_PNIE,
#   slfe_ci_PNIE,
#   fefe_ci_PNIE,
#   refe_ci_PNIE,
#   slre_ci_PNIE,
#   fere_ci_PNIE,
#   rere_ci_PNIE,
#   slre_cm_ci_PNIE,
#   fere_cm_ci_PNIE,
#   rere_cm_ci_PNIE
# )
# # rename columns & narrow df 
# pnie_df <- pnie_df |> 
#   dplyr::select("PS_model_type":"NDE_UL") |> 
#   dplyr::rename(PS = PS_model_type, 
#                 outcome = outcome_model_type, 
#                 PNIE = NIE_est, 
#                 PNIE_LL = NIE_LL, 
#                 PNIE_UL = NIE_UL, 
#                 TNDE = NDE_est, 
#                 TNDE_LL = NDE_LL, 
#                 TNDE_UL = NDE_UL) |> 
#   dplyr::select(!effect_type)
# # Create TNIE dataframe
# tnie_df <- rbind(
#   slsl_ci_TNIE,
#   fesl_ci_TNIE,
#   resl_ci_TNIE,
#   slfe_ci_TNIE,
#   fefe_ci_TNIE,
#   refe_ci_TNIE,
#   slre_ci_TNIE,
#   fere_ci_TNIE,
#   rere_ci_TNIE,
#   slre_cm_ci_TNIE,
#   fere_cm_ci_TNIE,
#   rere_cm_ci_TNIE
# )
# # rename columns & narrow df 
# tnie_df <- tnie_df |> 
#   dplyr::select("PS_model_type":"NDE_UL") |> 
#   dplyr::rename(PS = PS_model_type, 
#                 outcome = outcome_model_type, 
#                 TNIE = NIE_est, 
#                 TNIE_LL = NIE_LL, 
#                 TNIE_UL = NIE_UL, 
#                 PNDE = NDE_est, 
#                 PNDE_LL = NDE_LL, 
#                 PNDE_UL = NDE_UL) |> 
#   dplyr::select(!effect_type)
# 
# # Combine DFs effects 
# ci_df <- pnie_df |> 
#   left_join(tnie_df, by = c("PS", "outcome"))
# 
# 
# # Save the final dataframe of results
# write_rds(ci_df, file = "Application/Output/Estimates/Effect-Estimates_monte-carlo-CIs.rds")
# 

# Create PNIE dataframe
pnie_df <- rbind(
  slsl_ci_PNIE,
  fesl_ci_PNIE,
  resl_ci_PNIE,
  slfe_ci_PNIE,
  fefe_ci_PNIE,
  refe_ci_PNIE,
  slre_ci_PNIE,
  fere_ci_PNIE,
  rere_ci_PNIE,
  slre_cm_ci_PNIE,
  fere_cm_ci_PNIE,
  rere_cm_ci_PNIE
)
# Create TNIE dataframe
tnie_df <- rbind(
  slsl_ci_TNIE,
  fesl_ci_TNIE,
  resl_ci_TNIE,
  slfe_ci_TNIE,
  fefe_ci_TNIE,
  refe_ci_TNIE,
  slre_ci_TNIE,
  fere_ci_TNIE,
  rere_ci_TNIE,
  slre_cm_ci_TNIE,
  fere_cm_ci_TNIE,
  rere_cm_ci_TNIE
)

# Combine relevant columns from each dataset & adjust names 
monteCI_df <- pnie_df |> 
  dplyr::select("analysisCond":"PNIE_UCL") |> 
  left_join(dplyr::select(tnie_df, "analysisCond":"PNDE_UCL", "TNIE_est":"TNIE_UCL")) |> 
  rename_with(~ gsub("_LCL$", "_LL", .)) |> 
  rename_with(~ gsub("_UCL$", "_UL", .)) |> 
  rename_with(~ gsub("_est$", "", .))

# Save the final dataframe of results
write_rds(monteCI_df, file = "Application/Output/Estimates/Effect-Estimates_monte-carlo-CIs.rds")

# Save monte carlo CI objects
# Define the list of file names
list_names <- c("slsl_ci_TNIE", "slsl_ci_PNIE", 
                "fesl_ci_TNIE", "fesl_ci_PNIE",
                "resl_ci_TNIE", "resl_ci_PNIE", 
                "slfe_ci_TNIE", "slfe_ci_PNIE", 
                "fefe_ci_TNIE", "fefe_ci_PNIE",
                "refe_ci_TNIE", "refe_ci_PNIE", 
                "slre_ci_TNIE", "slre_ci_PNIE",
                "fere_ci_TNIE", "fere_ci_PNIE",
                "rere_ci_TNIE", "rere_ci_PNIE",
                "slre_cm_ci_TNIE", "slre_cm_ci_PNIE",
                "fere_cm_ci_TNIE", "fere_cm_ci_PNIE",
                "rere_cm_ci_TNIE", "rere_cm_ci_PNIE")
# Save outputs 
sapply(list_names, function(lstname) {
  write_rds(get(lstname), file = paste0("Application/Output/Temp-monte-carlo-CIs/", lstname, ".rds"))
})





################################## END ############################################



