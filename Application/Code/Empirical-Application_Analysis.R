################################################################################
##################### QP Empirical Application - Analysis ######################
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
# Last Updated: 08/13/2024 
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
source("Application/Functions/bootstrap_ci_paral_2.R")
source("Application/Functions/bootstrap_ci_re_paral_2.R")
source("Application/Functions/bootstrap_ci_re_mean_paral.R")



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

# Propensity Score (PS) & Inverse Probability of Treatment Weight (IPTW) Calculation ------
# This section calculates propensity scores and IPTW using three different models: 
# Standard Logistic (SL), Fixed Effects (FE), and Random Effects (RE).



# PS & IPTW Calculation ---------------------------------------------------------------
# This script calculates propensity scores and Inverse Probability of Treatment Weights (IPTW)
# for three models: SL (Single-Level), FE (Fixed-Effect), and RE (Random-Effect).

### SL Propensity Score Model -----------------------------------------------
# Single-level logistic regression model for propensity score estimation
psmod_sl <- glm(
  formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc",
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
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc +
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
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 + selfEst_w1_sc + 
                        (1 | CLUSTER2)",
  family = "binomial", 
  data = data
)

# Predict propensity scores and log-odds for the RE model
data$ps_re <- predict(psmod_re, type = "response") # Propensity scores
data$ps_re_logit <- predict(psmod_re, type = "link") # Log-odds

# Calculate IPTW for the RE model
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))



# Truncate Non-overlap Cases ---------------------------------------------------------

### SL (Single-Level) ------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_sl < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_sl > 0.99)))

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

# Adjust IPTW values to the 1st and 99th percentile thresholds
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < first_percentile, first_percentile, 
                          ifelse(iptw_sl > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_sl)))

### FE (Fixed-Effect) -------------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_fe < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_fe > 0.99)))

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

# Adjust IPTW values to the 1st and 99th percentile thresholds 
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < first_percentile, first_percentile, 
                          ifelse(iptw_fe > ninety_ninth_percentile, ninety_ninth_percentile, 
                                 iptw_fe)))

### RE (Random-Effect) -------------------------------------------------------------
# Identify and count instances of extreme PS values
paste0("Number of PSs < 0.01: ", sum(I(data$ps_re < 0.01)), 
       "; Number of PSs > 0.99: ", sum(I(data$ps_re > 0.99)))

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
                "healthInsur_w1", "selfEst_w1_sc")

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
                  "healthInsur_w1", "age_w1_sc", "sex_w1", 
                  "selfEst_w1_sc", "feelings_w1_sc", "parentalEdu_w1_sc")

# Define new labels for the y-axis to improve readability
new_labels <- c("Race: Black", "Race: White", "Family Structure", 
                "Health Insurance \n Coverage Gap", "Age", "Sex", 
                "Self-Esteem Score", "Feelings Scale Score", "Parental Education")


## Visualization: Love Plot for Paper --------------------------------

# Load Times New Roman font for publication-quality visuals
# loadfonts()

# Save the Love Plot visualization for the paper as a PDF
pdf("Application/Output/Covariate-Balance_QP-Doc.pdf")

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
  scale_color_manual(values = c("Unweighted" = "#1f77b4",  # Blue for unweighted
                                "Single-Level" = "#2ca02c",  # Green for single-level weighting
                                "Fixed-Effect" = "#ff7f0e",  # Orange for fixed-effect weighting
                                "Random-Effect" = "#9467bd"),  # Purple for random-effect weighting
                     name = NULL) +
  scale_shape_manual(values = c("Unweighted" = 16,  # Circle shape for unweighted
                                "Single-Level" = 17,  # Triangle shape for single-level
                                "Fixed-Effect" = 15,  # Square shape for fixed-effect
                                "Random-Effect" = 18),  # Diamond shape for random-effect
                     name = NULL) +
  scale_y_discrete(labels = new_labels)  # Apply new labels to the y-axis

# Close the PDF device
dev.off()

# Save the Love Plot visualization for the paper as a PNG
ggsave(filename = "Application/Output/Covariate-Balance_QP-Doc.pdf", 
       plot = last_plot(), 
       width = 6, 
       height = 7, 
       units = "in", 
       dpi = 300)

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

# cond       TNDE       PNDE        PNIE          TNIE
# 1     slsl -0.3317336 -0.3318426 -0.01230346  0.0247090071
# 2     fesl -0.2558156 -0.2563010 -0.09886261 -0.0261259555
# 3     resl -0.2959685 -0.2961819 -0.06592642  0.0012385894
# 4     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906
# 5     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050
# 6     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374
# 7     slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456
# 8     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596
# 10 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108
# 11 fere_cm -0.2425750 -0.2431934 -0.10177153 -0.0098265525
# 12 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984



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
  slre_ci_PNIE <- bootstrap_ci_re_paral_2(iterations = 1750, 
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
  rere_ci_PNIE <- bootstrap_ci_re_paral_2(iterations = 1750, 
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

#### Random-Effect with Cluster Means (RE-Mean) Med/Out Models ---------------------------------------
# SL PS Model 
execution_time_slre_cm <- system.time({ # Track computation time 
  slre_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_sl, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "PNIE")
})
# Print the execution time
print(execution_time_slre_cm)
#  

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_slre_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", slre_cm_ci_PNIE$mediator_converged_count, 
       " (", (slre_cm_ci_PNIE$mediator_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", slre_cm_ci_PNIE$outcome_converged_count, 
       " (", (slre_cm_ci_PNIE$outcome_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", slre_cm_ci_PNIE$both_converged_count, 
       " (", (slre_cm_ci_PNIE$both_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)")
# 

# FE PS Model
execution_time_fere_cm <- system.time({ # Track computation time 
  fere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_fe, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "PNIE")
})
# Print the execution time
print(execution_time_fere_cm)
#  

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_fere_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", fere_cm_ci_PNIE$mediator_converged_count, 
       " (", (fere_cm_ci_PNIE$mediator_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", fere_cm_ci_PNIE$outcome_converged_count, 
       " (", (fere_cm_ci_PNIE$outcome_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", fere_cm_ci_PNIE$both_converged_count, 
       " (", (fere_cm_ci_PNIE$both_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)")
# 

# RE PS Model
execution_time_rere_cm <- system.time({ # Track computation time 
  rere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_re, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "PNIE")
})
# Print the execution time
print(execution_time_rere_cm)
#  

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_rere_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", rere_cm_ci_PNIE$mediator_converged_count, 
       " (", (rere_cm_ci_PNIE$mediator_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", rere_cm_ci_PNIE$outcome_converged_count, 
       " (", (rere_cm_ci_PNIE$outcome_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", rere_cm_ci_PNIE$both_converged_count, 
       " (", (rere_cm_ci_PNIE$both_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)")
# 






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
  slre_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 2000, 
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
  fere_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 2000, 
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
  rere_ci_TNIE <- bootstrap_ci_re_paral_2(iterations = 2000, 
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

#### Random-Effect with Cluster Means (RE-Mean) Med/Out Models ---------------------------------------
# SL PS Model 
execution_time_slre_cm <- system.time({ # Track computation time 
  slre_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_sl, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_slre_cm)
#  

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_slre_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", slre_cm_ci_TNIE$mediator_converged_count, 
       " (", (slre_cm_ci_TNIE$mediator_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", slre_cm_ci_TNIE$outcome_converged_count, 
       " (", (slre_cm_ci_TNIE$outcome_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", slre_cm_ci_TNIE$both_converged_count, 
       " (", (slre_cm_ci_TNIE$both_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)")
# 

# FE PS Model
execution_time_fere_cm <- system.time({ # Track computation time 
  fere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_fe, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_fere_cm)
# 

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_fere_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", fere_cm_ci_TNIE$mediator_converged_count, 
       " (", (fere_cm_ci_TNIE$mediator_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", fere_cm_ci_TNIE$outcome_converged_count, 
       " (", (fere_cm_ci_TNIE$outcome_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", fere_cm_ci_TNIE$both_converged_count, 
       " (", (fere_cm_ci_TNIE$both_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)")
# 

# RE PS Model
execution_time_rere_cm <- system.time({ # Track computation time 
  rere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(iterations = 2000, 
                                          iptw = iptw_re, 
                                          data = data, 
                                          cores = 6, 
                                          core_seeds = c(4561:4566), 
                                          effect_type = "TNIE")
})
# Print the execution time
print(execution_time_rere_cm)
#  

# Print the elapsed time specifically
cat("Elapsed time:", execution_time_rere_cm["elapsed"], "seconds\n")
# 

# Print convergence statistics
paste0("Number of converged mediator models: ", rere_cm_ci_TNIE$mediator_converged_count, 
       " (", (rere_cm_ci_TNIE$mediator_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of converged outcome models: ", rere_cm_ci_TNIE$outcome_converged_count, 
       " (", (rere_cm_ci_TNIE$outcome_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)")
paste0("Number of iterations with both models converged: ", rere_cm_ci_TNIE$both_converged_count, 
       " (", (rere_cm_ci_TNIE$both_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)")
# 






################################################################################

# [NEXT TO WORK ON] -------------------------------------------------------





################################################################################


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


