################################################################################
##################### QP Empirical Application - Analysis ######################
################################################################################

############################ Script Description ################################
#
# Author: Cameron
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
# estimated effects to facilitate interpretation of the results.
# 
# Last Updated: 2025-03-10 08/15/2024 
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
#       cond       TNDE       PNDE        PNIE        TNIE
# 1     slsl -0.3259836 -0.3259562 -0.04471103 -0.01584023
# 2     fesl -0.2588860 -0.2591964 -0.12127187 -0.05814495
# 3     resl -0.2987412 -0.2987693 -0.09161688 -0.03054913
# 4     slfe -0.2325015 -0.2327627 -0.11816059 -0.08187151
# 5     fefe -0.2510769 -0.2514684 -0.11950564 -0.02166266
# 6     refe -0.2562645 -0.2566047 -0.11588065 -0.02853428
# 7     slre -0.2719431 -0.2720481 -0.08421877 -0.05117338
# 8     fere -0.2562957 -0.2566773 -0.11918980 -0.03342662
# 9     rere -0.2765615 -0.2767321 -0.10180340 -0.02434469
# 10 slre_cm -0.2243275 -0.2245039 -0.11099540 -0.08241447
# 11 fere_cm -0.2490276 -0.2494285 -0.12421280 -0.04336970
# 12 rere_cm -0.2557420 -0.2560127 -0.11442310 -0.04131232

## Clean Environment -------------------------------------------------------
# Remove all objects from the environment except for 'data', 'results_DF', and functions
# rm(list = setdiff(ls(), c("data", "results_DF", lsf.str())))



# Bootstrap Confidence Intervals (CI) -------------------------------------
# This section conducts the bootstrapping for confidence intervals for each effect 
# using each PS model (SL, FE, & RE) and mediator/outcome model (SL, FE, RE, & RE-Mean). 
# Results are saved and convergence statistics are printed.

# WARNING: Both RE and RE-Mean calculations require extended processing times.

## TNDE & PNIE -------------------------------------------------------------
# This subsection focuses on TNDE (Total Natural Direct Effect) and PNIE (Pure Natural Indirect Effect) across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify mediator & outcome model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "PNIE"      # Effect type: PNIE
)
saveRDS(slsl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slsl_ci_PNIE.rds")
rm(slsl_ci_PNIE)

# FE PS Model
fesl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fesl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fesl_ci_PNIE.rds")
rm(fesl_ci_PNIE)

# RE PS Model
resl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(resl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/resl_ci_PNIE.rds")
rm(resl_ci_PNIE)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",   
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(slfe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slfe_ci_PNIE.rds")
rm(slfe_ci_PNIE)

# FE PS Model
fefe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fefe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fefe_ci_PNIE.rds")
rm(fefe_ci_PNIE)

# RE PS Model
refe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(refe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/refe_ci_PNIE.rds")
rm(refe_ci_PNIE)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
execution_time <- system.time({ # Track computation time 
  slre_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,  
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(slre_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slre_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_ci_PNIE$mediator_converged_count,
    " (", (slre_ci_PNIE$mediator_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_ci_PNIE$outcome_converged_count,
    " (", (slre_ci_PNIE$outcome_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_ci_PNIE$both_converged_count,
    " (", (slre_ci_PNIE$both_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(slre_ci_PNIE)
# Elapsed time: 6652.806 seconds ( 111 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# FE PS Model
execution_time <- system.time({ 
  fere_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(fere_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fere_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_ci_PNIE$mediator_converged_count,
    " (", (fere_ci_PNIE$mediator_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_ci_PNIE$outcome_converged_count,
    " (", (fere_ci_PNIE$outcome_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_ci_PNIE$both_converged_count,
    " (", (fere_ci_PNIE$both_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(fere_ci_PNIE)
# Elapsed time: 7617.093 seconds ( 127 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# RE PS Model
execution_time <- system.time({ 
  rere_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(rere_ci_PNIE, file = "Application/Output/Bootstrap_Temp/rere_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_ci_PNIE$mediator_converged_count,
    " (", (rere_ci_PNIE$mediator_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_ci_PNIE$outcome_converged_count,
    " (", (rere_ci_PNIE$outcome_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_ci_PNIE$both_converged_count,
    " (", (rere_ci_PNIE$both_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(rere_ci_PNIE)
# Elapsed time: 10816.36 seconds ( 180 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
execution_time <- system.time({ 
  slre_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(slre_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slre_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_cm_ci_PNIE$mediator_converged_count,
    " (", (slre_cm_ci_PNIE$mediator_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_cm_ci_PNIE$outcome_converged_count,
    " (", (slre_cm_ci_PNIE$outcome_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_cm_ci_PNIE$both_converged_count,
    " (", (slre_cm_ci_PNIE$both_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(slre_cm_ci_PNIE)
# Elapsed time: 6615.656 seconds ( 110 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# FE PS Model
execution_time <- system.time({ 
  fere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(fere_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fere_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_cm_ci_PNIE$mediator_converged_count,
    " (", (fere_cm_ci_PNIE$mediator_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_cm_ci_PNIE$outcome_converged_count,
    " (", (fere_cm_ci_PNIE$outcome_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_cm_ci_PNIE$both_converged_count,
    " (", (fere_cm_ci_PNIE$both_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(fere_cm_ci_PNIE)
# Elapsed time: 6152.866 seconds ( 103 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# RE PS Model
execution_time <- system.time({ 
  rere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(rere_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/rere_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_cm_ci_PNIE$mediator_converged_count,
    " (", (rere_cm_ci_PNIE$mediator_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_cm_ci_PNIE$outcome_converged_count,
    " (", (rere_cm_ci_PNIE$outcome_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_cm_ci_PNIE$both_converged_count,
    " (", (rere_cm_ci_PNIE$both_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(rere_cm_ci_PNIE)
# Elapsed time: 5956.323 seconds ( 99 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

## PNDE & TNIE -------------------------------------------------------------
# This subsection focuses on PNDE (Pure Natural Direct Effect) and TNIE (Total Natural Indirect Effect)
# effects across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify PS model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "TNIE"      # Effect type: TNIE
)
saveRDS(slsl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slsl_ci_TNIE.rds")
rm(slsl_ci_TNIE)

# FE PS Model
fesl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fesl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fesl_ci_TNIE.rds")
rm(fesl_ci_TNIE)

# RE PS Model
resl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(resl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/resl_ci_TNIE.rds")
rm(resl_ci_TNIE)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(slfe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slfe_ci_TNIE.rds")
rm(slfe_ci_TNIE)

# FE PS Model
fefe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fefe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fefe_ci_TNIE.rds")
rm(fefe_ci_TNIE)

# RE PS Model
refe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(refe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/refe_ci_TNIE.rds")
rm(refe_ci_TNIE)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
execution_time <- system.time({ 
  slre_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(slre_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slre_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_ci_TNIE$mediator_converged_count,
    " (", (slre_ci_TNIE$mediator_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_ci_TNIE$outcome_converged_count,
    " (", (slre_ci_TNIE$outcome_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_ci_TNIE$both_converged_count,
    " (", (slre_ci_TNIE$both_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(slre_ci_TNIE)
# Elapsed time: 6310.98 seconds ( 105 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

# FE PS Model
execution_time <- system.time({ 
  fere_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(fere_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fere_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_ci_TNIE$mediator_converged_count,
    " (", (fere_ci_TNIE$mediator_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_ci_TNIE$outcome_converged_count,
    " (", (fere_ci_TNIE$outcome_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_ci_TNIE$both_converged_count,
    " (", (fere_ci_TNIE$both_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(fere_ci_TNIE)
# Elapsed time: 6343.685 seconds ( 106 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

# RE PS Model
execution_time <- system.time({ 
  rere_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(rere_ci_TNIE, file = "Application/Output/Bootstrap_Temp/rere_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_ci_TNIE$mediator_converged_count,
    " (", (rere_ci_TNIE$mediator_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_ci_TNIE$outcome_converged_count,
    " (", (rere_ci_TNIE$outcome_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_ci_TNIE$both_converged_count,
    " (", (rere_ci_TNIE$both_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(rere_ci_TNIE)
# Elapsed time: 10578.82 seconds ( 176 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
execution_time <- system.time({ 
  slre_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(slre_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slre_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_cm_ci_TNIE$mediator_converged_count,
    " (", (slre_cm_ci_TNIE$mediator_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_cm_ci_TNIE$outcome_converged_count,
    " (", (slre_cm_ci_TNIE$outcome_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_cm_ci_TNIE$both_converged_count,
    " (", (slre_cm_ci_TNIE$both_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(slre_cm_ci_TNIE)
# Elapsed time: 9869.804 seconds ( 164 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD:
# Elapsed time: 4365.37 seconds ( 73 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)

# FE PS Model
execution_time <- system.time({ 
  fere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(fere_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fere_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_cm_ci_TNIE$mediator_converged_count,
    " (", (fere_cm_ci_TNIE$mediator_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_cm_ci_TNIE$outcome_converged_count,
    " (", (fere_cm_ci_TNIE$outcome_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_cm_ci_TNIE$both_converged_count,
    " (", (fere_cm_ci_TNIE$both_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(fere_cm_ci_TNIE)
# Elapsed time: 5830.462 seconds ( 97 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD: 
# Elapsed time: 3742.066 seconds ( 62 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)

# RE PS Model
execution_time <- system.time({ 
  rere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(rere_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/rere_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_cm_ci_TNIE$mediator_converged_count,
    " (", (rere_cm_ci_TNIE$mediator_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_cm_ci_TNIE$outcome_converged_count,
    " (", (rere_cm_ci_TNIE$outcome_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_cm_ci_TNIE$both_converged_count,
    " (", (rere_cm_ci_TNIE$both_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(rere_cm_ci_TNIE)
# Elapsed time: 7174.466 seconds ( 120 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD: 
# Elapsed time: 3800.283 seconds ( 63 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)



# Store & Join Results ----------------------------------------------------
# This section handles the calculation and storage of confidence intervals (CIs). 

#### RE Mediator/Outcome Models CI -----------------------------------------
# This subsection focuses on obtaining CIs for RE & RE-Mean models and merging them with effect estimates.

###### Obtain 1,000 completed iterations -----------------------------------
# Function to get non-NA pairs (i.e., first 1,000 completed iterations)
get_non_na_pairs <- function(direct, indirect, n = 1000) {
  combined <- data.frame(direct = direct, indirect = indirect)  # Combine vectors into a dataframe
  combined <- na.omit(combined)  # Remove rows with any NA values
  return(head(combined, n))  # Return the first n rows (or all if less than n)
}

# Define the file paths for the RE models
files <- list(
  "Application/Output/Bootstrap_Temp/slre_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_ci_PNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_cm_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_cm_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_cm_ci_PNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_ci_TNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_cm_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_cm_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_cm_ci_TNIE.rds"
)

# Function to process each file, apply get_non_na_pairs, and remove the original data
process_file <- function(file_path, n = 1000) {
  data <- readRDS(file_path)  # Load the RDS file
  result <- get_non_na_pairs(data$direct_effects, data$indirect_effects, n = n)  # Apply function
  rm(data)  # Remove the original data to free up space
  return(result)  # Return the processed data
}

# Process each file and store results in corresponding dataframes
slre_ci_PNIE_DF <- process_file(files[[1]])
fere_ci_PNIE_DF <- process_file(files[[2]])
rere_ci_PNIE_DF <- process_file(files[[3]])

slre_cm_ci_PNIE_DF <- process_file(files[[4]])
fere_cm_ci_PNIE_DF <- process_file(files[[5]])
rere_cm_ci_PNIE_DF <- process_file(files[[6]])

slre_ci_TNIE_DF <- process_file(files[[7]])
fere_ci_TNIE_DF <- process_file(files[[8]])
rere_ci_TNIE_DF <- process_file(files[[9]])

slre_cm_ci_TNIE_DF <- process_file(files[[10]])
fere_cm_ci_TNIE_DF <- process_file(files[[11]])
rere_cm_ci_TNIE_DF <- process_file(files[[12]])

# Remove 'files' vector to save space 
rm(files)

###### Store RE Med/Outcome Model CIs --------------------------------------
# This block stores the calculated CIs for RE models in a structured dataframe.

# Create an empty dataframe to store the results
results_DF_RE <- data.frame(
  cond = c("slre", "fere", "rere", 
           "slre_cm", "fere_cm", "rere_cm"),
  PNIE_LL = numeric(6),
  PNIE_UL = numeric(6),
  TNDE_LL = numeric(6),
  TNDE_UL = numeric(6),
  TNIE_LL = numeric(6),
  TNIE_UL = numeric(6),
  PNDE_LL = numeric(6),
  PNDE_UL = numeric(6),
  stringsAsFactors = FALSE
)

# Lists of processed dataframes
df_list_PNIE <- list(slre_ci_PNIE_DF, fere_ci_PNIE_DF, rere_ci_PNIE_DF, 
                     slre_cm_ci_PNIE_DF, fere_cm_ci_PNIE_DF, rere_cm_ci_PNIE_DF)
df_list_TNIE <- list(slre_ci_TNIE_DF, fere_ci_TNIE_DF, rere_ci_TNIE_DF, 
                     slre_cm_ci_TNIE_DF, fere_cm_ci_TNIE_DF, rere_cm_ci_TNIE_DF)

# Calculate CIs and fill the dataframe
for (i in 1:6) {
  results_DF_RE[i, c("PNIE_LL", "PNIE_UL")] <- quantile(df_list_PNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("TNDE_LL", "TNDE_UL")] <- quantile(df_list_PNIE[[i]]$direct, probs = c(0.025, 0.975))
  
  results_DF_RE[i, c("TNIE_LL", "TNIE_UL")] <- quantile(df_list_TNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("PNDE_LL", "PNDE_UL")] <- quantile(df_list_TNIE[[i]]$direct, probs = c(0.025, 0.975))
}

# Display the RE model results
results_DF_RE
#       cond    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL
# 1    slre -0.2792701 0.09994341 -0.5000102 -0.055682131 -0.2819150 0.10238444 -0.5027937 -0.053682580
# 2    fere -0.3426149 0.07165350 -0.5128617 -0.019704421 -0.3590603 0.07522094 -0.5137532 -0.020466528
# 3    rere -0.3106561 0.08389139 -0.5176295 -0.050033639 -0.3216423 0.09003098 -0.5162965 -0.049349202
# 4 slre_cm -0.3235201 0.08291188 -0.4695924  0.003501834 -0.3145743 0.08801792 -0.4681990  0.003553328
# 5 fere_cm -0.3494602 0.07261279 -0.5102075 -0.008612200 -0.3630953 0.07850820 -0.5110733 -0.010316405
# 6 rere_cm -0.3300292 0.08252020 -0.5048648 -0.026267544 -0.3388155 0.08482095 -0.5054298 -0.026664633

# OLD: 
#       cond    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL
# 1    slre -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.044027732
# 2    fere -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.001785989
# 3    rere -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.034812727
# 4 slre_cm -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.004395259
# 5 fere_cm -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.010775095
# 6 rere_cm -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.014537828

###### Join CI & point estimates for RE med/outcome ------------------------
# This step merges the calculated RE model CIs with existing data.

# Merge RE results with existing results dataframe
results_DF_RE <- merge(results_DF[results_DF$cond %in% c("slre", "fere", "rere", "slre_cm", "fere_cm", "rere_cm"), ], 
                       results_DF_RE)

# Clean up environment (Remove all objects except 'data', 'results_DF', 'results_DF_RE', and functions)
rm(list = setdiff(ls(), c("data", "results_DF", "results_DF_RE", lsf.str())))

#### SL & FE Mediator/Outcome Models CI ------------------------------------
# This section focuses on the calculation and joining of CIs for SL & FE models.

###### Join CI & point estimates for SL & FE med/outcome -------------------
# This block imports data for SL and FE models, calculates CIs, and combines results.

# Define the list of file names
list_names <- c("slsl_ci_TNIE", "slsl_ci_PNIE", 
                "fesl_ci_TNIE", "fesl_ci_PNIE",
                "resl_ci_TNIE", "resl_ci_PNIE", 
                "slfe_ci_TNIE", "slfe_ci_PNIE", 
                "fefe_ci_TNIE", "fefe_ci_PNIE",
                "refe_ci_TNIE", "refe_ci_PNIE")

# Import RDS files and store them in a list of lists
lists <- lapply(list_names, function(name) {
  readRDS(paste0("Application/Output/Bootstrap_Temp/", name, ".rds"))
})

# Name the elements of the list for easier access
names(lists) <- list_names

# Extract 'indirect_ci' and store in dataframe
indirect_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$indirect_ci
  data.frame(
    list_name = list_names[i],
    indirect_ci_LL = ci_values[1],  # Lower bound
    indirect_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Extract 'direct_ci' and store in dataframe
direct_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$direct_ci
  data.frame(
    list_name = list_names[i],
    direct_ci_LL = ci_values[1],  # Lower bound
    direct_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Combine both 'indirect' and 'direct' CIs into a single dataframe
combined_df <- merge(indirect_df, direct_df, by = "list_name")

# Modify dataframe: extract effect and cond labels, and reshape to wide format
combined_df <- combined_df %>%
  mutate(effect = substr(list_name, nchar(list_name) - 3, nchar(list_name))) %>%  # Extract effect label
  mutate(list_name = substr(list_name, 1, 4)) %>%  # Extract cond label
  pivot_longer(cols = c(indirect_ci_LL, indirect_ci_UL, direct_ci_LL, direct_ci_UL), 
               names_to = "variable", values_to = "value") %>%
  unite("variable", effect, variable, sep = "_") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  as.data.frame()

# Update column names for clarity
colnames(combined_df) <- c("cond", 
                           "PNIE_LL", "PNIE_UL", "TNDE_LL", "TNDE_UL", 
                           "TNIE_LL", "TNIE_UL", "PNDE_LL", "PNDE_UL")

# Merge the SL & FE results with the existing results dataframe
results_DF_noRE <- merge(results_DF[1:6, ], combined_df, by = "cond")

# Clean up environment (Remove all objects except 'data', 'results_DF', 'results_DF_RE', 'results_DF_noRE', and functions)
rm(list = setdiff(ls(), c("data", "results_DF", "results_DF_RE", "results_DF_noRE", lsf.str())))

#### Create dataframe of final estimates -----------------------------------
# This final step combines the RE and non-RE results and exports/saves the final estimates.

# Combine RE and non-RE results
results_DF <- rbind(results_DF_noRE, results_DF_RE)

# Display the combined results
# print(results_DF)

# Add labels for PS model & Mediator/Outcome model
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
      endsWith(cond, "re_cm") ~ "Random-Effect with Cluster Means",
      TRUE ~ NA_character_  # Default case
    )
  )

# Display the final dataframe with PS model & Mediator/Outcome labels
print(results_DF)
#       cond       TNDE       PNDE        PNIE        TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL       PNDE_UL            PS                            Model
# 1     fefe -0.2510769 -0.2514684 -0.11950564 -0.02166266 -0.3255577 0.07700916 -0.5060783 -0.005091589 -0.3352327 0.08488664 -0.5072594 -0.0063891720  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2588860 -0.2591964 -0.12127187 -0.05814495 -0.3246048 0.07121469 -0.5103306 -0.003194003 -0.3211064 0.07322311 -0.5112406  0.0009317933  Fixed-Effect                     Single-Level
# 3     refe -0.2562645 -0.2566047 -0.11588065 -0.02853428 -0.3126527 0.07480969 -0.5032962 -0.012264249 -0.3184659 0.07940688 -0.5031124 -0.0093673542 Random-Effect                     Fixed-Effect
# 4     resl -0.2987412 -0.2987693 -0.09161688 -0.03054913 -0.2646908 0.08650299 -0.5280205 -0.062530159 -0.2778580 0.09129741 -0.5305844 -0.0650180918 Random-Effect                     Single-Level
# 5     slfe -0.2325015 -0.2327627 -0.11816059 -0.08187151 -0.3050444 0.07659792 -0.4810479  0.006200502 -0.3079430 0.07819370 -0.4822922  0.0065073266  Single-Level                     Fixed-Effect
# 6     slsl -0.3259836 -0.3259562 -0.04471103 -0.01584023 -0.2067034 0.13493934 -0.5633203 -0.085134771 -0.2155183 0.13529510 -0.5648185 -0.0870138644  Single-Level                     Single-Level
# 7     fere -0.2562957 -0.2566773 -0.11918980 -0.03342662 -0.3426149 0.07165350 -0.5128617 -0.019704421 -0.3590603 0.07522094 -0.5137532 -0.0204665279  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2490276 -0.2494285 -0.12421280 -0.04336970 -0.3494602 0.07261279 -0.5102075 -0.008612200 -0.3630953 0.07850820 -0.5110733 -0.0103164046  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2765615 -0.2767321 -0.10180340 -0.02434469 -0.3106561 0.08389139 -0.5176295 -0.050033639 -0.3216423 0.09003098 -0.5162965 -0.0493492021 Random-Effect                    Random-Effect
# 10 rere_cm -0.2557420 -0.2560127 -0.11442310 -0.04131232 -0.3300292 0.08252020 -0.5048648 -0.026267544 -0.3388155 0.08482095 -0.5054298 -0.0266646326 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2719431 -0.2720481 -0.08421877 -0.05117338 -0.2792701 0.09994341 -0.5000102 -0.055682131 -0.2819150 0.10238444 -0.5027937 -0.0536825803  Single-Level                    Random-Effect
# 12 slre_cm -0.2243275 -0.2245039 -0.11099540 -0.08241447 -0.3235201 0.08291188 -0.4695924  0.003501834 -0.3145743 0.08801792 -0.4681990  0.0035533282  Single-Level Random-Effect with Cluster Means

# OLD: 
#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 1     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050 -0.3007998 0.09570006 -0.5084805  0.004777044 -0.3077056 0.10218566 -0.5085235  0.006911447  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2558156 -0.2563010 -0.09886261 -0.0261259555 -0.2921923 0.08181015 -0.5118420 -0.010035295 -0.2961822 0.08569538 -0.5117962 -0.011418954  Fixed-Effect                     Single-Level
# 3     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374 -0.2797301 0.09890207 -0.5069342 -0.000239732 -0.2829443 0.10403690 -0.5077968  0.002558890 Random-Effect                     Fixed-Effect
# 4     resl -0.2959685 -0.2961819 -0.06592642  0.0012385894 -0.2334787 0.10900586 -0.5444086 -0.064872212 -0.2386974 0.11433888 -0.5457297 -0.064180945 Random-Effect                     Single-Level
# 5     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906 -0.2637974 0.10503455 -0.4788932  0.010653283 -0.2650747 0.10698408 -0.4789265  0.014197515  Single-Level                     Fixed-Effect
# 6     slsl -0.3317336 -0.3318426 -0.01230346  0.0247090071 -0.1747149 0.16068737 -0.5693915 -0.084157900 -0.1768113 0.16116459 -0.5715905 -0.083434503  Single-Level                     Single-Level
# 7     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.001785989  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2425750 -0.2431934 -0.10177153 -0.0098265525 -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.010775095  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.034812727 Random-Effect                    Random-Effect
# 10 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.014537828 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.044027732  Single-Level                    Random-Effect
# 12 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.004395259  Single-Level Random-Effect with Cluster Means

#       cond       TNDE       PNDE        PNIE        TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL       PNDE_UL            PS                            Model
# 1     fefe -0.2510769 -0.2514684 -0.11950564 -0.02166266 -0.3255577 0.07700916 -0.5060783 -0.005091589 -0.3352327 0.08488664 -0.5072594 -0.0063891720  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2588860 -0.2591964 -0.12127187 -0.05814495 -0.3246048 0.07121469 -0.5103306 -0.003194003 -0.3211064 0.07322311 -0.5112406  0.0009317933  Fixed-Effect                     Single-Level
# 3     refe -0.2562645 -0.2566047 -0.11588065 -0.02853428 -0.3126527 0.07480969 -0.5032962 -0.012264249 -0.3184659 0.07940688 -0.5031124 -0.0093673542 Random-Effect                     Fixed-Effect
# 4     resl -0.2987412 -0.2987693 -0.09161688 -0.03054913 -0.2646908 0.08650299 -0.5280205 -0.062530159 -0.2778580 0.09129741 -0.5305844 -0.0650180918 Random-Effect                     Single-Level
# 5     slfe -0.2325015 -0.2327627 -0.11816059 -0.08187151 -0.3050444 0.07659792 -0.4810479  0.006200502 -0.3079430 0.07819370 -0.4822922  0.0065073266  Single-Level                     Fixed-Effect
# 6     slsl -0.3259836 -0.3259562 -0.04471103 -0.01584023 -0.2067034 0.13493934 -0.5633203 -0.085134771 -0.2155183 0.13529510 -0.5648185 -0.0870138644  Single-Level                     Single-Level
# 7     fere -0.2562957 -0.2566773 -0.11918980 -0.03342662 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.0017859895  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2490276 -0.2494285 -0.12421280 -0.04336970 -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.0107750952  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2765615 -0.2767321 -0.10180340 -0.02434469 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.0348127272 Random-Effect                    Random-Effect
# 10 rere_cm -0.2557420 -0.2560127 -0.11442310 -0.04131232 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.0145378283 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2719431 -0.2720481 -0.08421877 -0.05117338 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.0440277319  Single-Level                    Random-Effect
# 12 slre_cm -0.2243275 -0.2245039 -0.11099540 -0.08241447 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.0043952591  Single-Level Random-Effect with Cluster Means

# Save the final dataframe of results
write_rds(results_DF, file = "Application/Output/Estimates/Effect-Estimates.rds")



# Result visuals ----------------------------------------------------------

# TNDE 
## Save visual 
pdf("Application/Output/Visuals/TNDE-Estimates.pdf")
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
ggsave(filename = "Application/Output/Visuals/TNDE-Estimates.png", plot = last_plot())


# PNDE 
## Save visual 
pdf("Application/Output/Visuals/PNDE-Estimates.pdf")
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
ggsave(filename = "Application/Output/Visuals/PNDE-Estimates.png", plot = last_plot())




# TNIE 
## Save visual 
pdf("Application/Output/Visuals/TNIE-Estimates.pdf")
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
ggsave(filename = "Application/Output/Visuals/TNIE-Estimates.png", plot = last_plot())


# PNIE 
## Save visual 
pdf("Application/Output/Visuals/PNIE-Estimates.pdf")
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
ggsave(filename = "Application/Output/Visuals/PNIE-Estimates.png", plot = last_plot())




# Result Table ------------------------------------------------------------

results_DF[, c("PS", "Model", 
               "TNDE", "TNDE_LL", "TNDE_UL", 
               "PNDE", "PNDE_LL", "PNDE_UL", 
               "TNIE", "TNIE_LL", "TNIE_UL", 
               "PNIE", "PNIE_LL", "PNIE_UL")]


