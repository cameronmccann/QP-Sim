################################################################################
############### QP Empirical Application - Confidence Intervals ################
################################################################################

############################ Script Description ################################
#
# Author: 
# 
# Date Created: 2025-06-12
#
#
# Script Description:
#   This R script obtains confidence intervals for the empirical application.
#   It first computes percentile bootstrap confidence intervals (which may take several hours).
#   Then, it computes Monte Carlo confidence intervals. The results are 
#   saved as RDS files. Confidence intervals are displayed in the 
#   'Empirical-Application_Results.R' script.

# 
# Last Updated: 2025-06-12 
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
source("Application/Functions/bootstrapCIb.R")    
# source("Application/Functions/bootstrap_ci_paral_2.R")      # SL / FE (parallel)
# source("Application/Functions/bootstrap_ci_re_paral_2.R")   # RE helpers (parallel)
# source("Application/Functions/bootstrap_ci_re_mean_paral.R")# RE-Mean helpers


# Import Data -----------------------------------------------------
# Load clean dataset 
data <- read_rds(file = "Application/Data/Cleaned/Empirical-Application-Data.rds")


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








# Bootstrap CI ---------------------------------------------------

# ══════════════════════════════
#     WARNING: Both RE and RE-Mean calculations require extended processing times.
# ══════════════════════════════

# create temp folder 
if (!dir.exists("Application/Output/Temp-bootstrap-CIs")) {
  dir.create(path = "Application/Output/Temp-bootstrap-CIs", recursive = TRUE)
}

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
if ("iptw_sl" %in% names(data)) {
  # run bootstrap 
  slsl_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_sl,
    data = data,
    model = "SL",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(slsl_ci, file = "Application/Output/Temp-bootstrap-CIs/slsl_ci.rds")
} else {
  print("The column 'iptw_sl' does not exist in the 'data' data frame.")
}

rm(slsl_ci)

# FE PS Model
if ("iptw_fe" %in% names(data)) {
  # run bootstrap 
  fesl_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_fe,
    data = data,
    model = "SL",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(fesl_ci, file = "Application/Output/Temp-bootstrap-CIs/fesl_ci.rds")
} else {
  print("The column 'iptw_fe' does not exist in the 'data' data frame.")
}

rm(fesl_ci)

# RE PS Model
if ("iptw_re" %in% names(data)) {
  # run bootstrap 
  resl_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_re,
    data = data,
    model = "SL",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(resl_ci, file = "Application/Output/Temp-bootstrap-CIs/resl_ci.rds")
} else {
  print("The column 'iptw_re' does not exist in the 'data' data frame.")
}

rm(resl_ci)


### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
if ("iptw_sl" %in% names(data)) {
  # run bootstrap 
  slfe_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_sl,
    data = data,
    model = "FE",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(slfe_ci, file = "Application/Output/Temp-bootstrap-CIs/slfe_ci.rds")
} else {
  print("The column 'iptw_sl' does not exist in the 'data' data frame.")
}

rm(slfe_ci)

# FE PS Model
if ("iptw_fe" %in% names(data)) {
  # run bootstrap 
  fefe_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_fe,
    data = data,
    model = "FE",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(fefe_ci, file = "Application/Output/Temp-bootstrap-CIs/fefe_ci.rds")
} else {
  print("The column 'iptw_fe' does not exist in the 'data' data frame.")
}

rm(fefe_ci)

# RE PS Model
if ("iptw_re" %in% names(data)) {
  # run bootstrap 
  refe_ci <- bootstrapCIb(
    iterations = 1000,
    iptw = iptw_re,
    data = data,
    model = "FE",
    cores = 6,
    core_seeds = c(4561:4566),
    # Seeds for reproducibility
    effect_type = "TNIE"
  )
  # save 
  saveRDS(refe_ci, file = "Application/Output/Temp-bootstrap-CIs/refe_ci.rds")
} else {
  print("The column 'iptw_re' does not exist in the 'data' data frame.")
}

rm(refe_ci)


### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
if ("iptw_sl" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    slre_ci <- bootstrapCIb(
      iterations = 1000, #1700,
      iptw = iptw_sl,
      data = data,
      model = "RE",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(slre_ci, file = "Application/Output/Temp-bootstrap-CIs/slre_ci.rds")
} else {
  print("The column 'iptw_sl' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", slre_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_ci$mediator_converged_count,
    " (", (slre_ci$mediator_converged_count / length(slre_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_ci$outcome_converged_count,
    " (", (slre_ci$outcome_converged_count / length(slre_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_ci$both_converged_count,
    " (", (slre_ci$both_converged_count / length(slre_ci$PNDE_est)) * 100, "%)\n")
rm(slre_ci)
# Elapsed time: 6652.806 seconds ( 111 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# FE PS Model
if ("iptw_fe" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    fere_ci <- bootstrapCIb(
      iterations = 1000, #1700, 
      iptw = iptw_fe,
      data = data,
      model = "RE",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(fere_ci, file = "Application/Output/Temp-bootstrap-CIs/fere_ci.rds")
} else {
  print("The column 'iptw_fe' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", fere_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_ci$mediator_converged_count,
    " (", (fere_ci$mediator_converged_count / length(fere_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_ci$outcome_converged_count,
    " (", (fere_ci$outcome_converged_count / length(fere_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_ci$both_converged_count,
    " (", (fere_ci$both_converged_count / length(fere_ci$PNDE_est)) * 100, "%)\n")
rm(fere_ci)
# Elapsed time: 7617.093 seconds ( 127 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# RE PS Model
if ("iptw_re" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    rere_ci <- bootstrapCIb(
      iterations = 1000, #1700, 
      iptw = iptw_re,
      data = data,
      model = "RE",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(rere_ci, file = "Application/Output/Temp-bootstrap-CIs/rere_ci.rds")
} else {
  print("The column 'iptw_re' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", rere_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_ci$mediator_converged_count,
    " (", (rere_ci$mediator_converged_count / length(rere_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_ci$outcome_converged_count,
    " (", (rere_ci$outcome_converged_count / length(rere_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_ci$both_converged_count,
    " (", (rere_ci$both_converged_count / length(rere_ci$PNDE_est)) * 100, "%)\n")
rm(rere_ci)
# Elapsed time: 10816.36 seconds ( 180 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)


### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
if ("iptw_sl" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    slre_cm_ci <- bootstrapCIb(
      iterations = 1000, #1750, 
      iptw = iptw_sl,
      data = data,
      model = "RECM",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(slre_cm_ci, file = "Application/Output/Temp-bootstrap-CIs/slre_cm_ci.rds")
} else {
  print("The column 'iptw_sl' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", slre_cm_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_cm_ci$mediator_converged_count,
    " (", (slre_cm_ci$mediator_converged_count / length(slre_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_cm_ci$outcome_converged_count,
    " (", (slre_cm_ci$outcome_converged_count / length(slre_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_cm_ci$both_converged_count,
    " (", (slre_cm_ci$both_converged_count / length(slre_cm_ci$PNDE_est)) * 100, "%)\n")
rm(slre_cm_ci)
# Elapsed time: 6615.656 seconds ( 110 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# FE PS Model
if ("iptw_fe" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    fere_cm_ci <- bootstrapCIb(
      iterations = 1000, #1750, 
      iptw = iptw_fe,
      data = data,
      model = "RECM",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(fere_cm_ci, file = "Application/Output/Temp-bootstrap-CIs/fere_cm_ci.rds")
} else {
  print("The column 'iptw_fe' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", fere_cm_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_cm_ci$mediator_converged_count,
    " (", (fere_cm_ci$mediator_converged_count / length(fere_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_cm_ci$outcome_converged_count,
    " (", (fere_cm_ci$outcome_converged_count / length(fere_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_cm_ci$both_converged_count,
    " (", (fere_cm_ci$both_converged_count / length(fere_cm_ci$PNDE_est)) * 100, "%)\n")
rm(fere_cm_ci)
# Elapsed time: 6152.866 seconds ( 103 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# RE PS Model
if ("iptw_re" %in% names(data)) {
  execution_time <- system.time({# track computation time
    # run bootstrap
    rere_cm_ci <- bootstrapCIb(
      iterations = 1000, #1750, 
      iptw = iptw_re,
      data = data,
      model = "RECM",
      cores = 6,
      core_seeds = c(4561:4566),
      # Seeds for reproducibility
      effect_type = "TNIE"
    )
  })
  # save
  saveRDS(rere_cm_ci, file = "Application/Output/Temp-bootstrap-CIs/rere_cm_ci.rds")
} else {
  print("The column 'iptw_re' does not exist in the 'data' data frame.")
}

# Print elapsed time and convergence statistics
cat("~~~~~~~~~", rere_cm_ci$analysisCond, "~~~~~~~~~")
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_cm_ci$mediator_converged_count,
    " (", (rere_cm_ci$mediator_converged_count / length(rere_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_cm_ci$outcome_converged_count,
    " (", (rere_cm_ci$outcome_converged_count / length(rere_cm_ci$PNDE_est)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_cm_ci$both_converged_count,
    " (", (rere_cm_ci$both_converged_count / length(rere_cm_ci$PNDE_est)) * 100, "%)\n")
rm(rere_cm_ci)
# Elapsed time: 5956.323 seconds ( 99 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)





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







