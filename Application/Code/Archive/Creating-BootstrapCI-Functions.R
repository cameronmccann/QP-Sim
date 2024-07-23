################################################################################
############## QP Empirical Application - Create Bootstrap function ############
################################################################################

############################ Script Description ################################
#
# Author: Cameron
# 
# Date Created: 07/23/24
#
#
# Script Description:
#           This is a draft document to help create the bootstrap function for the empirical example. 
#
# Last Updated: 07/23/2024 
#
#
# Notes:
#   To-Do:
#     # Try to parallelize the RE function 
# 
#   Next: 
#     
#   Done: 
#     # make 2 bootstrap functions: 1) SL & FE med/out models & 2) RE med/out models
#     # make each function return direct & indirect CIs
# 
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






# Create RERE bootstrap CI function  --------------------------------------

## Initial Draft of bootstrap function for RERE ----------------------------


# RERE
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
    WeMix::mix(formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
                 white_w1 + black_w1 + 
                 parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  outcome_rere <- tryCatch({
    WeMix::mix(formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
                 white_w1 + black_w1 + 
                 parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
               data = data_boot, 
               weights = c("iptw_re", "L2weight"))
  }, error = function(e) NULL)
  
  # Check convergence and calculate indirect effect if both models converged
  mediator_converged[i] <- !is.null(mediator_rere)
  outcome_converged[i] <- !is.null(outcome_rere)
  
  if (mediator_converged[i] && outcome_converged[i]) {
    indirect_effects_boot_rere[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
      summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
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





## Initial AI assisted output  ---------------------------------------------


bootstrap_ci <- function(seed = 456, iterations = 50, iptw, data) {
  # Initialize vectors to store the bootstrap indirect effects and convergence statuses
  indirect_effects_boot_rere <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    # Add column of just 1s for level-2 weight
    data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
    
    # Fit the RE models with RE ps for the bootstrap sample
    mediator_rere <- tryCatch({
      WeMix::mix(
        formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    outcome_rere <- tryCatch({
      WeMix::mix(
        formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    # Check convergence and calculate indirect effect if both models converged
    mediator_converged[i] <- !is.null(mediator_rere)
    outcome_converged[i] <- !is.null(outcome_rere)
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects_boot_rere[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
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
  
  # Return the results as a list
  list(
    boot_ci_rere = boot_ci_rere,
    indirect_effects_boot_rere = indirect_effects_boot_rere,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci(seed = 456, iterations = 50, iptw = iptw_re, data = data)
print(result$boot_ci_rere)

result$mediator_converged_count
result$outcome_converged_count
result$both_converged_count

head(result$indirect_effects_boot_rere)

# This works lets attempt to add direct effect to it 




## Add direct effect CI to RE bootstrap function  --------------------------

bootstrap_ci_re <- function(seed = 456, iterations = 50, iptw, data) {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    # Add column of just 1s for level-2 weight
    data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
    
    # Fit the RE models with RE ps for the bootstrap sample
    mediator_rere <- tryCatch({
      WeMix::mix(
        formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    outcome_rere <- tryCatch({
      WeMix::mix(
        formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator_rere)
    outcome_converged[i] <- !is.null(outcome_rere)
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome_rere)$coef["sportPartic_w1", "Estimate"]
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci_re(seed = 456, iterations = 10, iptw = iptw_re, data = data)
# print(result$boot_ci_rere)
# print(result$direct_ci_rere)

print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count
result$outcome_converged_count



## Parallelize the RE bootstrap function  ----------------------------------




cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

par_time <- system.time(
  
  cond_results_DF <- foreach::foreach(
    i = 1:reps,
    .combine = rbind, 
    .export = c("bootstrap_ci_re")
  ) %dopar% {
    
  }
)

parallel::stopCluster(cl)




# Code from conducting simulation 
for (condition in 1:nrow(cond)) {
  
  # set condition number
  cond_num <- condition
  
  # Make/Register cores
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl) # 7 cores (personal); 19 (campus/lab)
  
  
  # Conduct Simulation
  par_time <- system.time(
    # Track computation time
    
    cond_Results_DF <- foreach::foreach(
      i = 1:reps,
      .combine = rbind,
      .export = c("genOneData_Sim2", "AnalysisFunc_Sim2")
    ) %dopar% {
      set.seed(paste0(135, i))
      
      # Generate data set
      data <- genOneData_Sim2(
        num_clust = cond[cond_num, "num_clust"],
        clust_size = cond[cond_num, "clust_size"],
        num_x = cond[cond_num, "num_x"],
        icct = cond[cond_num, "icc"],
        iccm = cond[cond_num, "icc"],
        iccy = cond[cond_num, "icc"]
      )
      
      # Analyze data set
      temp_DF <- NULL
      full_DF <- NULL
      
      for (PSmod in c("FE", "RE")) { # PS Models
        for (Outmod in c("SL", "FE", "RE", "RE-Mean")) { # Med/Outcome Models
          
          temp_DF <- AnalysisFunc_Sim2(
            PSmodel = PSmod,
            Medmodel = Outmod,
            Outcomemodel = Outmod,
            data = data,
            condition = cond,
            condition_num = cond_num
          )
          
          # Combine results for condition
          full_DF <- as.data.frame(rbind(full_DF, temp_DF))
          
          # print(paste0("results for ", PSmod, "_", Outmod, " Done!"))
        }
      }
      
      results <- as.data.frame(cbind(
        seed = rep(paste0(135, i), nrow(full_DF)),
        rep = rep(i, nrow(full_DF)),
        full_DF, 
        true_ps_10pctle = rep(quantile(data$ps_true, probs = c(0.1)), nrow(full_DF)), 
        true_ps_15pctle = rep(quantile(data$ps_true, probs = c(0.15)), nrow(full_DF)), 
        true_ps_50pctle = rep(quantile(data$ps_true, probs = c(0.5)), nrow(full_DF)), 
        true_ps_85pctle = rep(quantile(data$ps_true, probs = c(0.85)), nrow(full_DF)), 
        true_ps_90pctle = rep(quantile(data$ps_true, probs = c(0.9)), nrow(full_DF))
      ))
    }
  )
  
  parallel::stopCluster(cl)
  
  # Save conditions results/simulation output 
  saveRDS(
    cond_Results_DF,
    file = paste0(
      "Output/S2_Simulation-Output_", 
      date, 
      "/S2_Condition-", 
      condition, 
      "-Estimates.rds"
    )
  )
  
  # cond_Results_DF
  
  # Print message
  print(paste0("Condition ", condition, " Done!"))
  
  
  # Log computation time 
  OverallPar_time <- rbind(OverallPar_time,
                           cbind(t(as.data.frame(par_time)), cond_num))
  
}






## 1st AI assistance  ------------------------------------------------------


library(parallel)
library(doParallel)
library(foreach)
library(WeMix)
library(txtProgressBar)

# Define the bootstrap_ci_re function
bootstrap_ci_re <- function(seed = 456, iterations = 50, iptw, data) {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    # Add column of just 1s for level-2 weight
    data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
    
    # Fit the RE models with RE ps for the bootstrap sample
    mediator_rere <- tryCatch({
      WeMix::mix(
        formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    outcome_rere <- tryCatch({
      WeMix::mix(
        formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator_rere)
    outcome_converged[i] <- !is.null(outcome_rere)
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome_rere)$coef["sportPartic_w1", "Estimate"]
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Parallelization setup
cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

# Parallel bootstrap sampling
par_time <- system.time(
  cond_results_DF <- foreach::foreach(
    i = 1:iterations,
    .combine = rbind,
    .export = c("bootstrap_ci_re", "data", "iptw_str")
  ) %dopar% {
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
    
    mediator_rere <- tryCatch({
      WeMix::mix(
        formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    outcome_rere <- tryCatch({
      WeMix::mix(
        formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)
    
    mediator_converged <- !is.null(mediator_rere)
    outcome_converged <- !is.null(outcome_rere)
    
    if (mediator_converged && outcome_converged) {
      indirect_effects <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects <- summary(outcome_rere)$coef["sportPartic_w1", "Estimate"]
    } else {
      indirect_effects <- NA
      direct_effects <- NA
    }
    
    c(indirect_effects, direct_effects, mediator_converged, outcome_converged)
  }
)

parallel::stopCluster(cl)

# Print the results
cond_results_DF








#####
# Load necessary libraries
library(parallel)
library(doParallel)
library(foreach)

# Number of iterations
iterations <- 100

# Number of random numbers in each iteration
n <- 1000

# Set up the parallel backend
cl <- makeCluster(detectCores() - 1)
cl <- makeCluster(3)
registerDoParallel(cl)

# Parallel computation using foreach
result <- foreach(i = 1:iterations, .combine = c) %dopar% {
  # Generate random numbers and calculate their mean
  mean(rnorm(n))
}

# Stop the cluster
stopCluster(cl)

# Print the result
print(result)





















# Create SL & FE bootstrap CI function ------------------------------------





## Initial start  ----------------------------------------------------------

bootstrap_ci <- function(seed = 456, iterations = 50, iptw, data) {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  # iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    # Fit the SL models for the bootstrap sample
    mediator <- tryCatch({
      glm(
        formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
        data = data_boot,
        weights = iptw
      )
    }, error = function(e) NULL)
    
    outcome <- tryCatch({
      glm(
        formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
        data = data_boot,
        weights = iptw
      )
    }, error = function(e) NULL)
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator)
    outcome_converged[i] <- !is.null(outcome)
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome)$coef["sportPartic_w1", "Estimate"]
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci(seed = 456, iterations = 10, iptw = iptw_sl, data = data)
# print(result$boot_ci_rere)
# print(result$direct_ci_rere)

print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count
result$outcome_converged_count




## Debugging function  -----------------------------------------------------

bootstrap_ci <- function(seed = 456, iterations = 50, iptw, data) {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    cat("\nIteration:", i, "\n")
    
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    cat("Resampled clusters:", unique(cluster_boot), "\n")
    cat("Resampled data size:", nrow(data_boot), "\n")
    
    # Fit the SL models for the bootstrap sample
    mediator <- tryCatch({
      glm(
        formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
        data = data_boot,
        weights = iptw
      )
    }, error = function(e) {
      cat("Mediator model failed: ", e$message, "\n")
      NULL
    })
    
    outcome <- tryCatch({
      glm(
        formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
        data = data_boot,
        weights = iptw
      )
    }, error = function(e) {
      cat("Outcome model failed: ", e$message, "\n")
      NULL
    })
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator)
    outcome_converged[i] <- !is.null(outcome)
    cat("Mediator model converged:", mediator_converged[i], "\n")
    cat("Outcome model converged:", outcome_converged[i], "\n")
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome)$coef["sportPartic_w1", "Estimate"]
      cat("Indirect effect estimate:", indirect_effects[i], "\n")
      cat("Direct effect estimate:", direct_effects[i], "\n")
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
      cat("One or both models did not converge. Setting effects to NA.\n")
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci(seed = 456, iterations = 10, iptw = iptw_sl, data = data)
# print(result$indirect_ci)
# print(result$direct_ci)


######################


bootstrap_ci <- function(seed = 456, iterations = 50, iptw, data, model = "SL") {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # cat("\nIteration:", i, "\n")
    
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    # cat("Resampled clusters:", unique(cluster_boot), "\n")
    # cat("Resampled data size:", nrow(data_boot), "\n")
    
    if (model == "SL") {
      # Fit the SL models for the bootstrap sample
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Mediator model failed: ", e$message, "\n")
        NULL
      })
      
      outcome <- tryCatch({
        glm(
          formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Outcome model failed: ", e$message, "\n")
        NULL
      })
    } else if (model == "FE") {
      # Fit the FE models for the bootstrap sample
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2), 
          data = data_boot, 
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Mediator FE model failed: ", e$message, "\n")
        NULL
      })
      
      outcome <- tryCatch({
        glm(
          formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2),
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Outcome FE model failed: ", e$message, "\n")
        NULL
      })
    } else {
      stop("Invalid model type. Choose 'SL' or 'FE'.")
    }
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator)
    outcome_converged[i] <- !is.null(outcome)
    # cat("Mediator model converged:", mediator_converged[i], "\n")
    # cat("Outcome model converged:", outcome_converged[i], "\n")
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome)$coef["sportPartic_w1", "Estimate"]
      # cat("Indirect effect estimate:", indirect_effects[i], "\n")
      # cat("Direct effect estimate:", direct_effects[i], "\n")
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
      # cat("One or both models did not converge. Setting effects to NA.\n")
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_re", data = data, model = "SL")
print(result$indirect_ci)
print(result$direct_ci)


result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_re", data = data, model = "FE")
print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count





# Add FE as an option to function  ----------------------------------------

bootstrap_ci <- function(seed = 456, iterations = 50, iptw, data, model = "SL") {
  # Initialize vectors to store the bootstrap indirect and direct effects and convergence statuses
  indirect_effects <- numeric(iterations)
  direct_effects <- numeric(iterations)
  mediator_converged <- logical(iterations)
  outcome_converged <- logical(iterations)
  
  # Convert iptw to string
  iptw_str <- deparse(substitute(iptw))
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3, width = 50, char = "=")
  
  # Perform bootstrap resampling
  for (i in 1:iterations) {
    # cat("\nIteration:", i, "\n")
    
    # Resample with replacement at the cluster level
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    # cat("Resampled clusters:", unique(cluster_boot), "\n")
    # cat("Resampled data size:", nrow(data_boot), "\n")
    
    if (model == "SL") {
      # Fit the SL models for the bootstrap sample
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Mediator model failed: ", e$message, "\n")
        NULL
      })
      
      outcome <- tryCatch({
        glm(
          formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Outcome model failed: ", e$message, "\n")
        NULL
      })
    } else if (model == "FE") {
      # Fit the FE models for the bootstrap sample
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2), 
          data = data_boot, 
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Mediator FE model failed: ", e$message, "\n")
        NULL
      })
      
      outcome <- tryCatch({
        glm(
          formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
            parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2),
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) {
        # cat("Outcome FE model failed: ", e$message, "\n")
        NULL
      })
    } else {
      stop("Invalid model type. Choose 'SL' or 'FE'.")
    }
    
    # Check convergence and calculate effects if both models converged
    mediator_converged[i] <- !is.null(mediator)
    outcome_converged[i] <- !is.null(outcome)
    # cat("Mediator model converged:", mediator_converged[i], "\n")
    # cat("Outcome model converged:", outcome_converged[i], "\n")
    
    if (mediator_converged[i] && outcome_converged[i]) {
      indirect_effects[i] <- summary(mediator)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome)$coef["selfEst_w3_sc", "Estimate"]
      direct_effects[i] <- summary(outcome)$coef["sportPartic_w1", "Estimate"]
      # cat("Indirect effect estimate:", indirect_effects[i], "\n")
      # cat("Direct effect estimate:", direct_effects[i], "\n")
    } else {
      indirect_effects[i] <- NA
      direct_effects[i] <- NA
      # cat("One or both models did not converge. Setting effects to NA.\n")
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar & print completion message
  close(pb)
  cat("\nBootstrap resampling completed!\n")
  
  # Calculate the percentile bootstrap CI
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return the results as a list
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

# Example usage
result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_sl", data = data, model = "SL")
print(result$indirect_ci)
print(result$direct_ci)


result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_re", data = data, model = "FE")
print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count


result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_fe", data = data, model = "FE")
print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count
result$outcome_converged_count


result <- bootstrap_ci(seed = 456, iterations = 10, iptw = "iptw_sl", data = data, model = "FE")
print(result$indirect_ci)
print(result$direct_ci)
result$mediator_converged_count
result$outcome_converged_count
















