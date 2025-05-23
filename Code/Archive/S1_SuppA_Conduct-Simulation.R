################################################################################
######################## QP Simulation 1 Supplemental A ########################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 2025-02-08
#
#
# Script Description: This code runs the supplemental simulation for the 
#                       first simulation study (i.e., generates 
#                       & analyzes data) and stores the estimates for each 
#                       iteration in the relevant Simulation-Output folder. 
#
#
# Last Updated: 2025-02-28
#
#
# Notes:
#   To-Do
#     # 
#
#   Done: 
# 
# 
################################################################################


# ---------------------------- Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  doParallel, 
  foreach,
  parallel, 
  WeMix
)

# Load Functions 
# source("Functions/AnalysisFunc_Sim1.R")
# source("Functions/AnalysisFunc_Sim1b.R")
source("Functions/AnalysisFunc_Sim1c.R")
source("Functions/genOneData_Sim1.R")

# ---------------------------- Simulation Conditions  --------------------------------------------------

cond <- expand.grid(num_clust = c(70, 100),
                    clust_size = c(20, 40, 100),
                    num_x = 3,
                    icc = c(0.05, 0.2, 0.5))
# cond <- cond[1:5, ]

# ---------------------------- Set Parameters ----------------------------------------------------------

OverallPar_time <- NULL           # To log computation times
reps <- 1000                      # Total number of replications per condition
dir.create(path = "Output/S1_SuppA_Simulation-Output")
dir.create(path = "Output/S1_SuppA_Simulation-Output/2025-02-28_1000-reps", showWarnings = FALSE)
dir.create(path = "Output/S1_SuppA_Simulation-Output/2025-02-28_1000-reps/interim", showWarnings = FALSE)
path <- "Output/S1_SuppA_Simulation-Output/2025-02-28_1000-reps"

# ---------------------------- Supp Simulation 1 ------------------------------------------------------------

for (condition in 1:nrow(cond)) {
  
  cond_num <- condition
  
  # Register parallel backend (using all but one core)
  doParallel::registerDoParallel(parallel::detectCores() - 1)
  
  # This list will hold the results for each 100-replication block.
  cond_results_list <- list()
  total_condition_time <- 0
  
  # Process replications in blocks of 100 (adjust the block size as needed)
  for (block in seq(1, reps, by = 100)) {
    block_end <- min(block + 99, reps)
    
    # Run the replications in this block in parallel
    block_time <- system.time({
      block_results <- foreach::foreach(
        i = block:block_end,
        .combine = rbind,
        .export = c("genOneData_Sim1", "AnalysisFunc_Sim1c")
      ) %dopar% {
        # Set a seed (here using a combination of 135 and i; adjust as needed)
        set.seed(as.numeric(paste0(135, i)))
        
        # Generate data for the current replication
        data <- genOneData_Sim1(
          num_clust = cond[cond_num, "num_clust"],
          clust_size = cond[cond_num, "clust_size"],
          num_x = cond[cond_num, "num_x"],
          icct = cond[cond_num, "icc"],
          iccm = cond[cond_num, "icc"],
          iccy = cond[cond_num, "icc"]
        )
        
        # Run all models for this replication
        full_DF <- NULL
        for (PSmod in c("SL", "FE", "RE", "RE-Mean")) {      # Propensity score models
          for (Medmod in c("SL", "FE", "RE", "RE-Mean")) {
            for (Outmod in c("SL", "FE", "RE", "RE-Mean")) {   # Outcome/mediation models
              
              # create copy to avoid issues with creating new columns for some models
              data_copy <- data
              
              temp_DF <- AnalysisFunc_Sim1c(
                PSmodel = PSmod,
                Medmodel = Medmod,
                Outcomemodel = Outmod,
                data = data_copy, #data,
                condition = cond,
                condition_num = cond_num
              )
              full_DF <- rbind(full_DF, temp_DF)
            }
          }
        }

        # Add extra info (seed, replication number, quantiles of ps_true) to the results
        results <- as.data.frame(cbind(
          seed = paste0(135, i),
          rep = i,
          full_DF,
          true_ps_10pctle = rep(quantile(data$ps_true, probs = 0.1), nrow(full_DF)), 
          true_ps_15pctle = rep(quantile(data$ps_true, probs = 0.15), nrow(full_DF)), 
          true_ps_50pctle = rep(quantile(data$ps_true, probs = 0.5), nrow(full_DF)), 
          true_ps_85pctle = rep(quantile(data$ps_true, probs = 0.85), nrow(full_DF)), 
          true_ps_90pctle = rep(quantile(data$ps_true, probs = 0.9), nrow(full_DF))
        ))
        results
      }
    })
    
    total_condition_time <- total_condition_time + block_time["elapsed"]
    
    # Convert block time into minutes and seconds
    mins <- floor(block_time["elapsed"] / 60)
    secs <- round(block_time["elapsed"] %% 60, 2)
    block_time_display <- paste(mins, "mins", secs, "secs", sep = " ")
    
    # Construct a filename that includes:
    # - The condition number
    # - The condition parameters (clust_size, icc, num_clust, num_x)
    # - The replication range for this block (e.g., reps-1-100)
    file_name <- paste0(
      path, "/interim/", "/S1_Supp_Condition-", condition,
      "-Estimates_clust_size-", cond[condition, "clust_size"],
      "_icc-", cond[condition, "icc"],
      "_num_clust-", cond[condition, "num_clust"],
      "_num_x-", cond[condition, "num_x"],
      "_reps-", block, "-", block_end, ".rds"
    )
    
    # Save the block’s results to file
    saveRDS(block_results, file = file_name)
    message("Condition ", condition, ": replications ", block, " to ", block_end, " saved. (Block time: ", block_time_display, ")")
    
    # Store this block’s results in the list (if you want to combine later)
    cond_results_list[[length(cond_results_list) + 1]] <- block_results
  }
  
  # Optionally combine all block results for the current condition into one data frame
  cond_Results_DF <- do.call(rbind, cond_results_list)
  
  # Save an overall file for the condition (if desired)
  overall_file_name <- paste0(
    path, "/S1_SuppA_Condition-", condition,
    "-Overall_Estimates_clust_size-", cond[condition, "clust_size"],
    "_icc-", cond[condition, "icc"],
    "_num_clust-", cond[condition, "num_clust"],
    "_num_x-", cond[condition, "num_x"],
    "_reps-1-", reps, ".rds"
  )
  saveRDS(cond_Results_DF, file = overall_file_name)
  
  # Log the total time for the condition (with condition information)
  OverallPar_time <- rbind(OverallPar_time,
                           data.frame(condition = condition,
                                      condition_params = paste0("clust_size-", cond[condition, "clust_size"], 
                                                                "_icc-", cond[condition, "icc"], 
                                                                "_num_clust-", cond[condition, "num_clust"],
                                                                "_num_x-", cond[condition, "num_x"]),
                                      elapsed = total_condition_time))
  
  print(paste0("Condition ", condition, " Done! (Progress: ", condition, "/", nrow(cond), " = ", 
               round((condition/nrow(cond))*100), "% Complete)"))
  
  if (condition == nrow(cond)) {
    print("~~~~~ Simulation Complete ~~~~~")
  }
}

# Save the overall computation time log (with minutes computed)
OverallPar_time <- as.data.frame(OverallPar_time)
OverallPar_time$mins <- as.numeric(OverallPar_time$elapsed) / 60

# Check if the file already exists
output_file <- paste0(path, "/S1_SuppA_Computation-Time.rds")
if (file.exists(output_file)) {
  existing_data <- readRDS(output_file)
  OverallPar_time <- rbind(existing_data, OverallPar_time)
}

# Save the updated computation time log
saveRDS(OverallPar_time, file = output_file)




