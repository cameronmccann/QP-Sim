################################################################################
############################### QP Simulation 2 ################################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 03/17/2024
#
#
# Script Description: This code runs the second simulation study (i.e., generates 
#                       & analyzes data) and stores the estimates for each 
#                       iteration in the relevant Simulation-Output folder. 
#
#
# Last Updated: 03/21/2024 
#
#
# Notes:
#   To-Do
#
# 
#   Done: 
# 
# 
################################################################################


# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  doParallel, 
  foreach,
  parallel
)


# Load Functions 
source("Functions/AnalysisFunc_Sim2.R")
source("Functions/genOneData_Sim2.R")


# Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = 100,
                    clust_size = c(20, 40, 100),
                    num_x = 6,
                    icc = c(0.05, 0.2, 0.5))



# Set Parameters ----------------------------------------------------------

## Initialize DF to store results 
OverallPar_time <- NULL

## Set number of replications/repetitions 
reps <- 1000 

## Create directory to store results 
dir.create(path = "Output/S2_Simulation-Output")
path <- "Output/S2_Simulation-Output"



# Simulation 2 ------------------------------------------------------------

for (condition in 1:nrow(cond)) { 
  
  # set condition number
  cond_num <- condition
  
  # Make/Register cores
  doParallel::registerDoParallel(parallel::detectCores() - 1)
  
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
      
      for (PSmod in c("SL", "FE", "RE")) { # PS Models
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
          
        }
      }
      
      # Add extra info to DF 
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
  
  # Save conditions results/simulation output 
  saveRDS(
    cond_Results_DF,
    file = paste0(
      path, 
      "/S2_Condition-", 
      condition, 
      "-Estimates.rds"
    )
  )

  # Print message
  print(paste0("Condition ", condition, " Done! ", 
               "(Progress: ", condition, "/", nrow(cond), " = ", 
               round((condition/nrow(cond))*100), "% Complete)"))
  
  if(condition == nrow(cond)) {
    print("~~~~~ Simulation Complete ~~~~~")
  }
  
  # Log computation time 
  OverallPar_time <- rbind(OverallPar_time,
                           cbind(t(as.data.frame(par_time)), cond_num))
  
}


# Add mins to computation time log & save DF 
OverallPar_time <- as.data.frame(OverallPar_time)
OverallPar_time <- cbind(OverallPar_time, 
                         mins = OverallPar_time[, "elapsed"] / 60)

saveRDS(OverallPar_time,
        file = paste0(path, "/S2_Computation-Time.rds"))



