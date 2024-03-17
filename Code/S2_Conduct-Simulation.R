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
# Last Updated: 03/17/2024 
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
  parallel, 
  BRRR
)


# Load Functions 
source("Functions/AnalysisFunc_Sim2.R")
source("Functions/genOneData_Sim2.R")


# Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = 100,
                    clust_size = c(20, 40, 100),
                    num_x = 6,
                    icc = c(0.05, 0.2, 0.5))


# Analysis models 
# analysisCond <- expand.grid(
#   "PS model" = c("SL", "FE", "RE"),
#   "Mediator & Outcome model" = c("SL", "FE", "RE")
# )



# Set Parameters ----------------------------------------------------------

## Initialize DF to store results 
# OverallCondResultsDF <- NULL
# OverallPerfMeasuresDF <- NULL
OverallPar_time <- NULL

## Set number of replications/repetitions 
reps <- 200 #400 #10 
## Set Date 
date <- "2024-02-17"#Sys.Date()
## Create directory to store results 
dir.create(path = paste0("Output/S2_Simulation-Output_", date))


# Simulation 2 ------------------------------------------------------------

for (condition in 1:nrow(cond)) { #nrow(cond[1, ])) { 
  
  # set condition number
  cond_num <- condition
  
  # Make/Register cores
  # cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # doParallel::registerDoParallel(cl) # 7 cores (personal); 19 (campus/lab)
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
      
      # print(paste0("condition ", condition, ": ", (i/reps)*100, "%"))
      # print(paste0("condition ", condition, ": ", i))
      # print(paste0("Replication ", i, " Done!")) #new code 
    }
  )
  
  # parallel::stopCluster(cl)
  
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
  print(" ")
  print(" ")
  
  
  # Log computation time 
  OverallPar_time <- rbind(OverallPar_time,
                           cbind(t(as.data.frame(par_time)), cond_num))
  
}


# Indicate simulation is finished 
# BRRR::skrrrahh("bigshaq1")
BRRR::skrrrahh("gucci2")


# Add mins to time DF & save DF 
OverallPar_time <- as.data.frame(OverallPar_time)
OverallPar_time <- cbind(OverallPar_time, 
                         mins = OverallPar_time[, "elapsed"] / 60)

saveRDS(
  OverallPar_time, 
  file = paste0(
    "Output/S2_Simulation-Output_", 
    date, 
    "/S2_Computation-Time.rds"
  )
)




