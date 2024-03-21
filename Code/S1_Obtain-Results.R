################################################################################
########################## QP Simulation 1 Results #############################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 03/20/2024
#
#
# Script Description: This code summarizes and reports the results for the 
#                       first simulation study (i.e., obtains performance measures).  
#                       This is stored in the relevant Results folder.
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
  tidyverse, 
  ggplot2, 
  flextable, 
  stringr,
  ggdag, 
  dagitty, 
  huxtable
)



# Set Parameters & Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = 100, 
                    clust_size = c(20, 40, 100), 
                    num_x = 3, 
                    icc = c(0.05, 0.2, 0.5)) 


# set direct & indirect effects 
treat_m <- 1.17 # trt on med   
treat_y <- 1.35 # trt on outcome
med_y <- 1.2 # med on outcome  
TNDE <- treat_y
PNIE <- treat_m * med_y

# Create directory to store reporting of results 
dir.create(path = "Output/S1_Results")
dir.create(path = paste0("Output/S1_Results/Data"))
dir.create(path = paste0("Output/S1_Results/Tables"))
dir.create(path = paste0("Output/S1_Results/Figures"))


# Import data  ------------------------------------------------------------

sim1_data <- NULL
# Loop to store in single dataframe
for(i in 1:nrow(cond)) {
  temp_data <-
    readRDS(
      paste0(
        "Output/S1_Simulation-Output/S1_Condition-", i, "-Estimates.rds"
      )
    )
  
  sim1_data <- as.data.frame(rbind(sim1_data,
                                   temp_data))
  
}

rm(temp_data)



# Clean Data --------------------------------------------------------------

# Change row names 
rownames(sim1_data) <- 1:nrow(sim1_data)
# Change to numeric & properly name effects 
sim1_data <- sim1_data %>% 
  mutate_at(.vars = c("NDE_est", "NIE_est", "ICC", "clust_size", "conditionNum", 
                      "a_path_est", "a_path_se", "b_path_est", "b_path_se", 
                      "direct_est", "direct_se"), 
            .funs = as.numeric) %>% 
  rename(TNDE_est = NDE_est, 
         PNIE_est = NIE_est)

# Check 
head(sim1_data)



# Compute Performance Measures --------------------------------------------

# Performance measures summary DF 
perf_measure_DF <- sim1_data %>% 
  group_by(ICC, clust_size, 
           conditionNum, analysisCond) %>% 
  summarize(TNDE_relBias = (mean(TNDE_est) / TNDE) - 1, 
            NIE_relBias = (mean(PNIE_est) / PNIE) - 1, 
            
            # TNDE_MSE = (mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2), 
            # PNIE_MSE = (mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2), 
            
            TNDE_RMSE = sqrt((mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2)), 
            PNIE_RMSE = sqrt((mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2))
  ) 



# Export Performance Measures & Simulation Data ---------------------------------------------

write_rds(perf_measure_DF, 
          file = paste0("Output/S1_Results/Data/", "S1_Performance-Measures.rds"))

write_rds(sim1_data, 
          file = paste0("Output/S1_Results/Data/", "S1_Simulation-Data.rds"))

