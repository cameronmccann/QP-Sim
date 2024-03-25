################################################################################
########################## QP Simulation 2B Results #############################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 03/21/2024
#
#
# Script Description: This code summarizes and reports the results for the 
#                       second simulation study (i.e., obtains performance measures), 
#                       version B (where the Total Natural Indirect Effect & Pure 
#                       Natural Direct Effect are estimated). 
#                       This is stored in the relevant Results folder.
#
#
# Last Updated: 03/25/2024
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



# Set parameters & simulation conditions ----------------------------------

cond <- expand.grid(num_clust = 100, 
                    clust_size = c(20, 40, 100), 
                    num_x = 6, 
                    icc = c(0.05, 0.2, 0.5)) 


# set direct & indirect effects 
treat_m <- 1 # trt on med   
treat_y <- 1.3 # trt on outcome
med_y <- 1 # med on outcome  
treat_med_y <- 1.15 # trt-med interaction on outcome 
PNDE <- treat_y
PNIE <- treat_m * med_y + treat_med_y

# Create directory to store reporting of results 
dir.create(path = "Output/S2B_Results")
dir.create(path = paste0("Output/S2B_Results/Data"))
dir.create(path = paste0("Output/S2B_Results/Tables"))
dir.create(path = paste0("Output/S2B_Results/Figures"))



# Import data  ------------------------------------------------------------

sim2B_data <- NULL

for (i in 1:nrow(cond)) {
  temp_data <-
    readRDS(paste0(
      "Output/S2B_Simulation-Output/S2B_Condition-",
      i,
      "-Estimates.rds"
    ))
  
  sim2B_data <- as.data.frame(rbind(sim2B_data,
                                   temp_data))
  
}

rm(temp_data)



# Clean Data --------------------------------------------------------------

# Change row names 
rownames(sim2B_data) <- 1:nrow(sim2B_data)
# Change to numeric & properly name effects 
sim2B_data <- sim2B_data %>% 
  mutate_at(.vars = c("NDE_est", "NIE_est", "ICC", "clust_size", "conditionNum", 
                      "a_path_est", "a_path_se", "b_path_est", "b_path_se", 
                      "direct_est", "direct_se"), 
            .funs = as.numeric) %>% 
  rename(PNDE_est = NDE_est, 
         TNIE_est = NIE_est)

# Check 
head(sim2B_data)



# Compute Performance Measures --------------------------------------------

# Performance measures summary DF 
perf_measure_DF <- sim2B_data %>% 
  group_by(ICC, clust_size, 
           conditionNum, analysisCond) %>% 
  summarize(PNDE_relBias = (mean(PNDE_est) / PNDE) - 1, 
            TNIE_relBias = (mean(TNIE_est) / TNIE) - 1, 
            
            b_relBias = (mean(b_path_est) / med_y) - 1,
            
            # TNDE_MSE = (mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2), 
            # PNIE_MSE = (mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2), 
            
            PNDE_RMSE = sqrt((mean(PNDE_est) - PNDE)^2 + (sd(PNDE_est)^2)), 
            TNIE_RMSE = sqrt((mean(TNIE_est) - TNIE)^2 + (sd(TNIE_est)^2))
  ) 



# Export Performance Measures & Simulation Data ---------------------------------------------

write_rds(perf_measure_DF, 
          file = paste0("Output/S2B_Results/Data/", "S2B_Performance-Measures.rds"))

write_rds(sim2B_data, 
          file = paste0("Output/S2B_Results/Data/", "S2B_Simulation-Data.rds"))



# TNDE Relative Bias Table ------------------------------------------------

# Pivot wide 
relBias_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNDE_relBias",
                      "TNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNDE_relBias", "TNIE_relBias")
)

# drop common peice of column names 
colnames(relBias_Tbl) <-
  stringr::str_remove(string = colnames(relBias_Tbl),
                      pattern = "_relBias")

# Relative Bias TNDE Table   
relBias_PNDE_Table <- relBias_Tbl %>% 
  select(analysisCond, starts_with("PNDE_")) %>% 
  separate(col = c("analysisCond"), into = c("PS Model", "Mediator/Outcome Model"), sep = "_") %>% 
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_PNDE_Table) %>% 
  set_number_format(row = 1:nrow(relBias_PNDE_Table) + 1, col = -c(1:2), value = 3)



# PNIE Relative Bias Table ------------------------------------------------

# Pivot wide
relBias_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNDE_relBias",
                      "TNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNDE_relBias", "TNIE_relBias")
)

# drop common peice of column names
colnames(relBias_Tbl) <-
  stringr::str_remove(string = colnames(relBias_Tbl),
                      pattern = "_relBias")

# Relative Bias PNIE Table   
relBias_TNIE_Table <- relBias_Tbl %>%
  select(analysisCond, starts_with("TNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_TNIE_Table) %>%
  set_number_format(
    row = 1:nrow(relBias_TNIE_Table) + 1,
    col = -c(1:2),
    value = 3
  )



