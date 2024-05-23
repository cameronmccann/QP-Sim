################################################################################
########################## QP Simulation 2 Results #############################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 03/21/2024
#
#
# Script Description: This code summarizes and reports the results for the 
#                       second simulation study (i.e., obtains performance measures).  
#                       This is stored in the relevant Results folder.
#
#
# Last Updated: 05/23/2024
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
TNDE <- treat_y
PNIE <- treat_m * med_y

# Create directory to store reporting of results 
dir.create(path = "Output/S2_Results")
dir.create(path = paste0("Output/S2_Results/Data"))
dir.create(path = paste0("Output/S2_Results/Tables"))
dir.create(path = paste0("Output/S2_Results/Figures"))



# Import data  ------------------------------------------------------------

sim2_data <- NULL

for (i in 1:nrow(cond)) {
  temp_data <-
    readRDS(paste0(
      "Output/S2_Simulation-Output/S2_Condition-",
      i,
      "-Estimates.rds"
    ))
  
  sim2_data <- as.data.frame(rbind(sim2_data,
                                   temp_data))
  
}

rm(temp_data)



# Clean Data --------------------------------------------------------------

# Change row names 
rownames(sim2_data) <- 1:nrow(sim2_data)
# Change to numeric & properly name effects 
sim2_data <- sim2_data %>% 
  mutate_at(.vars = c("NDE_est", "NIE_est", "ICC", "clust_size", "conditionNum", 
                      "a_path_est", "a_path_se", "b_path_est", "b_path_se", 
                      "direct_est", "direct_se"), 
            .funs = as.numeric) %>% 
  rename(TNDE_est = NDE_est, 
         PNIE_est = NIE_est)

# Check 
head(sim2_data)



# Compute Performance Measures --------------------------------------------

# Performance measures summary DF 
perf_measure_DF <- sim2_data %>% 
  group_by(ICC, clust_size, 
           conditionNum, analysisCond) %>% 
  summarize(TNDE_relBias = (mean(TNDE_est) / TNDE) - 1, 
            PNIE_relBias = (mean(PNIE_est) / PNIE) - 1, 
            
            b_relBias = (mean(b_path_est) / med_y) - 1,
            
            # TNDE_MSE = (mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2), 
            # PNIE_MSE = (mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2), 
            
            TNDE_RMSE = sqrt((mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2)), 
            PNIE_RMSE = sqrt((mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2))
  ) 



# Export Performance Measures & Simulation Data ---------------------------------------------

write_rds(perf_measure_DF, 
          file = paste0("Output/S2_Results/Data/", "S2_Performance-Measures.rds"))

write_rds(sim2_data, 
          file = paste0("Output/S2_Results/Data/", "S2_Simulation-Data.rds"))



# TNDE Relative Bias Table ------------------------------------------------

# Pivot wide 
relBias_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "TNDE_relBias",
                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names 
colnames(relBias_Tbl) <-
  stringr::str_remove(string = colnames(relBias_Tbl),
                      pattern = "_relBias")

# Relative Bias TNDE Table   
relBias_TNDE_Table <- relBias_Tbl %>% 
  select(analysisCond, starts_with("TNDE_")) %>% 
  separate(col = c("analysisCond"), into = c("PS Model", "Mediator/Outcome Model"), sep = "_") %>% 
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_TNDE_Table) %>% 
  set_number_format(row = 1:nrow(relBias_TNDE_Table) + 1, col = -c(1:2), value = 3)



# PNIE Relative Bias Table ------------------------------------------------

# Pivot wide
relBias_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "TNDE_relBias",
                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names
colnames(relBias_Tbl) <-
  stringr::str_remove(string = colnames(relBias_Tbl),
                      pattern = "_relBias")

# Relative Bias PNIE Table   
relBias_PNIE_Table <- relBias_Tbl %>%
  select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_PNIE_Table) %>%
  set_number_format(
    row = 1:nrow(relBias_PNIE_Table) + 1,
    col = -c(1:2),
    value = 3
  )



# TNDE RMSE Table ------------------------------------------------

# Pivot wide
rmse_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNIE_RMSE",
                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl) <-
  stringr::str_remove(string = colnames(rmse_Tbl),
                      pattern = "_RMSE")

# RMSE TNDE Table   
rmse_TNDE_Table <- rmse_Tbl %>%
  select(analysisCond, starts_with("TNDE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_TNDE_Table) %>%
  set_number_format(
    row = 1:nrow(rmse_TNDE_Table) + 1,
    col = -c(1:2),
    value = 3
  )



# PNIE RMSE Table ------------------------------------------------

# Pivot wide
rmse_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNIE_RMSE",
                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl) <-
  stringr::str_remove(string = colnames(rmse_Tbl),
                      pattern = "_RMSE")

# RMSE PNIE Table   
rmse_PNIE_Table <- rmse_Tbl %>%
  select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_PNIE_Table) %>%
  set_number_format(
    row = 1:nrow(rmse_PNIE_Table) + 1,
    col = -c(1:2),
    value = 3
  )


