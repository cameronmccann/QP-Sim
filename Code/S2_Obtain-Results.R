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

# drop common piece of column names 
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


# TNDE Relative Bias Visual -----------------------------------------------

# set values 
treat_m <- 1 # trt on med   
treat_y <- 1.3 # trt on outcome
med_y <- 1 # med on outcome  
TNDE <- treat_y
PNIE <- treat_m * med_y

# visual settings
gglayer_theme <- list(theme_bw(),
                      scale_fill_manual(values = c("#BF5700", #Fixed-effect
                                                   "#A6CD57", #Random-effect 
                                                   "#333F48")), #Single-level 
                      #"#9CADB7" <-- light gray 
                      # Used following website with university colors: https://projects.susielu.com/viz-palette?
                      theme(text = element_text(family = "Times New Roman", size = 12), 
                            axis.title = element_text(size = 12),  # Adjust axis title size
                            axis.text = element_text(size = 10),  # Adjust axis text size
                            legend.title = element_text(size = 12),  # Legend title size
                            legend.text = element_text(size = 10),  # Legend text size
                            strip.text = element_text(size = 12),  # Facet labels
                            line = element_line(linewidth = 0.5),  # APA recommends thin lines
                            legend.position = "top"  
                      )) 

# Labels 
gglayer_labs <- list(
  labs(
    x = "\n Absolute Relative Bias",
    y = "Mediator and Outcome Model \n",
    color = "PS Model",
    linetype = "PS Model",
    shape = "PS Model"
  ),
  guides(y.sec = guide_none("Cluster Size"),
         x.sec = guide_none("Residual ICC"))
)

pdf("Output/S2_Results/Figures/TNDE-Relative-Bias-Boxplot.pdf")
readRDS(file = "Output/S2_Results/Data/S2_Simulation-Data.rds") |>
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           "ERROR")))) %>% 
  
  # visual 
  ggplot(aes(x = abs((TNDE_est / TNDE) - 1), y = outModel, fill = `PS Model`)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.4) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.25, 
               outlier.size = 0.5, 
               outlier.alpha = 0.5) +
  facet_grid(ICC ~ clust_size) +
  gglayer_labs +
  gglayer_theme 

dev.off()
# Save visual 
ggsave(filename = "Output/S2_Results/Figures/TNDE-Relative-Bias-Boxplot.png", plot = last_plot())


