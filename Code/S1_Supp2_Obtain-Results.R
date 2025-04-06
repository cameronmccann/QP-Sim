################################################################################
#################### QP Simulation 1 Supplemental 2 Results ####################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 2025-03-09
#
#
# Script Description: This code summarizes and reports the results for 
#                       Simulation Study 1 Supplemental 2.  
#                       This is stored in the relevant Results folder.
#
#
# Last Updated: 2025-03-09
#
#
# Notes:
#   To-Do
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
  huxtable, 
  dplyr
)

# Set Parameters & Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = c(100),
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
dir.create(path = "Output/S1_Supp2_Results")
path <- "Output/S1_Supp2_Results/2025-03-02_1000-reps"
dir.create(path = path)
dir.create(path = paste0(path, "/Data"))
dir.create(path = paste0(path, "/Tables"))
dir.create(path = paste0(path, "/Figures"))
retrieval_path <- "Output/S1_Supp2_Simulation-Output/2025-03-02_1000-reps"


# Import data  ------------------------------------------------------------

# List all files matching the pattern.
file_list <- list.files(path = retrieval_path,
                        pattern = "^S1_Supp2_Condition-[0-9]+-Overall_Estimates_.*\\.rds$",
                        full.names = TRUE)
# Read each file and combine into one data frame
sim1_data <- do.call(rbind, lapply(file_list, readRDS))


# Clean Data --------------------------------------------------------------

# Change row names 
rownames(sim1_data) <- 1:nrow(sim1_data)
# Change to numeric & properly name effects 
sim1_data <- sim1_data %>% 
  mutate_at(.vars = c("NDE_est", "NIE_est", 
                      "ICC", "clust_size", "num_clust", "conditionNum", 
                      "a_path_est", "a_path_se", "b_path_est", "b_path_se", 
                      "direct_est", "direct_se"), 
            .funs = as.numeric) %>% 
  rename(TNDE_est = NDE_est, # change TNDE to PNDE with TNIE (default)
         PNIE_est = NIE_est)

# Check 
head(sim1_data)

# drop ".2.5%" & ".97.5%" from column names
names(sim1_data) <- gsub("\\.(2\\.5%|97\\.5%)", "", names(sim1_data))


# Compute Performance Measures --------------------------------------------

# Performance measures summary DF 
perf_measure_DF <- sim1_data |> 
  mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
         if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE)) |> 
  group_by(ICC, clust_size, num_clust, 
           conditionNum, analysisCond) |> 
  summarize(TNDE_relBias = (mean(TNDE_est) / TNDE) - 1, 
            PNIE_relBias = (mean(PNIE_est) / PNIE) - 1, 
            
            # TNDE_MSE = (mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2), 
            # PNIE_MSE = (mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2), 
            
            TNDE_RMSE = sqrt((mean(TNDE_est) - TNDE)^2 + (sd(TNDE_est)^2)), 
            PNIE_RMSE = sqrt((mean(PNIE_est) - PNIE)^2 + (sd(PNIE_est)^2)), 
            # CI Coverage Rate
            coverage_TNDE = mean(if_cover_TNDE), 
            coverage_PNIE = mean(if_cover_PNIE)
  ) 


# Export Performance Measures & Simulation Data ---------------------------------------------

write_rds(perf_measure_DF, 
          file = paste0(path, "/Data/", "S1_Supp2_Performance-Measures.rds"))

write_rds(sim1_data, 
          file = paste0(path, "/Data/", "S1_Supp2_Simulation-Data.rds"))



# TNDE Relative Bias Table ------------------------------------------------

# Pivot wide 
relBias_Tbl100 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 100, c("ICC",
                                                      "clust_size",
                                                      "analysisCond",
                                                      "TNDE_relBias",
                                                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names 
colnames(relBias_Tbl100) <-
  stringr::str_remove(string = colnames(relBias_Tbl100),
                      pattern = "_relBias")

# Relative Bias TNDE Table   
relBias_TNDE_Table100 <- relBias_Tbl100 %>% 
  dplyr::select(analysisCond, starts_with("TNDE_")) %>% 
  separate(
    col = c("analysisCond"),
    into = c("PS Model","Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# separate(col = c("analysisCond"), into = c("PS Model", "Mediator Model", "Outcome Model"), sep = "_") %>% 
# arrange("Outcome Model")

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_TNDE_Table100) %>% 
  set_number_format(row = 1:nrow(relBias_TNDE_Table100) + 1, col = -c(1:2), value = 3)

# Save csv 
write.csv(relBias_TNDE_Table100,
          file = paste0(path, "/Tables/", "TNDE-Relative-Bias_num_clust-100.csv"), row.names = FALSE)


# PNIE Relative Bias Table ------------------------------------------------

# Pivot wide
relBias_Tbl100 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 100, c("ICC",
                                                      "clust_size",
                                                      "analysisCond",
                                                      "TNDE_relBias",
                                                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names
colnames(relBias_Tbl100) <-
  stringr::str_remove(string = colnames(relBias_Tbl100),
                      pattern = "_relBias")

# Relative Bias PNIE Table   
relBias_PNIE_Table100 <- relBias_Tbl100 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_PNIE_Table100) %>%
  set_number_format(
    row = 1:nrow(relBias_PNIE_Table100) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(relBias_PNIE_Table100,
          file = paste0(path, "/Tables/", "PNIE-Relative-Bias_num_clust-100.csv"), row.names = FALSE)


# TNDE RMSE Table ------------------------------------------------

# Pivot wide
rmse_Tbl100 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 100, c("ICC",
                                                      "clust_size",
                                                      "analysisCond",
                                                      "PNIE_RMSE",
                                                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl100) <-
  stringr::str_remove(string = colnames(rmse_Tbl100),
                      pattern = "_RMSE")

# RMSE TNDE Table   
rmse_TNDE_Table100 <- rmse_Tbl100 %>%
  dplyr::select(analysisCond, starts_with("TNDE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_TNDE_Table100) %>%
  set_number_format(
    row = 1:nrow(rmse_TNDE_Table100) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(rmse_TNDE_Table100,
          file = paste0(path, "/Tables/", "TNDE-RMSE_num_clust-100.csv"), row.names = FALSE)


# PNIE RMSE Table ------------------------------------------------

# Pivot wide
rmse_Tbl100 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 100, c("ICC",
                                                      "clust_size",
                                                      "analysisCond",
                                                      "PNIE_RMSE",
                                                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl100) <-
  stringr::str_remove(string = colnames(rmse_Tbl100),
                      pattern = "_RMSE")

# RMSE PNIE Table   
rmse_PNIE_Table100 <- rmse_Tbl100 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_PNIE_Table100) %>%
  set_number_format(
    row = 1:nrow(rmse_PNIE_Table100) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(rmse_PNIE_Table100,
          file = paste0(path, "/Tables/", "PNIE-RMSE_num_clust-100.csv"), row.names = FALSE)


# PNIE Coverage Table -----------------------------------------------------

# Pivot wide
cover_Tbl100 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 100, c("ICC",
                                                      "clust_size",
                                                      "analysisCond",
                                                      "coverage_PNIE",
                                                      "coverage_TNDE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("coverage_PNIE", "coverage_TNDE")
)

# drop common piece of column names
colnames(cover_Tbl100) <-
  stringr::str_remove(string = colnames(cover_Tbl100),
                      pattern = "coverage_")

# RMSE PNIE Table   
cover_PNIE_Table100 <- cover_Tbl100 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(cover_PNIE_Table100) %>%
  set_number_format(
    row = 1:nrow(cover_PNIE_Table100) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(cover_PNIE_Table100,
          file = paste0(path, "/Tables/", "PNIE-Coverage_num_clust-100.csv"), row.names = FALSE)

# PNIE Relative Bias Visual -----------------------------------------------

# set direct & indirect effects 
treat_m <- 1.17 # trt on med   
treat_y <- 1.35 # trt on outcome
med_y <- 1.2 # med on outcome  
TNDE <- treat_y
PNIE <- treat_m * med_y

# visual settings
gglayer_theme <- list(theme_bw(),
                      scale_fill_manual(values = c("#BF5700", #Fixed-effect
                                                   "#A6CD57", #Random-effect 
                                                   "#00a9b7", #Random-effect means 
                                                   "#333F48" )),#Single-level 
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

## Absolute Relative Bias Visual -----------------------------------------------

readRDS(file = paste0(path, "/Data/S1_Supp2_Simulation-Data.rds")) |> 
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                                  "ERROR"))))) %>% 
  # visual 
  ggplot(aes(x = abs((PNIE_est / PNIE) - 1), y = outModel, fill = `PS Model`)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.8, linewidth = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.3, 
               outlier.size = 0.5, 
               outlier.alpha = 0.5) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Absolute-Relative-Bias-Boxplot_num_clust-100.png"), plot = last_plot())



## Relative Bias Visual -----------------------------------------------

readRDS(file = paste0(path, "/Data/S1_Supp2_Simulation-Data.rds")) |> 
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                                  "ERROR"))))) %>% 
  # visual 
  ggplot(aes(x = (PNIE_est / PNIE) - 1, y = outModel, fill = `PS Model`)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.8, linewidth = 0.7) +
  geom_vline(xintercept = -0.10, color = "red", alpha = 0.8, linewidth = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.3, 
               outlier.size = 0.5, 
               outlier.alpha = 0.5) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot_num_clust-100.png"), plot = last_plot())


# PNIE MC CI Coverage Rate Visual -----------------------------------------

# add formatting & color scheme 
gglayer_theme_line <- list(theme_bw(),
                           scale_color_manual(values = c("#BF5700", #Fixed-effect
                                                         "#A6CD57", #Random-effect 
                                                         "#00a9b7", #Random-effect means 
                                                         "#333F48" )),#Single-level 
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
    x = "\n Cluster Size",
    y = "NIE Coverage Rate \n",
    color = "PS Model",
    linetype = "PS Model"#,
    # shape = "PS Model"
  ),
  guides(y.sec = guide_none("Mediator and Outcome Model"),
         x.sec = guide_none("Residual ICC")) #guide_none("Number of Clusters and Residual ICC"))
)

# separate PS & med/outcome models for labels 
perf_measure_DF |> 
  separate(analysisCond, sep = "_", into = c("PS", "Med", "Out")) |> 
  # rename PS models 
  mutate(PS = ifelse(PS == "FE", "Fixed-Effect", 
                     ifelse(PS == "RE", "Random-Effect", 
                            ifelse(PS == "SL", "Single-Level", 
                                   ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                          "ERROR"))))) |> 
  ggplot(aes(x = as.factor(clust_size), y = coverage_PNIE, color = PS, linetype = PS)) +
  geom_point() +
  geom_line(aes(group = PS)) +
  facet_grid(Out ~ ICC) +# facet_grid(Out ~ num_clust + ICC) +
  gglayer_theme_line +
  gglayer_labs

# Save plot 
ggsave(filename = paste0(path, "/Figures/", 
                         "PNIE-Coverage-Lineplot.png"), 
       plot = last_plot())

