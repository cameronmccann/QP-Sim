################################################################################
##################### QP Simulation 1 Supplemental Results #####################
################################################################################

############################ Script Description ################################
#
# Author: Cameron McCann
# 
# Date Created: 2025-02-11
#
#
# Script Description: This code summarizes and reports the results for the 
#                       supplemental for first simulation study (i.e., obtains performance measures).  
#                       This is stored in the relevant Results folder.
#
#
# Last Updated: 2025-02-13
#
#
# Notes:
#   To-Do
#       # Add coverage table & visual (started viz under "PNIE MC CI Coverage Rate Visual")
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
path <- "Output/S1_Supp_Results/2025-02-11-test_200-reps"
retrieval_path <- "Output/S1_Supp_Simulation-Output/2025-02-11-test_200-reps"
dir.create(path = "Output/S1_Supp_Results")
dir.create(path = path)
dir.create(path = paste0(path, "/Data"))
dir.create(path = paste0(path, "/Tables"))
dir.create(path = paste0(path, "/Figures"))

# Import data  ------------------------------------------------------------

# List all files matching the pattern.
file_list <- list.files(path = retrieval_path,
                        pattern = "^S1_Supp_Condition-[0-9]+-Overall_Estimates_.*\\.rds$",
                        full.names = TRUE)
# Read each file and combine into one data frame
sim1_data <- do.call(rbind, lapply(file_list, readRDS))


# Clean Data --------------------------------------------------------------

# Change row names 
rownames(sim1_data) <- 1:nrow(sim1_data)
# Change to numeric & properly name effects 
sim1_data <- sim1_data %>% 
  mutate_at(.vars = c("NDE_est", "NIE_est", "ICC", "clust_size", "conditionNum", 
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
  group_by(ICC, clust_size, 
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
          file = paste0(path, "/Data/", "S1_Supp_Performance-Measures.rds"))

write_rds(sim1_data, 
          file = paste0(path, "/Data/", "S1_Supp_Simulation-Data.rds"))


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
as_hux(relBias_TNDE_Table) %>% 
  set_number_format(row = 1:nrow(relBias_TNDE_Table) + 1, col = -c(1:2), value = 3)

# Save csv 
write.csv(relBias_TNDE_Table,
          file = paste0(path, "/Tables/", "TNDE-Relative-Bias.csv"), row.names = FALSE)


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
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
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

# Save csv 
write.csv(relBias_PNIE_Table,
          file = paste0(path, "/Tables/", "PNIE-Relative-Bias.csv"), row.names = FALSE)


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
  dplyr::select(analysisCond, starts_with("TNDE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
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

# Save csv 
write.csv(rmse_TNDE_Table,
          file = paste0(path, "/Tables/", "TNDE-RMSE.csv"), row.names = FALSE)


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
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
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

# Save csv 
write.csv(rmse_PNIE_Table,
          file = paste0(path, "/Tables/", "PNIE-RMSE.csv"), row.names = FALSE)


# PNIE Coverage Table -----------------------------------------------------

# Pivot wide
cover_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "coverage_PNIE",
                      "coverage_TNDE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("coverage_PNIE", "coverage_TNDE")
)

# drop common piece of column names
colnames(cover_Tbl) <-
  stringr::str_remove(string = colnames(cover_Tbl),
                      pattern = "coverage_")

# RMSE PNIE Table   
cover_PNIE_Table <- cover_Tbl %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(cover_PNIE_Table) %>%
  set_number_format(
    row = 1:nrow(cover_PNIE_Table) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(cover_PNIE_Table,
          file = paste0(path, "/Tables/", "PNIE-Coverage-Rate.csv"), row.names = FALSE)


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

# ═══════════════════
#    cluster size = 20 
# ═══════════════════
# Boxplot 
readRDS(file = paste0(path, "/Data/S1_Supp_Simulation-Data.rds")) |> 
  filter(clust_size == 20) |>
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                                  "ERROR"))))) |> 
  # visual 
  ggplot(aes(x = abs((PNIE_est / PNIE) - 1), y = interaction(medmodel, outModel, sep = " - "), fill = `PS Model`)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.6, linewidth = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.3, 
               outlier.size = 0.5, 
               outlier.alpha = 0.6) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot_clust-size-20.png"), plot = last_plot())

# ═══════════════════
#    cluster size = 40 
# ═══════════════════
# Boxplot 
readRDS(file = paste0(path, "/Data/S1_Supp_Simulation-Data.rds")) |> 
  filter(clust_size == 40) |>
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                                  "ERROR"))))) |> 
  # visual 
  ggplot(aes(x = abs((PNIE_est / PNIE) - 1), y = interaction(medmodel, outModel, sep = " - "), fill = `PS Model`)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.6, linewidth = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.3, 
               outlier.size = 0.5, 
               outlier.alpha = 0.6) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot_clust-size-40.png"), plot = last_plot())

# ═══════════════════
#    cluster size = 100 
# ═══════════════════
# Boxplot 
readRDS(file = paste0(path, "/Data/S1_Supp_Simulation-Data.rds")) |> 
  filter(clust_size == 100) |>
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           ifelse(PS == "RE-Mean", "Random-Effect Mean", 
                                                  "ERROR"))))) |> 
  # visual 
  ggplot(aes(x = abs((PNIE_est / PNIE) - 1), y = interaction(medmodel, outModel, sep = " - "), fill = `PS Model`)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.6, linewidth = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.3, 
               outlier.size = 0.5, 
               outlier.alpha = 0.6) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot_clust-size-100.png"), plot = last_plot())


# readRDS(file = paste0(path, "/Data/S1_Supp_Simulation-Data.rds")) |> 
#   # rename PS models 
#   mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
#                              ifelse(PS == "RE", "Random-Effect", 
#                                     ifelse(PS == "SL", "Single-Level", 
#                                            ifelse(PS == "RE-Mean", "Random-Effect Mean", 
#                                            "ERROR"))))) %>% 
#   # visual 
#   ggplot(aes(x = abs((PNIE_est / PNIE) - 1), y = outModel, fill = `PS Model`)) +
#   geom_vline(xintercept = 0, linewidth = 0.5) +
#   geom_vline(xintercept = 0.10, color = "red", alpha = 0.6, linewidth = 0.5) +
#   geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
#                linewidth = 0.3, 
#                outlier.size = 0.5, 
#                outlier.alpha = 0.6) +
#   facet_grid(clust_size ~ ICC) +
#   gglayer_labs +
#   gglayer_theme 
# 
# # Save visual 
# ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot.png"), plot = last_plot())




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
    x = "\n Absolute Relative Bias",
    y = "Mediator and Outcome Model \n",
    color = "PS Model",
    linetype = "PS Model",
    shape = "PS Model"
  ),
  guides(y.sec = guide_none("Cluster Size"),
         x.sec = guide_none("Residual ICC"))
)

# separate PS & med/outcome models for labels 
perf_measure_DF |> 
  separate(analysisCond, sep = "_", into = c("PS", "Med", "Out")) |> 
  filter(Out == "FE") |> 
  ggplot(aes(x = as.factor(clust_size), y = coverage_PNIE, color = PS, alpha = 0.9#, shape = Med #,linetype = Med
             )) +
  geom_point() +
  geom_line(aes(group = interaction(PS, Med))) +
  facet_grid(Med + Out ~ ICC) +
  gglayer_theme_line +
  guides(y.sec = guide_none("Mediator and Outcome Models"), x.sec = guide_none("Residual ICC")) +
  labs(y = "NIE Coverage Rate", 
       x = "Cluster Size")
  # gglayer_labs # NEED TO UPDATE LABELS 
# Save plot 
ggsave(filename = paste0(path, "/Figures/", 
                         "PNIE-Coverage-Lineplot.png"), 
       plot = last_plot())


## Examining SL Outcome further --------------------------------------------

readRDS(file = paste0(path, "/Data/S1_Supp_Simulation-Data.rds")) |> 
  filter(ICC == 0.2, clust_size == 40) |> 
  filter(PS %in% "SL", medmodel %in% "SL", outModel %in% "SL") |> 
  mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
  ggplot(aes(x = rep, y = PNIE_est)) +
  geom_hline(yintercept = PNIE) +
  geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL), width = 0.2, color = "darkgray") +
  geom_point(size = 2, color = "blue") 
# coverage is 0

# ## Examining FE Outcome further --------------------------------------------
# 
# # DF for coverage (used for text)
# coverDF <- perf_measure_DF |> 
#   separate(analysisCond, into = c("PS", "Med", "Out"), sep = "_") |> 
#   filter(ICC == 0.2, Out == "FE") |> 
#   mutate(coverage_PNIE = format(coverage_PNIE, digits = 3))
# 
# # FE
# coverFEDF <- readRDS(file = paste0(path, "/Data/S1_Simulation-Data.rds")) |> 
#     filter(ICC == 0.2) |> #, clust_size == 40) |> 
#     filter(PS %in% c("FE"), medmodel %in% "FE", outModel %in% "FE") |> 
#     mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
#     mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
#            if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
#            error_color = ifelse(if_cover_PNIE == TRUE, "#BF5700", "darkgray")) |>
#     arrange(PNIE_est) |> 
#     group_by(clust_size) |> 
#     mutate(order = row_number())
# coverPlotFE <- coverFEDF |>
#     ggplot(aes(x = order, y = PNIE_est)) +
#     geom_hline(yintercept = PNIE) +
#     # Map error_color to the error bars
#     geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
#     # Fix point color to match the “matched” color
#     geom_point(size = 0.5, color = "#BF5700") +
#     scale_color_identity() +
#     scale_x_continuous(breaks = seq(0, max(coverFEDF$order, na.rm = TRUE), by = 50)) +
#     # theme_minimal() +
#     gglayer_theme +
#     theme(axis.text.x.bottom = element_blank(), 
#           axis.ticks = element_blank()) +
#     labs(x = "", 
#          y = "") +
#     facet_grid(PS ~ clust_size) +
#     geom_text(
#       data = coverDF[coverDF$PS == "FE", ], aes(x = Inf, y = Inf, label = paste0(coverage_PNIE)), 
#       hjust = 4.5, vjust = 3, size = 3, color = "black") +
#   guides(x.sec = guide_none("Cluster Size"))
# 
# # RE
# coverREDF <- readRDS(file = paste0(path, "/Data/S1_Simulation-Data.rds")) |> 
#   filter(ICC == 0.2) |> #, clust_size == 40) |> 
#   filter(PS %in% c("RE"), medmodel %in% "FE", outModel %in% "FE") |> 
#   mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
#   mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
#          if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
#          error_color = ifelse(if_cover_PNIE == TRUE, "#A6CD57", "darkgray")) |>
#   arrange(PNIE_est) |> 
#   group_by(clust_size) |> 
#   mutate(order = row_number())
# coverPlotRE <- coverREDF |> 
#   ggplot(aes(x = order, y = PNIE_est)) +
#   geom_hline(yintercept = PNIE) +
#   # Map error_color to the error bars
#   geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
#   # Fix point color to match the “matched” color
#   geom_point(size = 0.5, color = "#A6CD57") +
#   scale_color_identity() +
#   scale_x_continuous(breaks = seq(0, max(coverREDF$order, na.rm = TRUE), by = 50)) +
#   # theme_minimal() +
#   gglayer_theme +
#   theme(axis.text.x.bottom = element_blank(), 
#         axis.ticks = element_blank()) +
#   labs(x = "", 
#        y = "NIE") +
#   facet_grid(PS ~ clust_size) +
#   geom_text(
#     data = coverDF[coverDF$PS == "RE", ], aes(x = Inf, y = Inf, label = paste0(coverage_PNIE)), 
#     hjust = 4.5, vjust = 3, size = 3, color = "black") +
#   guides(y.sec = guide_none("PS Model"))
# 
# # SL
# coverSLDF <- readRDS(file = paste0(path, "/Data/S1_Simulation-Data.rds")) |> 
#   filter(ICC == 0.2) |> #, clust_size == 40) |> 
#   filter(PS %in% c("SL"), medmodel %in% "FE", outModel %in% "FE") |> 
#   mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
#   mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
#          if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
#          error_color = ifelse(if_cover_PNIE == TRUE, "#333F48", "darkgray")) |>
#   arrange(PNIE_est) |> 
#   group_by(clust_size) |> 
#   mutate(order = row_number())
# coverPlotSL <- coverSLDF |> 
#   ggplot(aes(x = order, y = PNIE_est)) +
#   geom_hline(yintercept = PNIE) +
#   # Map error_color to the error bars
#   geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
#   # Fix point color to match the “matched” color
#   geom_point(size = 0.5, color = "#333F48") +
#   scale_color_identity() +
#   scale_x_continuous(breaks = seq(0, max(coverSLDF$order, na.rm = TRUE), by = 50)) +
#   # theme_minimal() +
#   gglayer_theme +
#   # theme(axis.text.x.bottom = element_blank()) +
#   labs(x = "", 
#        y = "") +
#   facet_grid(PS ~ clust_size) +
#   geom_text(
#     data = coverDF[coverDF$PS == "SL", ], aes(x = Inf, y = Inf, label = paste0(coverage_PNIE)), 
#     hjust = 4.5, vjust = 3, size = 3, color = "black", inherit.aes = TRUE) 
# 
# # Plot visual 
# coverPlotFE / coverPlotRE / coverPlotSL +
#   patchwork::plot_annotation(title = "Coverage rate across replications by PS & cluster size for ICC = 0.2 and med/outcome models = FE") 
# # Save plot 
# ggsave(filename = paste0(path, "/Figures/", 
#                          "S1_coverage-per-rep-by-PS-model-and-cluster-size-for-icc-0.2-and-outcome-FE.png"), 
#        plot = last_plot())


