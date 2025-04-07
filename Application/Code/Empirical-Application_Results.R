################################################################################
##################### QP Empirical Application - Results ######################
################################################################################

############################ Script Description ################################
#
# Author: 
# 
# Date Created: 2025-03-10
#
#
# Script Description:
#   This R script displays tables & visualizes the estimated effects to 
#     facilitate interpretation of the results.
# 
# Last Updated: 2025-03-19 
#
#
# Notes:
# 
################################################################################


# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  tidyverse, 
  ggplot2, 
  extrafont, 
  stringr, 
  mice, 
  cobalt, 
  WeightIt, 
  boot, 
  utils, 
  lme4, 
  WeMix, 
  parallel
)

# # Load Functions 
# source("Application/Functions/bootstrap_ci_paral_2.R")
# source("Application/Functions/bootstrap_ci_re_paral_2.R")
# source("Application/Functions/bootstrap_ci_re_mean_paral.R")
# source("Application/Functions/monteCarloCI.R")



# Import Data -----------------------------------------------------
# Load clean dataset 
# boot_ci_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_bootstrap-CIs.rds")
monte_ci_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_monte-carlo-CIs.rds")
estimates_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_noCIs.rds")

# 
# monte_ci_df <- monte_ci_df |> 
#   select("analysisCond":"PNDE_UL", "TNIE":"TNIE_UL")
# 
# boot_ci_df <- boot_ci_df |> 
#   select("cond", starts_with("PNDE"), starts_with("TNIE"), "PS", "Model") 



names(monte_ci_df)#; names(boot_ci_df)

head(monte_ci_df)#; head(boot_ci_df); head(estimates_df)

# Make model names match across dataframes 
monte_ci_df <- monte_ci_df |> 
  mutate(PS = ifelse(PS == "FE", "Fixed-Effect", 
                           ifelse(PS == "RE", "Random-Effect", 
                                  ifelse(PS == "SL", "Single-Level", 
                                         ifelse(PS == "RE-Mean", "Random-Effect with Cluster Means", 
                                                "ERROR"))))) |> 
  mutate(Model = ifelse(outModel == "FE", "Fixed-Effect", 
                             ifelse(outModel == "RE", "Random-Effect", 
                                    ifelse(outModel == "SL", "Single-Level", 
                                           ifelse(outModel == "RE-Mean", "Random-Effect with Cluster Means", 
                                                  "ERROR"))))) 

# 
monte_ci_df <- monte_ci_df |> 
  dplyr::select(!c("analysisCond", "medmodel", "outModel"))
  # select("PS", "Model", starts_with("PNDE"), starts_with("TNIE"))
# boot_ci_df <- boot_ci_df |> 
#   select(!c("cond"))
  # select("PS", "Model", starts_with("PNDE"), starts_with("TNIE"))

rownames(monte_ci_df) <- NULL
# rownames(boot_ci_df) <- NULL

# # Stack bootstrap CIs & monte carlo CIs
# results_df <- bind_rows(cbind(monte_ci_df, ci = "MC"),
#                         cbind(boot_ci_df, ci = "Bootstrap"))






# Result visuals ----------------------------------------------------------

# visual settings
gglayer_theme <- list(theme_bw(),
                      scale_color_manual(values = c("#BF5700", #Fixed-effect
                                                   "#A6CD57", #Random-effect 
                                                   # "#00a9b7", #Random-effect means 
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
    # x = "\n Absolute Relative Bias",
    y = "Propensity Score Model \n",
    color = "PS Model",
    linetype = "PS Model",
    shape = "PS Model"
  ),
  guides(#y.sec = guide_none("Cluster Size"),
         x.sec = guide_none("\n Mediator and Outcome Model"))
)

# ══════════════════════════════
#    TNDE 
# ══════════════════════════════
monte_ci_df |> 
  ggplot(aes(x = TNDE, y = PS, color = PS)) + #, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "TNDE Estimates with 95% CIs by PS and Outcome Model",
       x = "\n Total Natural Direct Effect (TNDE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  gglayer_theme +
  gglayer_labs +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Application/Output/Visuals/TNDE-Estimates.png", plot = last_plot())

# ══════════════════════════════
#    PNDE 
# ══════════════════════════════
monte_ci_df |> 
  ggplot(aes(x = PNDE, y = PS, color = PS)) + #, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "PNDE Estimates with 95% CIs by PS and Outcome Model",
       x = "\n Pure Natural Direct Effect (PNDE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  gglayer_theme +
  gglayer_labs +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Application/Output/Visuals/PNDE-Estimates.png", plot = last_plot())

# ══════════════════════════════
#    TNIE 
# ══════════════════════════════
monte_ci_df |> 
  ggplot(aes(x = TNIE, y = PS, color = PS)) + #, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "TNIE Estimates with 95% CIs by PS and Outcome Model",
       x = "\n Total Natural Indirect Effect (TNIE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  gglayer_theme +
  gglayer_labs +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Application/Output/Visuals/TNIE-Estimates.png", plot = last_plot())

# ══════════════════════════════
#    PNIE 
# ══════════════════════════════
monte_ci_df |> 
  ggplot(aes(x = PNIE, y = PS, color = PS)) + #, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "PNIE Estimates with 95% CIs by PS and Outcome Model",
       x = "\n Pure Natural Indirect Effect (PNIE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  gglayer_theme +
  gglayer_labs +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Application/Output/Visuals/PNIE-Estimates.png", plot = last_plot())




## compare bootstrap & MC --------------------------------------------------

# 
# # ══════════════════════════════
# #    TNDE 
# # ══════════════════════════════
# 
# # Create the TNDE plot
# results_df |>
#   filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
#   mutate(includes_zero = ifelse(TNDE_LL > 0 | TNDE_UL < 0, FALSE, TRUE)) |>
#   ggplot(aes(x = TNDE, y = PS, color = ci, group = ci, shape = includes_zero)) +
#   geom_point(position = position_dodge(width = 0.5), size = 2.75) +
#   geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL),
#                  position = position_dodge(width = 0.5),
#                  height = 0.25) +
#   facet_wrap(~ Model, ncol = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "TNDE Estimates with 95% CIs by PS and Outcome Model",
#        x = "Total Natural Direct Effect (TNDE)",
#        y = "Propensity Score Model",
#        color = "CI Method") +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12))
# 
# # Save 
# ggsave(filename = "Application/Output/Visuals/TNDE-CI-Comparison.png", plot = last_plot())
# 
# 
# # ══════════════════════════════
# #    PNIE 
# # ══════════════════════════════
# 
# # Create the PNIE plot
# results_df |> 
#   filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
#   mutate(includes_zero = ifelse(PNIE_LL > 0 | PNIE_UL < 0, FALSE, TRUE)) |> 
#   ggplot(aes(x = PNIE, y = PS, color = ci, group = ci, shape = includes_zero)) +
#   geom_point(position = position_dodge(width = 0.5), size = 2.75) +
#   geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL),
#                  position = position_dodge(width = 0.5),
#                  height = 0.25) +
#   facet_wrap(~ Model, ncol = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "PNIE Estimates with 95% CIs by PS and Outcome Model",
#        x = "Pure Natural Indirect Effect (PNIE)",
#        y = "Propensity Score Model",
#        color = "CI Method") +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12))
# 
# # Save 
# ggsave(filename = "Application/Output/Visuals/PNIE-CI-Comparison.png", plot = last_plot())
# 
# 
# 
# # ══════════════════════════════
# #    PNDE 
# # ══════════════════════════════
# 
# # Create the PNDE plot
# results_df |> 
#   filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
#   mutate(includes_zero = ifelse(PNDE_LL > 0 | PNDE_UL < 0, FALSE, TRUE)) |> 
#   ggplot(aes(x = PNDE, y = PS, color = ci, group = ci, shape = includes_zero)) +
#   geom_point(position = position_dodge(width = 0.5), size = 2.75) +
#   geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL),
#                  position = position_dodge(width = 0.5),
#                  height = 0.25) +
#   facet_wrap(~ Model, ncol = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "PNDE Estimates with 95% CIs by PS and Outcome Model",
#        x = "Pure Natural Direct Effect (PNDE)",
#        y = "Propensity Score Model",
#        color = "CI Method") +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12))
# 
# # Save 
# ggsave(filename = "Application/Output/Visuals/PNDE-CI-Comparison.png", plot = last_plot())
# 
# 
# 
# # ══════════════════════════════
# #    TNIE 
# # ══════════════════════════════
# 
# # Create the TNIE plot
# results_df |> 
#   filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
#   mutate(includes_zero = ifelse(TNIE_LL > 0 | TNIE_UL < 0, FALSE, TRUE)) |> 
#   ggplot(aes(x = TNIE, y = PS, color = ci, group = ci, shape = includes_zero)) +
#   geom_point(position = position_dodge(width = 0.5), size = 2.75) +
#   geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL),
#                  position = position_dodge(width = 0.5),
#                  height = 0.25) +
#   facet_wrap(~ Model, ncol = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "TNIE Estimates with 95% CIs by PS and Outcome Model",
#        x = "Total Natural Indirect Effect (TNIE)",
#        y = "Propensity Score Model",
#        color = "CI Method") +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12))
# 
# # Save 
# ggsave(filename = "Application/Output/Visuals/TNIE-CI-Comparison.png", plot = last_plot())


























# 
# # TNDE 
# ## Save visual 
# pdf("Application/Output/Visuals/TNDE-Estimates.pdf")
# ## Visual 
# boot_ci_df %>% 
#   mutate(
#     # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
#     Zero_Encompasses = ifelse(TNDE_LL > 0 | TNDE_UL < 0, "Below 0", "Includes 0")
#   ) %>% 
#   ggplot(aes(y = PS, x = TNDE)) +
#   geom_point(aes(color = Zero_Encompasses), size = 3) +
#   geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL, color = Zero_Encompasses), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "Total Natural Direct Effect (TNDE) with 95% Confidence Intervals",
#        x = "Total Natural Direct Effect (TNDE)",
#        y = "Propensity Score (PS)") +
#   facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
#   scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1),
#         legend.position = "none")  # Remove legend
# 
# ## 
# dev.off()
# # Save visual 
# ggsave(filename = "Application/Output/Visuals/TNDE-Estimates.png", plot = last_plot())
# 
# 
# # PNDE 
# ## Save visual 
# pdf("Application/Output/Visuals/PNDE-Estimates.pdf")
# ## Visual 
# boot_ci_df %>% 
#   mutate(
#     # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
#     Zero_Encompasses = ifelse(PNDE_LL > 0 | PNDE_UL < 0, "Below 0", "Includes 0")
#   ) %>% 
#   ggplot(aes(y = PS, x = PNDE)) +
#   geom_point(aes(color = Zero_Encompasses), size = 3) +
#   geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL, color = Zero_Encompasses), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "Pure Natural Direct Effect (PNDE) with 95% Confidence Intervals",
#        x = "Pure Natural Direct Effect (PNDE)",
#        y = "Propensity Score (PS)") +
#   facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
#   scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1),
#         legend.position = "none")  # Remove legend
# 
# ## 
# dev.off()
# # Save visual 
# ggsave(filename = "Application/Output/Visuals/PNDE-Estimates.png", plot = last_plot())
# 
# 
# 
# 
# # TNIE 
# ## Save visual 
# pdf("Application/Output/Visuals/TNIE-Estimates.pdf")
# ## Visual 
# boot_ci_df %>% 
#   mutate(
#     # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
#     Zero_Encompasses = ifelse(TNIE_LL > 0 | TNIE_UL < 0, "Below 0", "Includes 0")
#   ) %>% 
#   ggplot(aes(y = PS, x = TNIE)) +
#   geom_point(aes(color = Zero_Encompasses), size = 3) +
#   geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL, color = Zero_Encompasses), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "Total Natural Indirect Effect (TNIE) with 95% Confidence Intervals",
#        x = "Total Natural Indirect Effect (TNIE)",
#        y = "Propensity Score (PS)") +
#   facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
#   scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1),
#         legend.position = "none")  # Remove legend
# 
# ## 
# dev.off()
# # Save visual 
# ggsave(filename = "Application/Output/Visuals/TNIE-Estimates.png", plot = last_plot())
# 
# 
# # PNIE 
# ## Save visual 
# pdf("Application/Output/Visuals/PNIE-Estimates.pdf")
# ## Visual 
# boot_ci_df %>% 
#   mutate(
#     # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
#     Zero_Encompasses = ifelse(PNIE_LL > 0 | PNIE_UL < 0, "Below 0", "Includes 0")
#   ) %>% 
#   ggplot(aes(y = PS, x = PNIE)) +
#   geom_point(aes(color = Zero_Encompasses), size = 3) +
#   geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL, color = Zero_Encompasses), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
#   labs(title = "Pure Natural Indirect Effect (PNIE) with 95% Confidence Intervals",
#        x = "Pure Natural Indirect Effect (PNIE)",
#        y = "Propensity Score (PS)") +
#   facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
#   scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1),
#         legend.position = "none")  # Remove legend
# 
# ## 
# dev.off()
# # Save visual 
# ggsave(filename = "Application/Output/Visuals/PNIE-Estimates.png", plot = last_plot())
# 
