################################################################################
##################### QP Empirical Application - Results ######################
################################################################################

############################ Script Description ################################
#
# Author: Cameron
# 
# Date Created: 2025-03-10
#
#
# Script Description:
#   This R script displays tables & visualizes the estimated effects to 
#     facilitate interpretation of the results.
# 
# Last Updated: 2025-03-13 
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
boot_ci_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_bootstrap-CIs.rds")
monte_ci_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_monte-carlo-CIs.rds")
estimates_df <- read_rds(file = "Application/Output/Estimates/Effect-Estimates_noCIs.rds")

# 
# monte_ci_df <- monte_ci_df |> 
#   select("analysisCond":"PNDE_UL", "TNIE":"TNIE_UL")
# 
# boot_ci_df <- boot_ci_df |> 
#   select("cond", starts_with("PNDE"), starts_with("TNIE"), "PS", "Model") 



names(monte_ci_df); names(boot_ci_df)

head(monte_ci_df); head(boot_ci_df); head(estimates_df)

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
  select(!c("analysisCond", "medmodel", "outModel"))
  # select("PS", "Model", starts_with("PNDE"), starts_with("TNIE"))
boot_ci_df <- boot_ci_df |> 
  select(!c("cond"))
  # select("PS", "Model", starts_with("PNDE"), starts_with("TNIE"))

rownames(monte_ci_df) <- NULL
rownames(boot_ci_df) <- NULL

# Stack bootstrap CIs & monte carlo CIs
results_df <- bind_rows(cbind(monte_ci_df, ci = "MC"),
                        cbind(boot_ci_df, ci = "Bootstrap"))






# Result visuals ----------------------------------------------------------

# ══════════════════════════════
#    TNDE 
# ══════════════════════════════

# Create the TNDE plot
results_df |>
  filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
  mutate(includes_zero = ifelse(TNDE_LL > 0 | TNDE_UL < 0, FALSE, TRUE)) |>
  ggplot(aes(x = TNDE, y = PS, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "TNDE Estimates with 95% CIs by PS and Outcome Model",
       x = "Total Natural Direct Effect (TNDE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Save 
ggsave(filename = "Application/Output/Visuals/TNDE-CI-Comparison.png", plot = last_plot())


# ══════════════════════════════
#    PNIE 
# ══════════════════════════════

# Create the PNIE plot
results_df |> 
  filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
  mutate(includes_zero = ifelse(PNIE_LL > 0 | PNIE_UL < 0, FALSE, TRUE)) |> 
  ggplot(aes(x = PNIE, y = PS, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "PNIE Estimates with 95% CIs by PS and Outcome Model",
       x = "Pure Natural Indirect Effect (PNIE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Save 
ggsave(filename = "Application/Output/Visuals/PNIE-CI-Comparison.png", plot = last_plot())



# ══════════════════════════════
#    PNDE 
# ══════════════════════════════

# Create the PNDE plot
results_df |> 
  filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
  mutate(includes_zero = ifelse(PNDE_LL > 0 | PNDE_UL < 0, FALSE, TRUE)) |> 
  ggplot(aes(x = PNDE, y = PS, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "PNDE Estimates with 95% CIs by PS and Outcome Model",
       x = "Pure Natural Direct Effect (PNDE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Save 
ggsave(filename = "Application/Output/Visuals/PNDE-CI-Comparison.png", plot = last_plot())



# ══════════════════════════════
#    TNIE 
# ══════════════════════════════

# Create the TNIE plot
results_df |> 
  filter(Model %in% c("Single-Level", "Fixed-Effect")) |> # Drop RE & RE-Mean for now 
  mutate(includes_zero = ifelse(TNIE_LL > 0 | TNIE_UL < 0, FALSE, TRUE)) |> 
  ggplot(aes(x = TNIE, y = PS, color = ci, group = ci, shape = includes_zero)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.75) +
  geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL),
                 position = position_dodge(width = 0.5),
                 height = 0.25) +
  facet_wrap(~ Model, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "TNIE Estimates with 95% CIs by PS and Outcome Model",
       x = "Total Natural Indirect Effect (TNIE)",
       y = "Propensity Score Model",
       color = "CI Method") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Save 
ggsave(filename = "Application/Output/Visuals/TNIE-CI-Comparison.png", plot = last_plot())























results_df |> 
  filter(Model == "Single-Level") |> 
  select(PS, Model, ci, everything()) |> 
  arrange(PS)
#             PS        Model        ci       TNDE    TNDE_LL      TNDE_UL        PNIE    PNIE_LL    PNIE_UL      PNDE   PNDE_LL  PNDE_UL        TNIE    TNIE_LL    TNIE_UL
# 1  Fixed-Effect Single-Level        MC -0.2588860 -0.5676262  0.045215270 -0.12127187 -0.3890650 0.12039158 0.9966576 -4.764792 6.905444 -0.11841160 -0.3551257 0.12457045
# 2  Fixed-Effect Single-Level Bootstrap -0.2588860 -0.5103306 -0.003194003 -0.12127187 -0.3246048 0.07121469 0.9966576 -4.559389 5.770201 -0.11841160 -0.3137209 0.06910764
# 3 Random-Effect Single-Level        MC -0.2987412 -0.6059235  0.004602053 -0.09161689 -0.3569932 0.14711015 0.9023210 -4.799629 6.817355 -0.08952703 -0.3255360 0.15508453
# 4 Random-Effect Single-Level Bootstrap -0.2987412 -0.5280205 -0.062530157 -0.09161689 -0.2646909 0.08650299 0.9023210 -4.408213 5.261863 -0.08952703 -0.2615859 0.08346661
# 5  Single-Level Single-Level        MC -0.3259836 -0.6332298 -0.023838593 -0.04471103 -0.3052875 0.19400418 0.2311803 -5.464905 6.104233 -0.04422288 -0.2846380 0.20340117
# 6  Single-Level Single-Level Bootstrap -0.3259836 -0.5633203 -0.085134771 -0.04471103 -0.2067034 0.13493934 0.2311803 -4.993651 4.450986 -0.04422288 -0.2057100 0.13509154


results_df |> 
  filter(Model == "Fixed-Effect") |> 
  select(PS, Model, ci, everything()) |> 
  arrange(PS)
#             PS        Model        ci       TNDE    TNDE_LL      TNDE_UL       PNIE    PNIE_LL    PNIE_UL     PNDE   PNDE_LL  PNDE_UL       TNIE    TNIE_LL    TNIE_UL
# 1  Fixed-Effect Fixed-Effect        MC -0.2510769 -0.5683398  0.046112486 -0.1195056 -0.3838922 0.11772283 1.612312 -3.923283 7.679246 -0.1152346 -0.3578365 0.11800255
# 2  Fixed-Effect Fixed-Effect Bootstrap -0.2510769 -0.5060783 -0.005091589 -0.1195056 -0.3255577 0.07700916 1.612312 -3.852802 6.412009 -0.1152346 -0.3300916 0.07408897
# 3 Random-Effect Fixed-Effect        MC -0.2562645 -0.5766853  0.047863405 -0.1158807 -0.3804614 0.12229422 1.434864 -4.196621 7.567613 -0.1121335 -0.3539609 0.12542032
# 4 Random-Effect Fixed-Effect Bootstrap -0.2562645 -0.5032962 -0.012264249 -0.1158807 -0.3126527 0.07480969 1.434864 -3.736447 6.021261 -0.1121335 -0.3061472 0.07332897
# 5  Single-Level Fixed-Effect        MC -0.2325015 -0.5686607  0.083536920 -0.1181606 -0.3976647 0.13417959 0.475951 -5.302737 6.503549 -0.1165163 -0.3712595 0.13766733
# 6  Single-Level Fixed-Effect Bootstrap -0.2325015 -0.4810479  0.006200502 -0.1181606 -0.3050444 0.07659792 0.475951 -4.684792 4.942427 -0.1165163 -0.2928554 0.07527331





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
