################################################################################
################### QP Simulation 1 Supplemental 3 Results #####################
################################################################################

############################ Script Description ################################
#
# Author: 
# 
# Date Created: 2025-03-02
#
#
# Script Description: This code summarizes and reports the results for 
#                       Simulation Study 1 Supplemental 3, where the 
#                       unmeasured cluster-level confounding is weakened 
#                       relative to Simulation Study 1.  
#                       This is stored in the relevant Results folder.
#
#
# Last Updated: 2025-03-30
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


# Simulation 1 Explained Variance Table -----------------------------------

# Set simulation parameters
num_x = 3
# Coefficients for X variables' effects
a_x <- 0.25  # Effect of X on treatment (trt), corresponding R^2 = 0.02798
b_x <- 0.3   # Effect of X on mediator (med), corresponding R^2 = 0.02994
c_x <- 0.35  # Effect of X on outcome, corresponding R^2 = 0.03030

# Coefficients for Z variable's effects
a_z <- 1.03  # Effect of Z on treatment, corresponding R^2 = 0.19790
b_z <- 1.21  # Effect of Z on mediator, corresponding R^2 = 0.20292
c_z <- 1.4   # Effect of Z on outcome, corresponding R^2 = 0.20201

# Coefficients for treatment effects
treat_m <- 1.17  # Effect of treatment on mediator, corresponding R^2 = 0.15178
treat_y <- 1.35  # Effect of treatment on outcome, corresponding R^2 = 0.15027
med_y <- 1.2     # Effect of mediator on outcome, corresponding R^2 = 0.11873


# Create an empty data frame to store explained variance (R^2) values for different ICC levels
design_values_DF <- data.frame(cbind(
  var_expl = c("trt_x", "med_x", "out_x", 
               "trt_z", "med_z", "out_z", 
               "med_trt", "out_trt", "out_med"), 
  icc_0.05 = NA, 
  icc_0.2 = NA, 
  icc_0.5 = NA 
))

# Loop through ICC values and calculate the explained variance (R^2) for each variable
for (icc in c(0.05, 0.2, 0.5)) {
  # Variance explained by X variables (trt_x, med_x, out_x)
  design_values_DF[design_values_DF$var_expl == "trt_x", paste0("icc_", icc)] <-
    ((a_x^2)*num_x * (1 - icc)) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "med_x", paste0("icc_", icc)] <-
    ((b_x^2)*num_x * (1 - icc)) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_x", paste0("icc_", icc)] <-
    ((c_x^2)*num_x * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by Z variable (trt_z, med_z, out_z)
  design_values_DF[design_values_DF$var_expl == "trt_z", paste0("icc_", icc)] <-
    (a_z^2) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "med_z", paste0("icc_", icc)] <-
    (b_z^2) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_z", paste0("icc_", icc)] <-
    (c_z^2) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by treatment (med_trt, out_trt) and mediator (out_med)
  design_values_DF[design_values_DF$var_expl == "med_trt", paste0("icc_", icc)] <-
    (treat_m^2 * (1 - icc)) / (b_z^2 + treat_m^2 + (b_x^2)*num_x + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_trt", paste0("icc_", icc)] <-
    (treat_y^2 * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "out_med", paste0("icc_", icc)] <-
    (med_y^2 * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3))  
}

# Round the R^2 values for readability
design_values_DF <- design_values_DF %>% 
  mutate_at(.vars = c("icc_0.05", "icc_0.2", "icc_0.5"), 
            .funs = as.numeric) %>% 
  mutate_at(.vars = c("icc_0.05", "icc_0.2", "icc_0.5"), 
            .funs = round, 2)

# Remove rows related to Z variable as these are not needed for the final table
design_values_DF <- design_values_DF[!endsWith(design_values_DF$var_expl, suffix = "_z"), ]

# Create final table summarizing the explained variance for different factors
explained_var_table <- data.frame(
  cbind(
    Factor = c(
      "Total variance in the treatment explained by all X variables", 
      "Total variance in the mediator explained by all X variables", 
      "Total variance in the outcome explained by all X variables",
      
      "Total variance in the mediator explained by the treatment", 
      "Total variance in the outcome explained by the treatment", 
      "Total variance in the outcome explained by the mediator"
    ),
    "0.05" = design_values_DF$icc_0.05,
    "0.2" = design_values_DF$icc_0.2,
    "0.5" = design_values_DF$icc_0.5
  )
)

explained_var_table1 <- explained_var_table
explained_var_table1


# what new Z values to use ------------------------------------------------

# ####### NEW CODE #######
# # We'll keep everything else the same, just solve for new a_z, b_z, and c_z
# 
# # Desired R^2 for each Z effect:
# R2_target <- 0.05
# 
# # Known constants:
# num_x <- 3
# 
# # X-effect coefficients (unchanged):
# a_x <- 0.25
# b_x <- 0.3
# c_x <- 0.35
# 
# # Treatment effects (unchanged):
# treat_m <- 1.17
# treat_y <- 1.35
# med_y   <- 1.2
# 
# # For convenience, define the small constants used repeatedly:
# pi2_3     <- (pi^2 / 3)         # pi^2/3
# quarter   <- (pi^2 / 3) / 4     # (pi^2/3)/4
# 
# ###################
# # Solve for a_z so that R^2_{trt_z} ~ 0.05
# # R^2_{trt_z} = a_z^2 / [ (a_x^2)*num_x + a_z^2 + quarter + pi2_3 ]
# ###################
# {
#   denom_without_a_z <- (a_x^2)*num_x + quarter + pi2_3
#   
#   # Solve  R2_target = a_z^2 / (a_z^2 + denom_without_a_z)
#   # => a_z^2 = R2_target * (a_z^2 + denom_without_a_z)
#   # => a_z^2 - R2_target * a_z^2 = R2_target * denom_without_a_z
#   # => a_z^2 (1 - R2_target) = R2_target * denom_without_a_z
#   # => a_z^2 = [R2_target * denom_without_a_z] / (1 - R2_target)
#   a_z_sq <- (R2_target * denom_without_a_z) / (1 - R2_target)
#   a_z_new <- sqrt(a_z_sq)
# }
# 
# ###################
# # Solve for b_z so that R^2_{med_z} ~ 0.05
# # R^2_{med_z} = b_z^2 / [ (b_x^2)*num_x + b_z^2 + treat_m^2 + quarter + pi2_3 ]
# ###################
# {
#   denom_without_b_z <- (b_x^2)*num_x + treat_m^2 + quarter + pi2_3
#   
#   b_z_sq <- (R2_target * denom_without_b_z) / (1 - R2_target)
#   b_z_new <- sqrt(b_z_sq)
# }
# 
# ###################
# # Solve for c_z so that R^2_{out_z} ~ 0.05
# # R^2_{out_z} = c_z^2 / [ c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + quarter + pi2_3 ]
# ###################
# {
#   denom_without_c_z <- treat_y^2 + (c_x^2)*num_x + med_y^2 + quarter + pi2_3
#   
#   c_z_sq <- (R2_target * denom_without_c_z) / (1 - R2_target)
#   c_z_new <- sqrt(c_z_sq)
# }
# 
# # Print out the new coefficients that give R^2 ~ 0.05
# cat("New a_z =", round(a_z_new, 3), "\n")
# cat("New b_z =", round(b_z_new, 3), "\n")
# cat("New c_z =", round(c_z_new, 3), "\n")




# Simulation 1 Supp C Explained Variance Table -----------------------------------

# Set simulation parameters
num_x = 3
# Coefficients for X variables' effects
a_x <- 0.25  # Effect of X on treatment (trt), corresponding R^2 = 0.02798
b_x <- 0.3   # Effect of X on mediator (med), corresponding R^2 = 0.02994
c_x <- 0.35  # Effect of X on outcome, corresponding R^2 = 0.03030

# Coefficients for Z variable's effects
a_z <- 0.476 # Effect of Z on treatment, corresponding R^2 = 0.05006  # 1.03  # Effect of Z on treatment, corresponding R^2 = 0.19790
b_z <- 0.55  # Effect of Z on mediator, corresponding R^2 = 0.04997   # 1.21  # Effect of Z on mediator, corresponding R^2 = 0.20292
c_z <- 0.638 # Effect of Z on outcome, corresponding R^2 = 0.04995    # 1.4   # Effect of Z on outcome, corresponding R^2 = 0.20201

# Coefficients for treatment effects
treat_m <- 1.17  # Effect of treatment on mediator, corresponding R^2 = 0.15178
treat_y <- 1.35  # Effect of treatment on outcome, corresponding R^2 = 0.15027
med_y <- 1.2     # Effect of mediator on outcome, corresponding R^2 = 0.11873


# Create an empty data frame to store explained variance (R^2) values for different ICC levels
design_values_DF <- data.frame(cbind(
  var_expl = c("trt_x", "med_x", "out_x", 
               "trt_z", "med_z", "out_z", 
               "med_trt", "out_trt", "out_med"), 
  icc_0.05 = NA, 
  icc_0.2 = NA, 
  icc_0.5 = NA 
))

# Loop through ICC values and calculate the explained variance (R^2) for each variable
for (icc in c(0.05, 0.2, 0.5)) {
  # Variance explained by X variables (trt_x, med_x, out_x)
  design_values_DF[design_values_DF$var_expl == "trt_x", paste0("icc_", icc)] <-
    ((a_x^2)*num_x * (1 - icc)) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "med_x", paste0("icc_", icc)] <-
    ((b_x^2)*num_x * (1 - icc)) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_x", paste0("icc_", icc)] <-
    ((c_x^2)*num_x * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by Z variable (trt_z, med_z, out_z)
  design_values_DF[design_values_DF$var_expl == "trt_z", paste0("icc_", icc)] <-
    (a_z^2) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "med_z", paste0("icc_", icc)] <-
    (b_z^2) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_z", paste0("icc_", icc)] <-
    (c_z^2) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by treatment (med_trt, out_trt) and mediator (out_med)
  design_values_DF[design_values_DF$var_expl == "med_trt", paste0("icc_", icc)] <-
    (treat_m^2 * (1 - icc)) / (b_z^2 + treat_m^2 + (b_x^2)*num_x + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_trt", paste0("icc_", icc)] <-
    (treat_y^2 * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "out_med", paste0("icc_", icc)] <-
    (med_y^2 * (1 - icc)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3))  
}

# Round the R^2 values for readability
design_values_DF <- design_values_DF %>% 
  mutate_at(.vars = c("icc_0.05", "icc_0.2", "icc_0.5"), 
            .funs = as.numeric) %>% 
  mutate_at(.vars = c("icc_0.05", "icc_0.2", "icc_0.5"), 
            .funs = round, 2)

# Remove rows related to Z variable as these are not needed for the final table
design_values_DF <- design_values_DF[!endsWith(design_values_DF$var_expl, suffix = "_z"), ]

# Create final table summarizing the explained variance for different factors
explained_var_table <- data.frame(
  cbind(
    Factor = c(
      "Total variance in the treatment explained by all X variables", 
      "Total variance in the mediator explained by all X variables", 
      "Total variance in the outcome explained by all X variables",
      
      "Total variance in the mediator explained by the treatment", 
      "Total variance in the outcome explained by the treatment", 
      "Total variance in the outcome explained by the mediator"
    ),
    "0.05" = design_values_DF$icc_0.05,
    "0.2" = design_values_DF$icc_0.2,
    "0.5" = design_values_DF$icc_0.5
  )
)

explained_var_table



# ══════════════════════════════
#    NEED TO ADJUST EVERYTHING BELOW; AND COD SUMMARY AT TOP OF SCRIPT 
# ══════════════════════════════


# Simulation 1 Directed Acyclic Graph (DAG) Visualization -----------------

# Define the DAG structure and coordinates for each variable
dag1 <- dagify(
  Z ~ 1,              # Latent variable Z (no parent)
  x1 ~ 1, x2 ~ 1, x3 ~ 1,  # Independent variables x1, x2, x3 (no parents)
  T ~ 1,              # Treatment variable T (no parent)
  M ~ 1,              # Mediator variable M (no parent)
  Y ~ 1,              # Outcome variable Y (no parent)
  exposure = "T",     # Define T as the exposure
  outcome = "Y",      # Define Y as the outcome
  latent = "Z",       # Define Z as a latent variable
  coords = list(      # Coordinates for plotting the DAG nodes
    x = c(T = 1.5, M = 3.5, Y = 5.5, Z = 1, x1 = 2.75, x2 = 3.5, x3 = 4.25),
    y = c(T = 2, M = 3.5, Y = 2, Z = 6.5, x1 = 5, x2 = 5, x3 = 5)
  )
)


# Create DAG visualization
p1 <- ggdag_classic(dag1, size = 4) +
  
  # Set plot limits to ensure all elements fit well within the plot area
  xlim(0, 6) +
  ylim(1, 7) +
  coord_fixed() +  # Maintain fixed aspect ratio
  
  #----------------------
  # Z Section (Latent Variable Z)
  #----------------------
  
  # Add Z (latent variable) as a circle at specified coordinates
  ggforce::stat_circle(aes(x0 = 1, y0 = 6.5, r = 0.3)) +
  
  # Add dashed arrows from Z to T, M, and Y
  geom_curve(x = 0.65, y = 6.3, 
             xend = 1.25, yend = 2.3, 
             linetype = "dashed", 
             arrow = arrow(length = unit(0.2, "cm"), 
                           ends = "last", 
                           type = "closed"), 
             curvature = 0.4) + # Z -> T
  geom_curve(x = 1, y = 6.1, 
             xend = 3.2, yend = 3.5, 
             linetype = "dashed", 
             arrow = arrow(length = unit(0.2, "cm"), 
                           ends = "last", 
                           type = "closed"), 
             curvature = 0.4) + # Z -> M
  geom_curve(x = 1.3, y = 6.75, 
             xend = 5.65, yend = 2.3, 
             linetype = "dashed", 
             arrow = arrow(length = unit(0.2, "cm"), 
                           ends = "last", 
                           type = "closed"), 
             curvature = -0.65) + # Z -> Y
  
  #----------------------
  # T-M-Y Section (Treatment, Mediator, and Outcome)
  #----------------------

  # Draw rectangular boxes around T, M, and Y nodes

  # T: Treatment variable, with a box drawn around its coordinates
  geom_rect(aes(xmin = 1.25, ymin = 1.75, 
                xmax = 1.75, ymax = 2.25), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # M: Mediator variable, represented with a box at its defined coordinates
  geom_rect(aes(xmin = 3.25, ymin = 3.25, 
                xmax = 3.75, ymax = 3.75), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # Y: Outcome variable, drawn as a box at the Y-node coordinates
  geom_rect(aes(xmin = 5.25, ymin = 1.75, 
                xmax = 5.75, ymax = 2.25), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) + 
  
  # Draw the causal paths between T, M, and Y
  
  # Solid arrow from T (Treatment) to M (Mediator)
  geom_segment(aes(x = 1.8, y = 2.1, 
                   xend = 3.2, yend = 3.3), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  # Solid arrow from M (Mediator) to Y (Outcome)
  geom_segment(aes(x = 3.8, y = 3.3, 
                   xend = 5.2, yend = 2.1), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  # Direct solid arrow from T (Treatment) to Y (Outcome)
  geom_segment(aes(x = 1.8, y = 2, 
                   xend = 5.2, yend = 2), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  
  #----------------------
  # x1-x3 Section (Independent Variables)
  #----------------------
  
  # Draw rectangular boxes around x1, x2, and x3 variables

  # x1 variable box
  geom_rect(aes(xmin = 2.5, ymin = 5.25, 
                xmax = 3, ymax = 4.75), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # x2 variable box
  geom_rect(aes(xmin = 3.25, ymin = 5.25, 
                xmax = 3.75, ymax = 4.75), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # x3 variable box
  geom_rect(aes(xmin = 4, ymin = 5.25, 
                xmax = 4.5, ymax = 4.75), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # Draw a larger box encapsulating all three independent variables (x1-x3)
  geom_rect(aes(xmin = 2.25, ymin = 4.6, 
                xmax = 4.75, ymax = 5.4), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  
  # Add arrows showing causal paths from x1-x3 to T, M, and Y
  
  # Arrow from x1-x3 to T (Treatment)
  geom_segment(aes(x = 3.4, y = 4.55, 
                   xend = 1.5, yend = 2.3), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  # Arrow from x1-x3 to M (Mediator)
  geom_segment(aes(x = 3.5, y = 4.55, 
                   xend = 3.5, yend = 3.8), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  # Arrow from x1-x3 to Y (Outcome)
  geom_segment(aes(x = 3.6, y = 4.55, 
                   xend = 5.5, yend = 2.3), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  
  # Apply DAG-specific theme
  theme_dag()

# Display the DAG plot
p1



# Set Parameters & Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = c(70, 100),
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
dir.create(path = "Output/S1_Supp3_Results")
path <- "Output/S1_Supp3_Results/2025-03-30_1000-reps"
dir.create(path = path)
dir.create(path = paste0(path, "/Data"))
dir.create(path = paste0(path, "/Tables"))
dir.create(path = paste0(path, "/Figures"))
retrieval_path <- "Output/S1_Supp3_Simulation-Output/2025-03-30_1000-reps"


# Import data  ------------------------------------------------------------

# List all files matching the pattern.
file_list <- list.files(path = retrieval_path,
                        pattern = "^S1_Supp3_Condition-[0-9]+-Overall_Estimates_.*\\.rds$",
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

# ══════════════════════════════
#    TEMPORARY Add number of clusters to data  
# ══════════════════════════════
# sim1_data <- sim1_data |> 
#   left_join(mutate(cond[, c("num_clust", "num_x")], n = 1:nrow(cond)), by = c("conditionNum" = "n"))

# cond
# mutate(cond, n = 1:nrow(cond))


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
          file = paste0(path, "/Data/", "S1_Supp3_Performance-Measures.rds"))

write_rds(sim1_data, 
          file = paste0(path, "/Data/", "S1_Supp3_Simulation-Data.rds"))


# TNDE Relative Bias Table ------------------------------------------------

# ══════════════════════════════
#    num_clust = 70 
# ══════════════════════════════

# Pivot wide 
relBias_Tbl70 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 70, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "TNDE_relBias",
                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names 
colnames(relBias_Tbl70) <-
  stringr::str_remove(string = colnames(relBias_Tbl70),
                      pattern = "_relBias")

# Relative Bias TNDE Table   
relBias_TNDE_Table70 <- relBias_Tbl70 %>% 
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
as_hux(relBias_TNDE_Table70) %>% 
  set_number_format(row = 1:nrow(relBias_TNDE_Table70) + 1, col = -c(1:2), value = 3)

# Save csv 
write.csv(relBias_TNDE_Table70,
          file = paste0(path, "/Tables/", "TNDE-Relative-Bias_num_clust-70.csv"), row.names = FALSE)


# ══════════════════════════════
#    num_clust = 100 
# ══════════════════════════════

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

# ══════════════════════════════
#    num_clust = 70 
# ══════════════════════════════

# Pivot wide
relBias_Tbl70 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 70, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "TNDE_relBias",
                      "PNIE_relBias")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("TNDE_relBias", "PNIE_relBias")
)

# drop common peice of column names
colnames(relBias_Tbl70) <-
  stringr::str_remove(string = colnames(relBias_Tbl70),
                      pattern = "_relBias")

# Relative Bias PNIE Table   
relBias_PNIE_Table70 <- relBias_Tbl70 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_PNIE_Table70) %>%
  set_number_format(
    row = 1:nrow(relBias_PNIE_Table70) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(relBias_PNIE_Table70,
          file = paste0(path, "/Tables/", "PNIE-Relative-Bias_num_clust-70.csv"), row.names = FALSE)

# ══════════════════════════════
#    num_clust = 100 
# ══════════════════════════════

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

# ══════════════════════════════
#    num_clust = 70 
# ══════════════════════════════

# Pivot wide
rmse_Tbl70 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 70, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNIE_RMSE",
                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl70) <-
  stringr::str_remove(string = colnames(rmse_Tbl70),
                      pattern = "_RMSE")

# RMSE TNDE Table   
rmse_TNDE_Table70 <- rmse_Tbl70 %>%
  dplyr::select(analysisCond, starts_with("TNDE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_TNDE_Table70) %>%
  set_number_format(
    row = 1:nrow(rmse_TNDE_Table70) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(rmse_TNDE_Table70,
          file = paste0(path, "/Tables/", "TNDE-RMSE_num_clust-70.csv"), row.names = FALSE)

# ══════════════════════════════
#    num_clust = 100 
# ══════════════════════════════

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

# ══════════════════════════════
#    num_clust = 70 
# ══════════════════════════════

# Pivot wide
rmse_Tbl70 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 70, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNIE_RMSE",
                      "TNDE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNIE_RMSE", "TNDE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl70) <-
  stringr::str_remove(string = colnames(rmse_Tbl70),
                      pattern = "_RMSE")

# RMSE PNIE Table   
rmse_PNIE_Table70 <- rmse_Tbl70 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_PNIE_Table70) %>%
  set_number_format(
    row = 1:nrow(rmse_PNIE_Table70) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(rmse_PNIE_Table70,
          file = paste0(path, "/Tables/", "PNIE-RMSE_num_clust-70.csv"), row.names = FALSE)

# ══════════════════════════════
#    num_clust = 100 
# ══════════════════════════════

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

# ══════════════════════════════
#    num_clust = 70 
# ══════════════════════════════

# Pivot wide
cover_Tbl70 <- pivot_wider(
  perf_measure_DF[perf_measure_DF$num_clust == 70, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "coverage_PNIE",
                      "coverage_TNDE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("coverage_PNIE", "coverage_TNDE")
)

# drop common piece of column names
colnames(cover_Tbl70) <-
  stringr::str_remove(string = colnames(cover_Tbl70),
                      pattern = "coverage_")

# RMSE PNIE Table   
cover_PNIE_Table70 <- cover_Tbl70 %>%
  dplyr::select(analysisCond, starts_with("PNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Med", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(cover_PNIE_Table70) %>%
  set_number_format(
    row = 1:nrow(cover_PNIE_Table70) + 1,
    col = -c(1:2),
    value = 3
  )

# Save csv 
write.csv(cover_PNIE_Table70,
          file = paste0(path, "/Tables/", "PNIE-Coverage_num_clust-70.csv"), row.names = FALSE)

# ══════════════════════════════
#    num_clust = 100 
# ══════════════════════════════

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

# ══════════════════════════════
#    Boxplot (num_clust = 70) 
# ══════════════════════════════

# pdf("Output/S1_Results/Figures/PNIE-Absolute-Relative-Bias-Boxplot_num_clust-70.pdf")
readRDS(file = paste0(path, "/Data/S1_Supp3_Simulation-Data.rds")) |> 
  filter(num_clust == 70) |> 
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
# dev.off()
# Save visual 
ggsave(filename = paste0(path, "/Figures/PNIE-Absolute-Relative-Bias-Boxplot_num_clust-70.png"), 
       plot = last_plot())

# ══════════════════════════════
#    Boxplot (num_clust = 100) 
# ══════════════════════════════

readRDS(file = paste0(path, "/Data/S1_Supp3_Simulation-Data.rds")) |> 
  filter(num_clust == 100) |> 
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

# ══════════════════════════════
#    Boxplot (num_clust = 70) 
# ══════════════════════════════

readRDS(file = paste0(path, "/Data/S1_Supp3_Simulation-Data.rds")) |> 
  filter(num_clust == 70) |> 
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
ggsave(filename = paste0(path, "/Figures/PNIE-Relative-Bias-Boxplot_num_clust-70.png"), plot = last_plot())

# ══════════════════════════════
#    Boxplot (num_clust = 100) 
# ══════════════════════════════

readRDS(file = paste0(path, "/Data/S1_Supp3_Simulation-Data.rds")) |> 
  filter(num_clust == 100) |> 
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
         x.sec = guide_none("Number of Clusters and Residual ICC"))
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
  facet_grid(Out ~ num_clust + ICC) +
  gglayer_theme_line +
  gglayer_labs
# Save plot 
ggsave(filename = paste0(path, "/Figures/", 
                         "PNIE-Coverage-Lineplot.png"), 
       plot = last_plot())


## Examining SL Outcome further --------------------------------------------

readRDS(file = paste0(path, "/Data/S1_SuppC_Simulation-Data.rds")) |> 
  filter(ICC == 0.2, clust_size == 40) |> 
  filter(PS %in% "SL", medmodel %in% "SL", outModel %in% "SL") |> 
  mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
  ggplot(aes(x = rep, y = PNIE_est)) +
  geom_hline(yintercept = PNIE) +
  geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL), width = 0.2, color = "darkgray") +
  geom_point(size = 2, color = "blue") 
# coverage is 0

## Examining FE Outcome further --------------------------------------------

# DF for coverage (used for text)
coverDF <- perf_measure_DF |> 
  separate(analysisCond, into = c("PS", "Med", "Out"), sep = "_") |> 
  filter(ICC == 0.2, Out == "FE") |> 
  mutate(coverage_PNIE = format(coverage_PNIE, digits = 3))

# FE
coverFEDF <- readRDS(file = paste0(path, "/Data/S1_SuppC_Simulation-Data.rds")) |> 
    filter(ICC == 0.2) |> #, clust_size == 40) |> 
    filter(PS %in% c("FE"), medmodel %in% "FE", outModel %in% "FE") |> 
    mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
    mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
           if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
           error_color = ifelse(if_cover_PNIE == TRUE, "#BF5700", "darkgray")) |>
    arrange(PNIE_est) |> 
    group_by(clust_size) |> 
    mutate(order = row_number())
coverPlotFE <- coverFEDF |>
    ggplot(aes(x = order, y = PNIE_est)) +
    geom_hline(yintercept = PNIE) +
    # Map error_color to the error bars
    geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
    # Fix point color to match the “matched” color
    geom_point(size = 0.5, color = "#BF5700") +
    scale_color_identity() +
    scale_x_continuous(breaks = seq(0, max(coverFEDF$order, na.rm = TRUE), by = 50)) +
    # theme_minimal() +
    gglayer_theme +
    theme(axis.text.x.bottom = element_blank(), 
          axis.ticks = element_blank()) +
    labs(x = "", 
         y = "") +
    facet_grid(PS ~ clust_size) +
    geom_text(
      data = coverDF[coverDF$PS == "FE", ], aes(x = Inf, y = Inf, 
                                                label = paste0("coverage = ", coverage_PNIE, "\n",
                                                               "rel bias = ", format(round(PNIE_relBias, 5), scientific = FALSE))), 
      hjust = 2, vjust = 1.5, size = 3, color = "black") +
  guides(x.sec = guide_none("Cluster Size"))

# RE
coverREDF <- readRDS(file = paste0(path, "/Data/S1_SuppC_Simulation-Data.rds")) |> 
  filter(ICC == 0.2) |> #, clust_size == 40) |> 
  filter(PS %in% c("RE"), medmodel %in% "FE", outModel %in% "FE") |> 
  mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
  mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
         if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
         error_color = ifelse(if_cover_PNIE == TRUE, "#A6CD57", "darkgray")) |>
  arrange(PNIE_est) |> 
  group_by(clust_size) |> 
  mutate(order = row_number())
coverPlotRE <- coverREDF |> 
  ggplot(aes(x = order, y = PNIE_est)) +
  geom_hline(yintercept = PNIE) +
  # Map error_color to the error bars
  geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
  # Fix point color to match the “matched” color
  geom_point(size = 0.5, color = "#A6CD57") +
  scale_color_identity() +
  scale_x_continuous(breaks = seq(0, max(coverREDF$order, na.rm = TRUE), by = 50)) +
  # theme_minimal() +
  gglayer_theme +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks = element_blank()) +
  labs(x = "", 
       y = "NIE") +
  facet_grid(PS ~ clust_size) +
  geom_text(
    data = coverDF[coverDF$PS == "RE", ], aes(x = Inf, y = Inf, 
                                              label = paste0("coverage = ", coverage_PNIE, "\n",
                                                             "rel bias = ", format(round(PNIE_relBias, 5), scientific = FALSE))), 
    hjust = 2, vjust = 1.5, size = 3, color = "black", inherit.aes = TRUE)  +
  guides(y.sec = guide_none("PS Model"))

# SL
coverSLDF <- readRDS(file = paste0(path, "/Data/S1_SuppC_Simulation-Data.rds")) |> 
  filter(ICC == 0.2) |> #, clust_size == 40) |> 
  filter(PS %in% c("SL"), medmodel %in% "FE", outModel %in% "FE") |> 
  mutate(across(c(rep, PNIE_est, NIE_LCL, NIE_UCL), ~ as.numeric(.))) |> 
  mutate(if_cover_PNIE = (NIE_LCL < PNIE) & (NIE_UCL > PNIE),
         if_cover_TNDE = (NDE_LCL < TNDE) & (NDE_UCL > TNDE), 
         error_color = ifelse(if_cover_PNIE == TRUE, "#333F48", "darkgray")) |>
  arrange(PNIE_est) |> 
  group_by(clust_size) |> 
  mutate(order = row_number())
coverPlotSL <- coverSLDF |> 
  ggplot(aes(x = order, y = PNIE_est)) +
  geom_hline(yintercept = PNIE) +
  # Map error_color to the error bars
  geom_errorbar(aes(ymin = NIE_LCL, ymax = NIE_UCL, color = error_color), width = 0.2) +
  # Fix point color to match the “matched” color
  geom_point(size = 0.5, color = "#333F48") +
  scale_color_identity() +
  scale_x_continuous(breaks = seq(0, max(coverSLDF$order, na.rm = TRUE), by = 50)) +
  # theme_minimal() +
  gglayer_theme +
  # theme(axis.text.x.bottom = element_blank()) +
  labs(x = "", 
       y = "") +
  facet_grid(PS ~ clust_size) +
  geom_text(
    data = coverDF[coverDF$PS == "SL", ], aes(x = Inf, y = Inf, 
                                              label = paste0("coverage = ", coverage_PNIE, "\n",
                                                             "rel bias = ", format(round(PNIE_relBias, 5), scientific = FALSE))), 
    hjust = 2, vjust = 1.5, size = 3, color = "black", inherit.aes = TRUE) 

# Plot visual 
coverPlotFE / coverPlotRE / coverPlotSL +
  patchwork::plot_annotation(title = "Coverage rate across replications by PS & cluster size for ICC = 0.2 and med/outcome models = FE") 
# Save plot 
ggsave(filename = paste0(path, "/Figures/", 
                         "S1_SuppC_coverage-per-rep-by-PS-model-and-cluster-size-for-icc-0.2-and-outcome-FE.png"), 
       plot = last_plot())


