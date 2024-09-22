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
# Last Updated: 09/13/2024
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



# Simulation 2 Explained Variance Table -----------------------------------

# Set simulation parameters
num_x = 6

# Define coefficients for the effects of X variables
a_x <- 0.2  # Effect of X on treatment (trt)
b_x <- 0.8  # Effect of X on mediator (med)
c_x <- 0.8  # Effect of X on outcome

# Define coefficients for the effects of Z variable
a_z <- 0.5  # Effect of Z on treatment (trt)
c_z <- 1.4  # Effect of Z on outcome

# Define the effects of treatment and mediator on the outcome
treat_m <- 1        # Effect of treatment (trt) on mediator (med)
treat_y <- 1.3      # Effect of treatment (trt) on outcome
med_y <- 1          # Effect of mediator (med) on outcome
treat_med_y <- 1.15 # Effect of treatment-mediator interaction on outcome 

# Create an empty data frame to store explained variance (R^2) values for different ICC levels
design_values_DF <- data.frame(cbind(
  var_expl = c("trt_x", "med_x", "out_x", 
               "trt_z", "med_z", "out_z", 
               "med_trt", "out_trt", 
               "out_med", 
               "out_trt_med"), 
  icc_0.05 = NA, 
  icc_0.2 = NA, 
  icc_0.5 = NA))

# Loop through ICC values and calculate the explained variance (R^2) for each variable
for (icc in c(0.05, 0.2, 0.5)) {
  # Variance explained by X variables on treatment (trt_x), mediator (med_x), and outcome (out_x)
  design_values_DF[design_values_DF$var_expl == "trt_x", paste0("icc_", icc)] <-
    ((a_x^2)*num_x * (1 - icc)) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "med_x", paste0("icc_", icc)] <-
    ((b_x^2)*(num_x/2) * (1 - icc)) / ((b_x^2)*(num_x/2) + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  design_values_DF[design_values_DF$var_expl == "out_x", paste0("icc_", icc)] <-
    (c_x^2)*num_x * (1 - icc) / (treat_y^2 + (c_x^2)*num_x + c_z^2 + med_y^2 + treat_med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by Z variables on treatment (trt_z) and outcome (out_z)
  design_values_DF[design_values_DF$var_expl == "trt_z", paste0("icc_", icc)] <-
    (a_z^2) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "med_z", paste0("icc_", icc)] <-
    0 # No effect of Z on mediator (med) in this model
  
  design_values_DF[design_values_DF$var_expl == "out_z", paste0("icc_", icc)] <-
    c_z^2 / (treat_y^2 + (c_x^2)*num_x + c_z^2 + med_y^2 + treat_med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by the treatment variable on mediator (med_trt) and outcome (out_trt)
  design_values_DF[design_values_DF$var_expl == "med_trt", paste0("icc_", icc)] <-
    treat_m^2 * (1 - icc) / ((b_x^2)*(num_x/2) + treat_m^2 + (pi^2/3)/4 + (pi^2/3))
  
  design_values_DF[design_values_DF$var_expl == "out_trt", paste0("icc_", icc)] <-
    treat_y^2 * (1 - icc) / (treat_y^2 + (c_x^2)*num_x + c_z^2 + med_y^2 + treat_med_y^2 + (pi^2/3)/4 + (pi^2/3))
  
  # Variance explained by the mediator variable on outcome (out_med)
  design_values_DF[design_values_DF$var_expl == "out_med", paste0("icc_", icc)] <-
    med_y^2 * (1 - icc) / (treat_y^2 + (c_x^2)*num_x + c_z^2 + med_y^2 + treat_med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
  
  # Variance explained by the interaction between treatment and mediator on outcome (out_trt_med)
  design_values_DF[design_values_DF$var_expl == "out_trt_med", paste0("icc_", icc)] <- 
    treat_med_y^2 * (1 - icc) / (treat_y^2 + (c_x^2)*num_x + c_z^2 + med_y^2 + treat_med_y^2 + (pi^2/3)/4 + (pi^2/3))
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
      "Total variance in the treatment explained by all X variables (x1-x6)",
      "Total variance in the mediator explained by x1-x3 variables",
      "Total variance in the outcome explained by all X variables (x1-x6)",
      
      "Total variance in the mediator explained by the treatment",
      "Total variance in the outcome explained by the treatment",
      "Total variance in the outcome explained by the mediator",
      "Total variance in the outcome explained by the treatment-mediator interaction"
    ),
    "0.05" = design_values_DF$icc_0.05,
    L_0.2 = design_values_DF$icc_0.2,
    L_0.5 = design_values_DF$icc_0.5
  )
)

explained_var_table



# Simulation 2 Directed Acyclic Graph (DAG) Visualization -----------------

# Define the DAG structure and coordinates for each variable
dag2 <- dagify(
  Z ~ 1,              # Latent variable Z (no parent)
  x1 ~ 1, x2 ~ 1, x3 ~ 1,  # Independent variables x1, x2, x3 (no parents)
  x4 ~ 1, x5 ~ 1, x6 ~ 1,  # Independent variables x4, x5, x6 (no parents)
  T ~ 1,              # Treatment variable T (no parent)
  M ~ 1,              # Mediator variable M (no parent)
  Y ~ 1,              # Outcome variable Y (no parent)
  exposure = "T",     # Define T as the exposure
  outcome = "Y",      # Define Y as the outcome
  latent = "Z",       # Define Z as a latent variable
  coords = list(      # Coordinates for plotting the DAG nodes
    x = c(T = 1.5, M = 3.5, Y = 5.5, Z = 1,
          x1 = 2.75, x2 = 3.5, x3 = 4.25,
          x4 = 2.75, x5 = 3.5, x6 = 4.25),
    y = c(T = 2, M = 3.5, Y = 2, Z = 6.5,
          x1 = 5, x2 = 5, x3 = 5,
          x4 = 0, x5 = 0, x6 = 0)
  )
)


# Create DAG viz 
p2 <- ggdag_classic(dag2, size = 4) +
  
  # Set plot limits to ensure all elements fit well within the plot area
  xlim(0, 6) +
  ylim(-0.4, 6.8) +
  coord_fixed() +  # Maintain fixed aspect ratio
  
  #----------------------
  # Z Section (Latent Variable Z)
  #----------------------
  
  # Add Z (latent variable) as a circle at specified coordinates
  ggforce::stat_circle(aes(x0 = 1, y0 = 6.5, r = 0.3)) +
  
  # Add dashed arrows from Z to T and Y
  geom_curve(x = 0.65, y = 6.3, 
             xend = 1.25, yend = 2.3, 
             linetype = "dashed", 
             arrow = arrow(length = unit(0.2, "cm"), 
                           ends = "last", 
                           type = "closed"), 
             curvature = 0.4) + # Z -> T
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
  
  #----------------------
  # x4-x6 Section (Independent Variables)
  #----------------------
  
  # Draw rectangular boxes around x4, x5, and x6 variables
  
  # x4 variable box
  geom_rect(aes(xmin = 2.5, ymin = -0.25, 
                xmax = 3, ymax = 0.25), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # x5 variable box
  geom_rect(aes(xmin = 3.25, ymin = -0.25, 
                xmax = 3.75, ymax = 0.25), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # x6 variable box
  geom_rect(aes(xmin = 4, ymin = -0.25, 
                xmax = 4.5, ymax = 0.25), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  # Draw a larger box encapsulating all three independent variables (x4-x6)
  geom_rect(aes(xmin = 2.3, ymin = -0.4, 
                xmax = 4.7, ymax = 0.4), 
            fill = "transparent", color = "black", linewidth = 0.4, size = 1) +
  
  # Add arrows showing causal paths from x4-x6 to T and Y
  
  # Arrow from x4-x6 to T (Treatment)
  geom_segment(aes(x = 3.4, y = 0.45, 
                   xend = 1.5, yend = 1.7), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  # Arrow from x4-x6 to Y (Outcome)
  geom_segment(aes(x = 3.6, y = 0.45, 
                   xend = 5.5, yend = 1.7), 
               linewidth = 0.4, linetype = "solid", color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
  
  # Apply DAG-specific theme
  theme_dag()

# Display the DAG plot
p2



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
TNIE <- treat_m * med_y + treat_med_y

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
            
            # PNDE_MSE = (mean(PNDE_est) - PNDE)^2 + (sd(PNDE_est)^2), 
            # TNIE_MSE = (mean(TNIE_est) - TNIE)^2 + (sd(TNIE_est)^2), 
            
            PNDE_RMSE = sqrt((mean(PNDE_est) - PNDE)^2 + (sd(PNDE_est)^2)), 
            TNIE_RMSE = sqrt((mean(TNIE_est) - TNIE)^2 + (sd(TNIE_est)^2))
  ) 



# Export Performance Measures & Simulation Data ---------------------------------------------

write_rds(perf_measure_DF, 
          file = paste0("Output/S2B_Results/Data/", "S2B_Performance-Measures.rds"))

write_rds(sim2B_data, 
          file = paste0("Output/S2B_Results/Data/", "S2B_Simulation-Data.rds"))



# PNDE Relative Bias Table ------------------------------------------------

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

# Relative Bias PNDE Table   
relBias_PNDE_Table <- relBias_Tbl %>% 
  select(analysisCond, starts_with("PNDE_")) %>% 
  separate(col = c("analysisCond"), into = c("PS Model", "Mediator/Outcome Model"), sep = "_") %>% 
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(relBias_PNDE_Table) %>% 
  set_number_format(row = 1:nrow(relBias_PNDE_Table) + 1, col = -c(1:2), value = 3)



# TNIE Relative Bias Table ------------------------------------------------

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

# Relative Bias TNIE Table   
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



# PNDE Relative Bias Table ------------------------------------------------

# Pivot wide
rmse_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNDE_RMSE",
                      "TNIE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNDE_RMSE", "TNIE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl) <-
  stringr::str_remove(string = colnames(rmse_Tbl),
                      pattern = "_RMSE")

# RMSE PNDE Table   
rmse_PNDE_Table <- rmse_Tbl %>%
  select(analysisCond, starts_with("PNDE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_PNDE_Table) %>%
  set_number_format(
    row = 1:nrow(rmse_PNDE_Table) + 1,
    col = -c(1:2),
    value = 3
  )



# TNIE RMSE Table ------------------------------------------------

# Pivot wide
rmse_Tbl <- pivot_wider(
  perf_measure_DF[, c("ICC",
                      "clust_size",
                      "analysisCond",
                      "PNDE_RMSE",
                      "TNIE_RMSE")],
  names_from = c(ICC, clust_size),
  names_sep = "_",
  values_from = c("PNDE_RMSE", "TNIE_RMSE")
)

# drop common piece of column names
colnames(rmse_Tbl) <-
  stringr::str_remove(string = colnames(rmse_Tbl),
                      pattern = "_RMSE")

# RMSE TNIE Table   
rmse_TNIE_Table <- rmse_Tbl %>%
  select(analysisCond, starts_with("TNIE_")) %>%
  separate(
    col = c("analysisCond"),
    into = c("PS Model", "Mediator/Outcome Model"),
    sep = "_"
  ) %>%
  arrange(`Mediator/Outcome Model`)

# Display values (Copy-paste values into generate table website to obtain latex code)
as_hux(rmse_TNIE_Table) %>%
  set_number_format(
    row = 1:nrow(rmse_TNIE_Table) + 1,
    col = -c(1:2),
    value = 3
  )


# PNDE Relative Bias Visual -----------------------------------------------

# set direct & indirect effects 
treat_m <- 1 # trt on med   
treat_y <- 1.3 # trt on outcome
med_y <- 1 # med on outcome  
treat_med_y <- 1.15 # trt-med interaction on outcome 
PNDE <- treat_y
TNIE <- treat_m * med_y + treat_med_y

# # set values 
# treat_m <- 1 # trt on med   
# treat_y <- 1.3 # trt on outcome
# med_y <- 1 # med on outcome  
# TNDE <- treat_y
# PNIE <- treat_m * med_y

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

pdf("Output/S2B_Results/Figures/PNDE-Relative-Bias-Boxplot.pdf")
readRDS(file = "Output/S2B_Results/Data/S2B_Simulation-Data.rds") |>
  # rename PS models 
  mutate(`PS Model` = ifelse(PS == "FE", "Fixed-Effect", 
                             ifelse(PS == "RE", "Random-Effect", 
                                    ifelse(PS == "SL", "Single-Level", 
                                           "ERROR")))) %>% 
  
  # visual 
  ggplot(aes(x = abs((PNDE_est / PNDE) - 1), y = outModel, fill = `PS Model`)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.10, color = "red", alpha = 0.4) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, 
               linewidth = 0.25, 
               outlier.size = 0.5, 
               outlier.alpha = 0.5) +
  facet_grid(clust_size ~ ICC) +
  gglayer_labs +
  gglayer_theme 

dev.off()
# Save visual 
ggsave(filename = "Output/S2B_Results/Figures/PNDE-Relative-Bias-Boxplot.png", plot = last_plot())






