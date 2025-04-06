# This is old bootstrapping code that was in Empirical-Applicaiton_Analysis.R 



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


source("Application/Functions/bootstrap_ci_paral_2.R")
source("Application/Functions/bootstrap_ci_re_paral_2.R")
source("Application/Functions/bootstrap_ci_re_mean_paral.R")
# source("Application/Functions/monteCarloCI.R")
# source("Application/Functions/monteCarloCIb.R")
source("Application/Functions/bootstrapCI.R")
source("Application/Functions/bootstrapCIb.R")





# Bootstrap CI TESTING  ---------------------------------------------------

# WARNING: Both RE and RE-Mean calculations require extended processing times.


## PNDE & TNIE -------------------------------------------------------------
# This subsection focuses on PNDE (Pure Natural Direct Effect) and TNIE (Total Natural Indirect Effect)
# effects across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_TNIE <- bootstrapCIb(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify PS model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "TNIE"      # Effect type: TNIE
)
saveRDS(slsl_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/slsl_ci_TNIE.rds")
rm(slsl_ci_TNIE)

# FE PS Model
fesl_ci_TNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fesl_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/fesl_ci_TNIE.rds")
rm(fesl_ci_TNIE)

# RE PS Model
resl_ci_TNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(resl_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/resl_ci_TNIE.rds")
rm(resl_ci_TNIE)


### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_TNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(slfe_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/slfe_ci_TNIE.rds")
rm(slfe_ci_TNIE)

# FE PS Model
fefe_ci_TNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fefe_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/fefe_ci_TNIE.rds")
rm(fefe_ci_TNIE)

# RE PS Model
refe_ci_TNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(refe_ci_TNIE, file = "Application/Output/Temp-bootstrap-CIs/refe_ci_TNIE.rds")
rm(refe_ci_TNIE)


## TNDE & PNIE -------------------------------------------------------------
# This subsection focuses on TNDE (Total Natural Direct Effect) and PNIE (Pure Natural Indirect Effect) across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_PNIE <- bootstrapCIb(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify mediator & outcome model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "PNIE"      # Effect type: PNIE
)
saveRDS(slsl_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/slsl_ci_PNIE.rds")
rm(slsl_ci_PNIE)

# FE PS Model
fesl_ci_PNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fesl_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/fesl_ci_PNIE.rds")
rm(fesl_ci_PNIE)

# RE PS Model
resl_ci_PNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(resl_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/resl_ci_PNIE.rds")
rm(resl_ci_PNIE)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_PNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",   
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(slfe_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/slfe_ci_PNIE.rds")
rm(slfe_ci_PNIE)

# FE PS Model
fefe_ci_PNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fefe_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/fefe_ci_PNIE.rds")
rm(fefe_ci_PNIE)

# RE PS Model
refe_ci_PNIE <- bootstrapCIb(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(refe_ci_PNIE, file = "Application/Output/Temp-bootstrap-CIs/refe_ci_PNIE.rds")
rm(refe_ci_PNIE)





## Display & Save bootstrap CI Estimates -------------------------------------------------------



###### Join CI & point estimates for SL & FE med/outcome -------------------
# This block imports data for SL and FE models, calculates CIs, and combines results.

# Define the list of file names
list_names <- c("slsl_ci_TNIE", "slsl_ci_PNIE", 
                "fesl_ci_TNIE", "fesl_ci_PNIE",
                "resl_ci_TNIE", "resl_ci_PNIE", 
                "slfe_ci_TNIE", "slfe_ci_PNIE", 
                "fefe_ci_TNIE", "fefe_ci_PNIE",
                "refe_ci_TNIE", "refe_ci_PNIE")

# Import RDS files and store them in a list of lists
lists <- lapply(list_names, function(name) {
  readRDS(paste0("Application/Output/Temp-bootstrap-CIs/", name, ".rds"))
})

# Name the elements of the list for easier access
names(lists) <- list_names

# Extract 'indirect_ci' and store in dataframe
indirect_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$indirect_ci
  data.frame(
    list_name = list_names[i],
    indirect_ci_LL = ci_values[1],  # Lower bound
    indirect_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Extract 'direct_ci' and store in dataframe
direct_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$direct_ci
  data.frame(
    list_name = list_names[i],
    direct_ci_LL = ci_values[1],  # Lower bound
    direct_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Combine both 'indirect' and 'direct' CIs into a single dataframe
combined_df <- merge(indirect_df, direct_df, by = "list_name")

# Modify dataframe: extract effect and cond labels, and reshape to wide format
combined_df <- combined_df %>%
  mutate(effect = substr(list_name, nchar(list_name) - 3, nchar(list_name))) %>%  # Extract effect label
  mutate(list_name = substr(list_name, 1, 4)) %>%  # Extract cond label
  pivot_longer(cols = c(indirect_ci_LL, indirect_ci_UL, direct_ci_LL, direct_ci_UL), 
               names_to = "variable", values_to = "value") %>%
  unite("variable", effect, variable, sep = "_") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  as.data.frame()

# Update column names for clarity
colnames(combined_df) <- c("cond", 
                           "PNIE_LL", "PNIE_UL", "TNDE_LL", "TNDE_UL", 
                           "TNIE_LL", "TNIE_UL", "PNDE_LL", "PNDE_UL")

# Merge the SL & FE results with the existing results dataframe
results_DF_noRE <- merge(results_DF[1:6, ], combined_df, by = "cond")




# ═══════════════════
#    NOTE WE ARE SKIPPING RE & RE-MEAN MED/OUT MODELS FOR NOW 
# ═══════════════════



# Add labels for PS model & Mediator/Outcome model
results_DF_noRE <- results_DF_noRE %>%
  mutate(
    PS = case_when(
      startsWith(cond, "fe") ~ "Fixed-Effect",
      startsWith(cond, "sl") ~ "Single-Level",
      startsWith(cond, "re") ~ "Random-Effect",
      TRUE ~ NA_character_  # Default case
    ),
    Model = case_when(
      endsWith(cond, "fe") ~ "Fixed-Effect",
      endsWith(cond, "sl") ~ "Single-Level",
      endsWith(cond, "re") ~ "Random-Effect",
      endsWith(cond, "re_cm") ~ "Random-Effect with Cluster Means",
      TRUE ~ NA_character_  # Default case
    )
  )

# Display the final dataframe with PS model & Mediator/Outcome labels
print(results_DF_noRE)

# Save the final dataframe of results
write_rds(results_DF_noRE, file = "Application/Output/Estimates/Effect-Estimates_bootstrap-CIs.rds")

















# Bootstrap Confidence Intervals (CI) -------------------------------------
# This section conducts the bootstrapping for confidence intervals for each effect 
# using each PS model (SL, FE, & RE) and mediator/outcome model (SL, FE, RE, & RE-Mean). 
# Results are saved and convergence statistics are printed.

# WARNING: Both RE and RE-Mean calculations require extended processing times.

## TNDE & PNIE -------------------------------------------------------------
# This subsection focuses on TNDE (Total Natural Direct Effect) and PNIE (Pure Natural Indirect Effect) across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify mediator & outcome model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "PNIE"      # Effect type: PNIE
)
saveRDS(slsl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slsl_ci_PNIE.rds")
rm(slsl_ci_PNIE)

# FE PS Model
fesl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fesl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fesl_ci_PNIE.rds")
rm(fesl_ci_PNIE)

# RE PS Model
resl_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(resl_ci_PNIE, file = "Application/Output/Bootstrap_Temp/resl_ci_PNIE.rds")
rm(resl_ci_PNIE)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",   
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(slfe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slfe_ci_PNIE.rds")
rm(slfe_ci_PNIE)

# FE PS Model
fefe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(fefe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fefe_ci_PNIE.rds")
rm(fefe_ci_PNIE)

# RE PS Model
refe_ci_PNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "PNIE"
)
saveRDS(refe_ci_PNIE, file = "Application/Output/Bootstrap_Temp/refe_ci_PNIE.rds")
rm(refe_ci_PNIE)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
execution_time <- system.time({ # Track computation time 
  slre_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,  
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(slre_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slre_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_ci_PNIE$mediator_converged_count,
    " (", (slre_ci_PNIE$mediator_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_ci_PNIE$outcome_converged_count,
    " (", (slre_ci_PNIE$outcome_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_ci_PNIE$both_converged_count,
    " (", (slre_ci_PNIE$both_converged_count / length(slre_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(slre_ci_PNIE)
# Elapsed time: 6652.806 seconds ( 111 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# FE PS Model
execution_time <- system.time({ 
  fere_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(fere_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fere_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_ci_PNIE$mediator_converged_count,
    " (", (fere_ci_PNIE$mediator_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_ci_PNIE$outcome_converged_count,
    " (", (fere_ci_PNIE$outcome_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_ci_PNIE$both_converged_count,
    " (", (fere_ci_PNIE$both_converged_count / length(fere_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(fere_ci_PNIE)
# Elapsed time: 7617.093 seconds ( 127 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

# RE PS Model
execution_time <- system.time({ 
  rere_ci_PNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(rere_ci_PNIE, file = "Application/Output/Bootstrap_Temp/rere_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_ci_PNIE$mediator_converged_count,
    " (", (rere_ci_PNIE$mediator_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_ci_PNIE$outcome_converged_count,
    " (", (rere_ci_PNIE$outcome_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_ci_PNIE$both_converged_count,
    " (", (rere_ci_PNIE$both_converged_count / length(rere_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(rere_ci_PNIE)
# Elapsed time: 10816.36 seconds ( 180 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1692  ( 99.52941 %)
# Number of iterations with both models converged:  1231  ( 72.41176 %)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
execution_time <- system.time({ 
  slre_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(slre_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/slre_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_cm_ci_PNIE$mediator_converged_count,
    " (", (slre_cm_ci_PNIE$mediator_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_cm_ci_PNIE$outcome_converged_count,
    " (", (slre_cm_ci_PNIE$outcome_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_cm_ci_PNIE$both_converged_count,
    " (", (slre_cm_ci_PNIE$both_converged_count / length(slre_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(slre_cm_ci_PNIE)
# Elapsed time: 6615.656 seconds ( 110 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# FE PS Model
execution_time <- system.time({ 
  fere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(fere_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/fere_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_cm_ci_PNIE$mediator_converged_count,
    " (", (fere_cm_ci_PNIE$mediator_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_cm_ci_PNIE$outcome_converged_count,
    " (", (fere_cm_ci_PNIE$outcome_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_cm_ci_PNIE$both_converged_count,
    " (", (fere_cm_ci_PNIE$both_converged_count / length(fere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(fere_cm_ci_PNIE)
# Elapsed time: 6152.866 seconds ( 103 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

# RE PS Model
execution_time <- system.time({ 
  rere_cm_ci_PNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1750, 
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "PNIE"
  )
})
saveRDS(rere_cm_ci_PNIE, file = "Application/Output/Bootstrap_Temp/rere_cm_ci_PNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_cm_ci_PNIE$mediator_converged_count,
    " (", (rere_cm_ci_PNIE$mediator_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_cm_ci_PNIE$outcome_converged_count,
    " (", (rere_cm_ci_PNIE$outcome_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_cm_ci_PNIE$both_converged_count,
    " (", (rere_cm_ci_PNIE$both_converged_count / length(rere_cm_ci_PNIE$direct_effects)) * 100, "%)\n")
rm(rere_cm_ci_PNIE)
# Elapsed time: 5956.323 seconds ( 99 mins) 
# Number of converged mediator models:  1215  ( 69.42857 %)
# Number of converged outcome models:  1718  ( 98.17143 %)
# Number of iterations with both models converged:  1197  ( 68.4 %)

## PNDE & TNIE -------------------------------------------------------------
# This subsection focuses on PNDE (Pure Natural Direct Effect) and TNIE (Total Natural Indirect Effect)
# effects across different mediation/outcome models.

### Single-Level (SL) Mediation/Outcome Models -----------------------------
# SL PS Model
slsl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,        # Number of bootstrap iterations
  iptw = iptw_sl,           # IPTW weights from SL PS model
  data = data,              # Input data
  model = "SL",             # Specify PS model 
  cores = 6,                # Number of CPU cores for parallelization
  core_seeds = c(4561:4566),# Seeds for reproducibility
  effect_type = "TNIE"      # Effect type: TNIE
)
saveRDS(slsl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slsl_ci_TNIE.rds")
rm(slsl_ci_TNIE)

# FE PS Model
fesl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fesl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fesl_ci_TNIE.rds")
rm(fesl_ci_TNIE)

# RE PS Model
resl_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "SL",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(resl_ci_TNIE, file = "Application/Output/Bootstrap_Temp/resl_ci_TNIE.rds")
rm(resl_ci_TNIE)

### Fixed-Effect (FE) Mediation/Outcome Models -----------------------------
# SL PS Model
slfe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_sl,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(slfe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slfe_ci_TNIE.rds")
rm(slfe_ci_TNIE)

# FE PS Model
fefe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_fe,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(fefe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fefe_ci_TNIE.rds")
rm(fefe_ci_TNIE)

# RE PS Model
refe_ci_TNIE <- bootstrap_ci_paral_2(
  iterations = 1000,
  iptw = iptw_re,
  data = data,
  model = "FE",
  cores = 6,
  core_seeds = c(4561:4566),
  effect_type = "TNIE"
)
saveRDS(refe_ci_TNIE, file = "Application/Output/Bootstrap_Temp/refe_ci_TNIE.rds")
rm(refe_ci_TNIE)

### Random-Effect (RE) Mediation/Outcome Models ----------------------------
# SL PS Model
execution_time <- system.time({ 
  slre_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(slre_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slre_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_ci_TNIE$mediator_converged_count,
    " (", (slre_ci_TNIE$mediator_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_ci_TNIE$outcome_converged_count,
    " (", (slre_ci_TNIE$outcome_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_ci_TNIE$both_converged_count,
    " (", (slre_ci_TNIE$both_converged_count / length(slre_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(slre_ci_TNIE)
# Elapsed time: 6310.98 seconds ( 105 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

# FE PS Model
execution_time <- system.time({ 
  fere_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(fere_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fere_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_ci_TNIE$mediator_converged_count,
    " (", (fere_ci_TNIE$mediator_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_ci_TNIE$outcome_converged_count,
    " (", (fere_ci_TNIE$outcome_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_ci_TNIE$both_converged_count,
    " (", (fere_ci_TNIE$both_converged_count / length(fere_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(fere_ci_TNIE)
# Elapsed time: 6343.685 seconds ( 106 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

# RE PS Model
execution_time <- system.time({ 
  rere_ci_TNIE <- bootstrap_ci_re_paral_2(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(rere_ci_TNIE, file = "Application/Output/Bootstrap_Temp/rere_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_ci_TNIE$mediator_converged_count,
    " (", (rere_ci_TNIE$mediator_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_ci_TNIE$outcome_converged_count,
    " (", (rere_ci_TNIE$outcome_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_ci_TNIE$both_converged_count,
    " (", (rere_ci_TNIE$both_converged_count / length(rere_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(rere_ci_TNIE)
# Elapsed time: 10578.82 seconds ( 176 mins) 
# Number of converged mediator models:  1235  ( 72.64706 %)
# Number of converged outcome models:  1690  ( 99.41176 %)
# Number of iterations with both models converged:  1229  ( 72.29412 %)

### Random-Effect with Cluster Means (RE-Mean) Med/Out Models --------------
# SL PS Model
execution_time <- system.time({ 
  slre_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_sl,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(slre_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/slre_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", slre_cm_ci_TNIE$mediator_converged_count,
    " (", (slre_cm_ci_TNIE$mediator_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", slre_cm_ci_TNIE$outcome_converged_count,
    " (", (slre_cm_ci_TNIE$outcome_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", slre_cm_ci_TNIE$both_converged_count,
    " (", (slre_cm_ci_TNIE$both_converged_count / length(slre_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(slre_cm_ci_TNIE)
# Elapsed time: 9869.804 seconds ( 164 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD:
# Elapsed time: 4365.37 seconds ( 73 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)

# FE PS Model
execution_time <- system.time({ 
  fere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_fe,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(fere_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/fere_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", fere_cm_ci_TNIE$mediator_converged_count,
    " (", (fere_cm_ci_TNIE$mediator_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", fere_cm_ci_TNIE$outcome_converged_count,
    " (", (fere_cm_ci_TNIE$outcome_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", fere_cm_ci_TNIE$both_converged_count,
    " (", (fere_cm_ci_TNIE$both_converged_count / length(fere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(fere_cm_ci_TNIE)
# Elapsed time: 5830.462 seconds ( 97 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD: 
# Elapsed time: 3742.066 seconds ( 62 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)

# RE PS Model
execution_time <- system.time({ 
  rere_cm_ci_TNIE <- bootstrap_ci_re_mean_paral(
    iterations = 1700,
    iptw = iptw_re,
    data = data,
    cores = 6,
    core_seeds = c(4561:4566),
    effect_type = "TNIE"
  )
})
saveRDS(rere_cm_ci_TNIE, file = "Application/Output/Bootstrap_Temp/rere_cm_ci_TNIE.rds")

# Print elapsed time and convergence statistics
cat("Elapsed time:", execution_time["elapsed"], "seconds (", round(execution_time["elapsed"]/60), "mins) \n")
cat("Number of converged mediator models: ", rere_cm_ci_TNIE$mediator_converged_count,
    " (", (rere_cm_ci_TNIE$mediator_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of converged outcome models: ", rere_cm_ci_TNIE$outcome_converged_count,
    " (", (rere_cm_ci_TNIE$outcome_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
cat("Number of iterations with both models converged: ", rere_cm_ci_TNIE$both_converged_count,
    " (", (rere_cm_ci_TNIE$both_converged_count / length(rere_cm_ci_TNIE$direct_effects)) * 100, "%)\n")
rm(rere_cm_ci_TNIE)
# Elapsed time: 7174.466 seconds ( 120 mins) 
# Number of converged mediator models:  1181  ( 69.47059 %)
# Number of converged outcome models:  1664  ( 97.88235 %)
# Number of iterations with both models converged:  1161  ( 68.29412 %)

# OLD: 
# Elapsed time: 3800.283 seconds ( 63 mins) 
# Number of converged mediator models:  1219  ( 69.65714 %)
# Number of converged outcome models:  1722  ( 98.4 %)
# Number of iterations with both models converged:  1209  ( 69.08571 %)



# Store & Join Results ----------------------------------------------------
# This section handles the calculation and storage of confidence intervals (CIs). 

#### RE Mediator/Outcome Models CI -----------------------------------------
# This subsection focuses on obtaining CIs for RE & RE-Mean models and merging them with effect estimates.

###### Obtain 1,000 completed iterations -----------------------------------
# Function to get non-NA pairs (i.e., first 1,000 completed iterations)
get_non_na_pairs <- function(direct, indirect, n = 1000) {
  combined <- data.frame(direct = direct, indirect = indirect)  # Combine vectors into a dataframe
  combined <- na.omit(combined)  # Remove rows with any NA values
  return(head(combined, n))  # Return the first n rows (or all if less than n)
}

# Define the file paths for the RE models
files <- list(
  "Application/Output/Bootstrap_Temp/slre_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_ci_PNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_cm_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_cm_ci_PNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_cm_ci_PNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_ci_TNIE.rds",
  
  "Application/Output/Bootstrap_Temp/slre_cm_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/fere_cm_ci_TNIE.rds",
  "Application/Output/Bootstrap_Temp/rere_cm_ci_TNIE.rds"
)

# Function to process each file, apply get_non_na_pairs, and remove the original data
process_file <- function(file_path, n = 1000) {
  data <- readRDS(file_path)  # Load the RDS file
  result <- get_non_na_pairs(data$direct_effects, data$indirect_effects, n = n)  # Apply function
  rm(data)  # Remove the original data to free up space
  return(result)  # Return the processed data
}

# Process each file and store results in corresponding dataframes
slre_ci_PNIE_DF <- process_file(files[[1]])
fere_ci_PNIE_DF <- process_file(files[[2]])
rere_ci_PNIE_DF <- process_file(files[[3]])

slre_cm_ci_PNIE_DF <- process_file(files[[4]])
fere_cm_ci_PNIE_DF <- process_file(files[[5]])
rere_cm_ci_PNIE_DF <- process_file(files[[6]])

slre_ci_TNIE_DF <- process_file(files[[7]])
fere_ci_TNIE_DF <- process_file(files[[8]])
rere_ci_TNIE_DF <- process_file(files[[9]])

slre_cm_ci_TNIE_DF <- process_file(files[[10]])
fere_cm_ci_TNIE_DF <- process_file(files[[11]])
rere_cm_ci_TNIE_DF <- process_file(files[[12]])

# Remove 'files' vector to save space 
rm(files)

###### Store RE Med/Outcome Model CIs --------------------------------------
# This block stores the calculated CIs for RE models in a structured dataframe.

# Create an empty dataframe to store the results
results_DF_RE <- data.frame(
  cond = c("slre", "fere", "rere", 
           "slre_cm", "fere_cm", "rere_cm"),
  PNIE_LL = numeric(6),
  PNIE_UL = numeric(6),
  TNDE_LL = numeric(6),
  TNDE_UL = numeric(6),
  TNIE_LL = numeric(6),
  TNIE_UL = numeric(6),
  PNDE_LL = numeric(6),
  PNDE_UL = numeric(6),
  stringsAsFactors = FALSE
)

# Lists of processed dataframes
df_list_PNIE <- list(slre_ci_PNIE_DF, fere_ci_PNIE_DF, rere_ci_PNIE_DF, 
                     slre_cm_ci_PNIE_DF, fere_cm_ci_PNIE_DF, rere_cm_ci_PNIE_DF)
df_list_TNIE <- list(slre_ci_TNIE_DF, fere_ci_TNIE_DF, rere_ci_TNIE_DF, 
                     slre_cm_ci_TNIE_DF, fere_cm_ci_TNIE_DF, rere_cm_ci_TNIE_DF)

# Calculate CIs and fill the dataframe
for (i in 1:6) {
  results_DF_RE[i, c("PNIE_LL", "PNIE_UL")] <- quantile(df_list_PNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("TNDE_LL", "TNDE_UL")] <- quantile(df_list_PNIE[[i]]$direct, probs = c(0.025, 0.975))
  
  results_DF_RE[i, c("TNIE_LL", "TNIE_UL")] <- quantile(df_list_TNIE[[i]]$indirect, probs = c(0.025, 0.975))
  results_DF_RE[i, c("PNDE_LL", "PNDE_UL")] <- quantile(df_list_TNIE[[i]]$direct, probs = c(0.025, 0.975))
}

# Display the RE model results
results_DF_RE
#       cond    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL
# 1    slre -0.2792701 0.09994341 -0.5000102 -0.055682131 -0.2819150 0.10238444 -0.5027937 -0.053682580
# 2    fere -0.3426149 0.07165350 -0.5128617 -0.019704421 -0.3590603 0.07522094 -0.5137532 -0.020466528
# 3    rere -0.3106561 0.08389139 -0.5176295 -0.050033639 -0.3216423 0.09003098 -0.5162965 -0.049349202
# 4 slre_cm -0.3235201 0.08291188 -0.4695924  0.003501834 -0.3145743 0.08801792 -0.4681990  0.003553328
# 5 fere_cm -0.3494602 0.07261279 -0.5102075 -0.008612200 -0.3630953 0.07850820 -0.5110733 -0.010316405
# 6 rere_cm -0.3300292 0.08252020 -0.5048648 -0.026267544 -0.3388155 0.08482095 -0.5054298 -0.026664633

# OLD: 
#       cond    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL
# 1    slre -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.044027732
# 2    fere -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.001785989
# 3    rere -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.034812727
# 4 slre_cm -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.004395259
# 5 fere_cm -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.010775095
# 6 rere_cm -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.014537828

###### Join CI & point estimates for RE med/outcome ------------------------
# This step merges the calculated RE model CIs with existing data.

# Merge RE results with existing results dataframe
results_DF_RE <- merge(results_DF[results_DF$cond %in% c("slre", "fere", "rere", "slre_cm", "fere_cm", "rere_cm"), ], 
                       results_DF_RE)

# Clean up environment (Remove all objects except 'data', 'results_DF', 'results_DF_RE', and functions)
rm(list = setdiff(ls(), c("data", "results_DF", "results_DF_RE", lsf.str())))

#### SL & FE Mediator/Outcome Models CI ------------------------------------
# This section focuses on the calculation and joining of CIs for SL & FE models.

###### Join CI & point estimates for SL & FE med/outcome -------------------
# This block imports data for SL and FE models, calculates CIs, and combines results.

# Define the list of file names
list_names <- c("slsl_ci_TNIE", "slsl_ci_PNIE", 
                "fesl_ci_TNIE", "fesl_ci_PNIE",
                "resl_ci_TNIE", "resl_ci_PNIE", 
                "slfe_ci_TNIE", "slfe_ci_PNIE", 
                "fefe_ci_TNIE", "fefe_ci_PNIE",
                "refe_ci_TNIE", "refe_ci_PNIE")

# Import RDS files and store them in a list of lists
lists <- lapply(list_names, function(name) {
  readRDS(paste0("Application/Output/Bootstrap_Temp/", name, ".rds"))
})

# Name the elements of the list for easier access
names(lists) <- list_names

# Extract 'indirect_ci' and store in dataframe
indirect_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$indirect_ci
  data.frame(
    list_name = list_names[i],
    indirect_ci_LL = ci_values[1],  # Lower bound
    indirect_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Extract 'direct_ci' and store in dataframe
direct_df <- do.call(rbind, lapply(seq_along(lists), function(i) {
  ci_values <- lists[[i]]$direct_ci
  data.frame(
    list_name = list_names[i],
    direct_ci_LL = ci_values[1],  # Lower bound
    direct_ci_UL = ci_values[2],  # Upper bound
    stringsAsFactors = FALSE
  )
}))

# Combine both 'indirect' and 'direct' CIs into a single dataframe
combined_df <- merge(indirect_df, direct_df, by = "list_name")

# Modify dataframe: extract effect and cond labels, and reshape to wide format
combined_df <- combined_df %>%
  mutate(effect = substr(list_name, nchar(list_name) - 3, nchar(list_name))) %>%  # Extract effect label
  mutate(list_name = substr(list_name, 1, 4)) %>%  # Extract cond label
  pivot_longer(cols = c(indirect_ci_LL, indirect_ci_UL, direct_ci_LL, direct_ci_UL), 
               names_to = "variable", values_to = "value") %>%
  unite("variable", effect, variable, sep = "_") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  as.data.frame()

# Update column names for clarity
colnames(combined_df) <- c("cond", 
                           "PNIE_LL", "PNIE_UL", "TNDE_LL", "TNDE_UL", 
                           "TNIE_LL", "TNIE_UL", "PNDE_LL", "PNDE_UL")

# Merge the SL & FE results with the existing results dataframe
results_DF_noRE <- merge(results_DF[1:6, ], combined_df, by = "cond")

# Clean up environment (Remove all objects except 'data', 'results_DF', 'results_DF_RE', 'results_DF_noRE', and functions)
rm(list = setdiff(ls(), c("data", "results_DF", "results_DF_RE", "results_DF_noRE", lsf.str())))

#### Create dataframe of final estimates -----------------------------------
# This final step combines the RE and non-RE results and exports/saves the final estimates.

# Combine RE and non-RE results
results_DF <- rbind(results_DF_noRE, results_DF_RE)

# Display the combined results
# print(results_DF)

# Add labels for PS model & Mediator/Outcome model
results_DF <- results_DF %>%
  mutate(
    PS = case_when(
      startsWith(cond, "fe") ~ "Fixed-Effect",
      startsWith(cond, "sl") ~ "Single-Level",
      startsWith(cond, "re") ~ "Random-Effect",
      TRUE ~ NA_character_  # Default case
    ),
    Model = case_when(
      endsWith(cond, "fe") ~ "Fixed-Effect",
      endsWith(cond, "sl") ~ "Single-Level",
      endsWith(cond, "re") ~ "Random-Effect",
      endsWith(cond, "re_cm") ~ "Random-Effect with Cluster Means",
      TRUE ~ NA_character_  # Default case
    )
  )

# Display the final dataframe with PS model & Mediator/Outcome labels
print(results_DF)
#       cond       TNDE       PNDE        PNIE        TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL       PNDE_UL            PS                            Model
# 1     fefe -0.2510769 -0.2514684 -0.11950564 -0.02166266 -0.3255577 0.07700916 -0.5060783 -0.005091589 -0.3352327 0.08488664 -0.5072594 -0.0063891720  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2588860 -0.2591964 -0.12127187 -0.05814495 -0.3246048 0.07121469 -0.5103306 -0.003194003 -0.3211064 0.07322311 -0.5112406  0.0009317933  Fixed-Effect                     Single-Level
# 3     refe -0.2562645 -0.2566047 -0.11588065 -0.02853428 -0.3126527 0.07480969 -0.5032962 -0.012264249 -0.3184659 0.07940688 -0.5031124 -0.0093673542 Random-Effect                     Fixed-Effect
# 4     resl -0.2987412 -0.2987693 -0.09161688 -0.03054913 -0.2646908 0.08650299 -0.5280205 -0.062530159 -0.2778580 0.09129741 -0.5305844 -0.0650180918 Random-Effect                     Single-Level
# 5     slfe -0.2325015 -0.2327627 -0.11816059 -0.08187151 -0.3050444 0.07659792 -0.4810479  0.006200502 -0.3079430 0.07819370 -0.4822922  0.0065073266  Single-Level                     Fixed-Effect
# 6     slsl -0.3259836 -0.3259562 -0.04471103 -0.01584023 -0.2067034 0.13493934 -0.5633203 -0.085134771 -0.2155183 0.13529510 -0.5648185 -0.0870138644  Single-Level                     Single-Level
# 7     fere -0.2562957 -0.2566773 -0.11918980 -0.03342662 -0.3426149 0.07165350 -0.5128617 -0.019704421 -0.3590603 0.07522094 -0.5137532 -0.0204665279  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2490276 -0.2494285 -0.12421280 -0.04336970 -0.3494602 0.07261279 -0.5102075 -0.008612200 -0.3630953 0.07850820 -0.5110733 -0.0103164046  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2765615 -0.2767321 -0.10180340 -0.02434469 -0.3106561 0.08389139 -0.5176295 -0.050033639 -0.3216423 0.09003098 -0.5162965 -0.0493492021 Random-Effect                    Random-Effect
# 10 rere_cm -0.2557420 -0.2560127 -0.11442310 -0.04131232 -0.3300292 0.08252020 -0.5048648 -0.026267544 -0.3388155 0.08482095 -0.5054298 -0.0266646326 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2719431 -0.2720481 -0.08421877 -0.05117338 -0.2792701 0.09994341 -0.5000102 -0.055682131 -0.2819150 0.10238444 -0.5027937 -0.0536825803  Single-Level                    Random-Effect
# 12 slre_cm -0.2243275 -0.2245039 -0.11099540 -0.08241447 -0.3235201 0.08291188 -0.4695924  0.003501834 -0.3145743 0.08801792 -0.4681990  0.0035533282  Single-Level Random-Effect with Cluster Means

# OLD: 
#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 1     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050 -0.3007998 0.09570006 -0.5084805  0.004777044 -0.3077056 0.10218566 -0.5085235  0.006911447  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2558156 -0.2563010 -0.09886261 -0.0261259555 -0.2921923 0.08181015 -0.5118420 -0.010035295 -0.2961822 0.08569538 -0.5117962 -0.011418954  Fixed-Effect                     Single-Level
# 3     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374 -0.2797301 0.09890207 -0.5069342 -0.000239732 -0.2829443 0.10403690 -0.5077968  0.002558890 Random-Effect                     Fixed-Effect
# 4     resl -0.2959685 -0.2961819 -0.06592642  0.0012385894 -0.2334787 0.10900586 -0.5444086 -0.064872212 -0.2386974 0.11433888 -0.5457297 -0.064180945 Random-Effect                     Single-Level
# 5     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906 -0.2637974 0.10503455 -0.4788932  0.010653283 -0.2650747 0.10698408 -0.4789265  0.014197515  Single-Level                     Fixed-Effect
# 6     slsl -0.3317336 -0.3318426 -0.01230346  0.0247090071 -0.1747149 0.16068737 -0.5693915 -0.084157900 -0.1768113 0.16116459 -0.5715905 -0.083434503  Single-Level                     Single-Level
# 7     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.001785989  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2425750 -0.2431934 -0.10177153 -0.0098265525 -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.010775095  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.034812727 Random-Effect                    Random-Effect
# 10 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.014537828 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.044027732  Single-Level                    Random-Effect
# 12 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.004395259  Single-Level Random-Effect with Cluster Means

#       cond       TNDE       PNDE        PNIE        TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL       PNDE_UL            PS                            Model
# 1     fefe -0.2510769 -0.2514684 -0.11950564 -0.02166266 -0.3255577 0.07700916 -0.5060783 -0.005091589 -0.3352327 0.08488664 -0.5072594 -0.0063891720  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2588860 -0.2591964 -0.12127187 -0.05814495 -0.3246048 0.07121469 -0.5103306 -0.003194003 -0.3211064 0.07322311 -0.5112406  0.0009317933  Fixed-Effect                     Single-Level
# 3     refe -0.2562645 -0.2566047 -0.11588065 -0.02853428 -0.3126527 0.07480969 -0.5032962 -0.012264249 -0.3184659 0.07940688 -0.5031124 -0.0093673542 Random-Effect                     Fixed-Effect
# 4     resl -0.2987412 -0.2987693 -0.09161688 -0.03054913 -0.2646908 0.08650299 -0.5280205 -0.062530159 -0.2778580 0.09129741 -0.5305844 -0.0650180918 Random-Effect                     Single-Level
# 5     slfe -0.2325015 -0.2327627 -0.11816059 -0.08187151 -0.3050444 0.07659792 -0.4810479  0.006200502 -0.3079430 0.07819370 -0.4822922  0.0065073266  Single-Level                     Fixed-Effect
# 6     slsl -0.3259836 -0.3259562 -0.04471103 -0.01584023 -0.2067034 0.13493934 -0.5633203 -0.085134771 -0.2155183 0.13529510 -0.5648185 -0.0870138644  Single-Level                     Single-Level
# 7     fere -0.2562957 -0.2566773 -0.11918980 -0.03342662 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.5126628 -0.0017859895  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2490276 -0.2494285 -0.12421280 -0.04336970 -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.5117886 -0.0107750952  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2765615 -0.2767321 -0.10180340 -0.02434469 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.5246142 -0.0348127272 Random-Effect                    Random-Effect
# 10 rere_cm -0.2557420 -0.2560127 -0.11442310 -0.04131232 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.5094409 -0.0145378283 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2719431 -0.2720481 -0.08421877 -0.05117338 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.5121780 -0.0440277319  Single-Level                    Random-Effect
# 12 slre_cm -0.2243275 -0.2245039 -0.11099540 -0.08241447 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.4688228 -0.0043952591  Single-Level Random-Effect with Cluster Means

# Save the final dataframe of results
write_rds(results_DF, file = "Application/Output/Estimates/Effect-Estimates_bootstrap-CIs.rds")



# Result visuals ----------------------------------------------------------

# TNDE 
## Save visual 
pdf("Application/Output/Visuals/TNDE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(TNDE_LL > 0 | TNDE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = TNDE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = TNDE_LL, xmax = TNDE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Total Natural Direct Effect (TNDE) with 95% Confidence Intervals",
       x = "Total Natural Direct Effect (TNDE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/Visuals/TNDE-Estimates.png", plot = last_plot())


# PNDE 
## Save visual 
pdf("Application/Output/Visuals/PNDE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(PNDE_LL > 0 | PNDE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = PNDE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = PNDE_LL, xmax = PNDE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Pure Natural Direct Effect (PNDE) with 95% Confidence Intervals",
       x = "Pure Natural Direct Effect (PNDE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/Visuals/PNDE-Estimates.png", plot = last_plot())




# TNIE 
## Save visual 
pdf("Application/Output/Visuals/TNIE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(TNIE_LL > 0 | TNIE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = TNIE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = TNIE_LL, xmax = TNIE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Total Natural Indirect Effect (TNIE) with 95% Confidence Intervals",
       x = "Total Natural Indirect Effect (TNIE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/Visuals/TNIE-Estimates.png", plot = last_plot())


# PNIE 
## Save visual 
pdf("Application/Output/Visuals/PNIE-Estimates.pdf")
## Visual 
results_DF %>% 
  mutate(
    # Model = paste(Model, "Mediator/Outcome Model"),  # Append to Model variable
    Zero_Encompasses = ifelse(PNIE_LL > 0 | PNIE_UL < 0, "Below 0", "Includes 0")
  ) %>% 
  ggplot(aes(y = PS, x = PNIE)) +
  geom_point(aes(color = Zero_Encompasses), size = 3) +
  geom_errorbarh(aes(xmin = PNIE_LL, xmax = PNIE_UL, color = Zero_Encompasses), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Pure Natural Indirect Effect (PNIE) with 95% Confidence Intervals",
       x = "Pure Natural Indirect Effect (PNIE)",
       y = "Propensity Score (PS)") +
  facet_wrap(~ Model, ncol = 1) +  # Facet by updated Model variable
  scale_color_manual(values = c("Below 0" = "red", "Includes 0" = "black")) +  # Set colors
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none")  # Remove legend

## 
dev.off()
# Save visual 
ggsave(filename = "Application/Output/Visuals/PNIE-Estimates.png", plot = last_plot())




# Result Table ------------------------------------------------------------

results_DF[, c("PS", "Model", 
               "TNDE", "TNDE_LL", "TNDE_UL", 
               "PNDE", "PNDE_LL", "PNDE_UL", 
               "TNIE", "TNIE_LL", "TNIE_UL", 
               "PNIE", "PNIE_LL", "PNIE_UL")]


