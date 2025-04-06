
# TEST MONTE CARLO CI WITH NEW ANALYSIS FUNCTION (WHICH INCLUDE RE-MEAN PS MODEL) ON SINGLE DATASET








# Set Up (Load packages, functions, &/or data) ----------------------------

# Load Packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Packages 
  doParallel, 
  foreach,
  parallel
)

# Load Functions 
source("Functions/AnalysisFunc_Sim1.R")
source("Functions/AnalysisFunc_Sim1b.R")
source("Functions/AnalysisFunc_Sim1c.R")
source("Functions/AnalysisFunc_Sim1c_test.R")

source("Functions/genOneData_Sim1.R")



# Simulation conditions  --------------------------------------------------

cond <- expand.grid(num_clust = 100,
                    clust_size = c(20, 40, 100),
                    num_x = 3,
                    icc = c(0.05, 0.2, 0.5))



# Set Parameters ----------------------------------------------------------

# ## Initialize DF to store results 
# OverallPar_time <- NULL
# 
# ## Set number of replications/repetitions 
# reps <- 1000 
# 
# ## Create directory to store results & save path 
# dir.create(path = "Output/S1_Simulation-Output")
# path <- "Output/S1_Simulation-Output"


# generate data -----------------------------------------------------------

cond_num <- 2
# Generate data set
data <- genOneData_Sim1(
  num_clust = cond[cond_num, "num_clust"],
  clust_size = cond[cond_num, "clust_size"],
  num_x = cond[cond_num, "num_x"],
  icct = cond[cond_num, "icc"],
  iccm = cond[cond_num, "icc"],
  iccy = cond[cond_num, "icc"]
)




# test function -----------------------------------------------------------

AnalysisFunc_Sim1c_test(PSmodel = "FE", 
                        Medmodel = "FE", 
                        Outcomemodel = "FE", 
                        data = data, 
                        condition = cond, 
                        condition_num = cond_num, 
                        n_MC = 1000, 
                        seed_MC = 123456)

# ══════════════════════════════
#    testing again on single dataset
# ══════════════════════════════

#---------------------------
# 1) Simulate a small dataset
#---------------------------
set.seed(101)  # reproducible

n_clusters <- 5
cluster_size <- 20

df_list <- lapply(1:n_clusters, function(sch) {
  x1 <- rnorm(cluster_size)
  x2 <- rnorm(cluster_size)
  x3 <- rnorm(cluster_size)
  t  <- rbinom(cluster_size, 1, 0.4)       # treatment
  m  <- 0.3*t + 0.5*x1 + rnorm(cluster_size)
  y  <- 2.0*t + 0.4*m + x1 + x2 + rnorm(cluster_size)
  data.frame(
    school = sch,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    t  = t,
    m  = m,
    y  = y
  )
})
test_data <- do.call(rbind, df_list)

#---------------------------
# 2) Create a dummy condition dataframe
#   that your function expects
#---------------------------
cond <- data.frame(
  num_clust   = n_clusters,
  clust_size  = cluster_size,
  num_x       = 3,
  icc         = 0.05  # just a placeholder
)

# We'll say condition_num=1 for convenience
cond_num <- 1

#---------------------------
# 3) Run your function with
#    various model arguments
#---------------------------
res_FE_FE_FE <- AnalysisFunc_Sim1c_test(
  PSmodel       = "FE",
  Medmodel      = "FE",
  Outcomemodel  = "FE",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num,
  n_MC          = 500,     # e.g. 500 draws for speed
  seed_MC       = 999
)
res_FE_FE_FE
# You should see a data.frame with your
# point estimates and MC intervals for
# the direct & indirect effects.

# Try a single-level approach
res_SL_SL_SL <- AnalysisFunc_Sim1c_test(
  PSmodel       = "SL",
  Medmodel      = "SL",
  Outcomemodel  = "SL",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_SL_SL_SL

# And so on...

# Try a RE approach
res_RE_RE_RE <- AnalysisFunc_Sim1c_test(
  PSmodel       = "RE",
  Medmodel      = "RE",
  Outcomemodel  = "RE",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_RE_RE_RE


# Try a FE wtih RE approach
res_FE_FE_RE <- AnalysisFunc_Sim1c_test(
  PSmodel       = "FE",
  Medmodel      = "FE",
  Outcomemodel  = "RE",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_FE_FE_RE


# analyze data ------------------------------------------------------------



# Analyze data 
head(data)

AnalysisFunc_Sim1b(PSmodel = "SL", 
                   Medmodel = "SL", 
                   Outcomemodel = "SL", 
                   data = data, 
                   condition = cond, 
                   condition_num = cond_num)



## MC CI SL PS, med, outcome models ----------------------------------------
# TRYING MONTE CARLO CI FOR SL PS, MED, OUT SCENARIO 

# Fit models 
psmod <- glm(formula = "t ~ x1 + x2 + x3", family = "binomial", data = data)
data$ps <- predict(psmod, type = "response")
data$ps_logit <- predict(psmod, type = "link")
## IPTW
data <- cbind(data, iptw = with(data,
                                (t / ps) +
                                  (1 - t) / (1 - ps)))
# M
med <- glm(formula = m ~ t + x1 + x2 + x3,
             data = data,
             weights = iptw)
# Y 
out <- glm(
  formula = y ~ m + t + x1 + x2 + x3,
  data = data,
  weights = iptw
)

# Monte carlo CI 
# (1) Extract coefs & vcovs:
coefs_med <- coef(med)
coefs_out <- coef(out)
vcov_med  <- vcov(med)
vcov_out  <- vcov(out)

# (2) Simulate:
library(MASS)
n_sims <- 1000
sim_med <- mvrnorm(n = n_sims, mu = coefs_med, Sigma = vcov_med)
sim_out <- mvrnorm(n = n_sims, mu = coefs_out, Sigma = vcov_out)

# Identify columns:
idx_int_med <- which(names(coefs_med) == "(Intercept)")
idx_t_med   <- which(names(coefs_med) == "t")
idx_t_out   <- which(names(coefs_out) == "t")
idx_m_out   <- which(names(coefs_out) == "m")
idx_tm_out  <- which(names(coefs_out) == "m:t")  # check name carefully

# (3) Build effect draws:
PNDE_sim <- TNDE_sim <- TNIE_sim <- PNIE_sim <- numeric(n_sims)
NDE_sim <- NIE_sim <- numeric(n_sims)
for(i in seq_len(n_sims)) {
  alpha0_i <- sim_med[i, idx_int_med]
  alphaT_i <- sim_med[i, idx_t_med]
  betaT_i  <- sim_out[i, idx_t_out]
  betaM_i  <- sim_out[i, idx_m_out]
  # betaTM_i <- sim_out[i, idx_tm_out]
  
  # PNDE_sim[i] <- betaT_i + betaTM_i * alpha0_i
  # TNDE_sim[i] <- betaT_i + betaTM_i * (alpha0_i + alphaT_i)
  # TNIE_sim[i] <- alphaT_i * (betaM_i + betaTM_i)
  # PNIE_sim[i] <- alphaT_i * betaM_i
  
  NDE_sim[i] <- betaT_i 
  NIE_sim[i] <- alphaT_i * betaM_i
  # NDE_est = summary(out)$coef["t", "Estimate"],
  # NIE_est = summary(med)$coef["t", "Estimate"] *
  #   summary(out)$coef["m", "Estimate"],
  
}

# (4) Compute intervals:
# PNDE_CI <- quantile(PNDE_sim, c(0.025, 0.975))
# TNDE_CI <- quantile(TNDE_sim, c(0.025, 0.975))
# TNIE_CI <- quantile(TNIE_sim, c(0.025, 0.975))
# PNIE_CI <- quantile(PNIE_sim, c(0.025, 0.975))
NDE_CI <- quantile(NDE_sim, c(0.025, 0.975))
NIE_CI <- quantile(NIE_sim, c(0.025, 0.975))

# (5) Insert into results data frame
results <- NULL
# results$PNDE_lower  <- PNDE_CI[1]
# results$PNDE_upper  <- PNDE_CI[2]
# results$TNDE_lower  <- TNDE_CI[1]
# results$TNDE_upper  <- TNDE_CI[2]
# results$TNIE_lower  <- TNIE_CI[1]
# results$TNIE_upper  <- TNIE_CI[2]
# results$PNIE_lower  <- PNIE_CI[1]
# results$PNIE_upper  <- PNIE_CI[2]
results$NDE_lower  <- NDE_CI[1]
results$NDE_upper  <- NDE_CI[2]
results$NIE_lower  <- NIE_CI[1]
results$NIE_upper  <- NIE_CI[2]

# return(results)

## MC CI FE PS, med, outcome models ----------------------------------------
# TRYING MONTE CARLO CI FOR FE PS, MED, OUT SCENARIO 

# Fit models 
## Estimate PS model
psmod <- glm(
  formula = paste0("t ~ ",
                   paste0(colnames(data)[grep("^x\\d$", colnames(data))],
                          collapse = " + "),
                   " + as.factor(school)"),
  family = "binomial",
  data = data
)
data$ps <- predict(psmod, type = "response")
data$ps_logit <- predict(psmod, type = "link")
## IPTW
data <- cbind(data, iptw = with(data,
                                (t / ps) +
                                  (1 - t) / (1 - ps)))
# M
## Fixed-Effect
med <- glm(
  formula = m ~ t + x1 + x2 + x3 + as.factor(school),
  data = data,
  weights = iptw
)
# Y 
## Fixed-Effect
out <- glm(
  formula = y ~ m + t + x1 + x2 + x3 + as.factor(school),
  data = data,
  weights = iptw
)
out <- glm(
  formula = y ~ m + t + x1 + x2 + x3,
  data = data,
  weights = iptw
)

# Monte carlo CI 
# (1) Extract coefs & vcovs:
coefs_med <- coef(med)
lme4::fixef(med)
coefs_out <- coef(out)
vcov_med  <- vcov(med)
vcov_out  <- vcov(out)

# (2) Simulate:
library(MASS)
n_sims <- 2000
sim_med <- mvrnorm(n = n_sims, mu = coefs_med, Sigma = vcov_med)
sim_out <- mvrnorm(n = n_sims, mu = coefs_out, Sigma = vcov_out)

# Identify columns:
idx_int_med <- which(names(coefs_med) == "(Intercept)")
idx_t_med   <- which(names(coefs_med) == "t")
idx_t_out   <- which(names(coefs_out) == "t")
idx_m_out   <- which(names(coefs_out) == "m")
# idx_tm_out  <- which(names(coefs_out) == "m:t")  # check name carefully

# (3) Build effect draws:
# PNDE_sim <- TNDE_sim <- TNIE_sim <- PNIE_sim <- numeric(n_sims)
NDE_sim <- NIE_sim <- numeric(n_sims)
for(i in seq_len(n_sims)) {
  alpha0_i <- sim_med[i, idx_int_med]
  alphaT_i <- sim_med[i, idx_t_med]
  betaT_i  <- sim_out[i, idx_t_out]
  betaM_i  <- sim_out[i, idx_m_out]
  # betaTM_i <- sim_out[i, idx_tm_out]
  
  # PNDE_sim[i] <- betaT_i + betaTM_i * alpha0_i
  # TNDE_sim[i] <- betaT_i + betaTM_i * (alpha0_i + alphaT_i)
  # TNIE_sim[i] <- alphaT_i * (betaM_i + betaTM_i)
  # PNIE_sim[i] <- alphaT_i * betaM_i
  
  NDE_sim[i] <- betaT_i 
  NIE_sim[i] <- alphaT_i * betaM_i
  # NDE_est = summary(out)$coef["t", "Estimate"],
  # NIE_est = summary(med)$coef["t", "Estimate"] *
  #   summary(out)$coef["m", "Estimate"],
  
}

# (4) Compute intervals:
# PNDE_CI <- quantile(PNDE_sim, c(0.025, 0.975))
# TNDE_CI <- quantile(TNDE_sim, c(0.025, 0.975))
# TNIE_CI <- quantile(TNIE_sim, c(0.025, 0.975))
# PNIE_CI <- quantile(PNIE_sim, c(0.025, 0.975))
NDE_CI <- quantile(NDE_sim, c(0.025, 0.975))
NIE_CI <- quantile(NIE_sim, c(0.025, 0.975))

# (5) Insert into results data frame
results <- NULL
# results$PNDE_lower  <- PNDE_CI[1]
# results$PNDE_upper  <- PNDE_CI[2]
# results$TNDE_lower  <- TNDE_CI[1]
# results$TNDE_upper  <- TNDE_CI[2]
# results$TNIE_lower  <- TNIE_CI[1]
# results$TNIE_upper  <- TNIE_CI[2]
# results$PNIE_lower  <- PNIE_CI[1]
# results$PNIE_upper  <- PNIE_CI[2]
results$NDE_lower  <- NDE_CI[1]
results$NDE_upper  <- NDE_CI[2]
results$NIE_lower  <- NIE_CI[1]
results$NIE_upper  <- NIE_CI[2]

# return(results)
# compare 
AnalysisFunc_Sim1b(PSmodel = "FE", Medmodel = "FE", Outcomemodel = "FE", 
                   data = data, condition = cond, condition_num = cond_num)
results





??mediation

## Estimate PS model
psmod <- glm(
  formula = paste0("t ~ ",
                   paste0(colnames(data)[grep("^x\\d$", colnames(data))],
                          collapse = " + ")),
  family = "binomial",
  data = data
)
data$ps <- predict(psmod, type = "response")
data$ps_logit <- predict(psmod, type = "link")
## IPTW
data <- cbind(data, iptw = with(data,
                                (t / ps) +
                                  (1 - t) / (1 - ps)))
# SL MEd
med <- glm(formula = m ~ t + x1 + x2 + x3,
           data = data,
           weights = iptw)
# RE Out 
data$L2weight <- 1 #data <- cbind(data, L2weight = rep(1, nrow(data)))
out <-
  WeMix::mix(
    formula = y ~ m + t + x1 + x2 + x3 + (1 | school),
    data = data,
    weights = c("iptw", "L2weight")
  )


# 
vcov(out)


out$cov_mat

res_SL_SL_SL <- AnalysisFunc_Sim1c(
  PSmodel = "SL", 
  Medmodel = "SL",
  Outcomemodel = "SL", 
  data = data,
  condition = cond, 
  condition_num = cond_num
)

res_FE_FE_FE <- AnalysisFunc_Sim1c(
  PSmodel = "FE", 
  Medmodel = "FE",
  Outcomemodel = "FE", 
  data = data,
  condition = cond, 
  condition_num = cond_num
)

res_RE_RE_RE <- AnalysisFunc_Sim1c(
  PSmodel = "RE", 
  Medmodel = "RE",
  Outcomemodel = "RE", 
  data = data,
  condition = cond, 
  condition_num = cond_num
)

res_REMean <- AnalysisFunc_Sim1c(
  PSmodel = "RE-Mean", 
  Medmodel = "RE-Mean",
  Outcomemodel = "RE-Mean", 
  data = data,
  condition = cond, 
  condition_num = cond_num
)


res_SL_SL_SL
res_FE_FE_FE
res_RE_RE_RE
res_REMean





# Try a single-level approach
res_SL_SL_SL <- AnalysisFunc_Sim1c_test(
  PSmodel       = "SL",
  Medmodel      = "SL",
  Outcomemodel  = "SL",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_SL_SL_SL

# And so on...

# Try a RE approach
res_RE_RE_RE <- AnalysisFunc_Sim1c_test(
  PSmodel       = "RE",
  Medmodel      = "RE",
  Outcomemodel  = "RE",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_RE_RE_RE


# Try a FE wtih RE approach
res_FE_FE_RE <- AnalysisFunc_Sim1c_test(
  PSmodel       = "FE",
  Medmodel      = "FE",
  Outcomemodel  = "RE",
  data          = test_data,
  condition     = cond,
  condition_num = cond_num
)
res_FE_FE_RE


# Simulation 1 ------------------------------------------------------------

for (condition in 1:nrow(cond)) {
  
  # set condition number
  cond_num <- condition
  
  # Make/Register cores
  doParallel::registerDoParallel(parallel::detectCores() - 1)
  
  # Conduct Simulation
  par_time <- system.time(
    # Track computation time
    
    cond_Results_DF <- foreach::foreach(
      i = 1:reps,
      .combine = rbind,
      .export = c("genOneData_Sim1", "AnalysisFunc_Sim1")
    ) %dopar% {
      set.seed(paste0(135, i))
      
      # Generate data set
      data <- genOneData_Sim1(
        num_clust = cond[cond_num, "num_clust"],
        clust_size = cond[cond_num, "clust_size"],
        num_x = cond[cond_num, "num_x"],
        icct = cond[cond_num, "icc"],
        iccm = cond[cond_num, "icc"],
        iccy = cond[cond_num, "icc"]
      )
      
      # Analyze data set
      temp_DF <- NULL
      full_DF <- NULL
      
      for (PSmod in c("SL", "FE", "RE")) { # PS Models
        for (Outmod in c("SL", "FE", "RE", "RE-Mean")) { # Med/Outcome Models
          
          temp_DF <- AnalysisFunc_Sim1(
            PSmodel = PSmod,
            Medmodel = Outmod,
            Outcomemodel = Outmod,
            data = data,
            condition = cond,
            condition_num = cond_num
          )
          
          # Combine results for condition
          full_DF <- as.data.frame(rbind(full_DF, temp_DF))
          
        }
      }
      
      # Add extra info to DF 
      results <- as.data.frame(cbind(
        seed = rep(paste0(135, i), nrow(full_DF)),
        rep = rep(i, nrow(full_DF)),
        full_DF,
        true_ps_10pctle = rep(quantile(data$ps_true, probs = c(0.1)), nrow(full_DF)), 
        true_ps_15pctle = rep(quantile(data$ps_true, probs = c(0.15)), nrow(full_DF)), 
        true_ps_50pctle = rep(quantile(data$ps_true, probs = c(0.5)), nrow(full_DF)), 
        true_ps_85pctle = rep(quantile(data$ps_true, probs = c(0.85)), nrow(full_DF)), 
        true_ps_90pctle = rep(quantile(data$ps_true, probs = c(0.9)), nrow(full_DF))
      ))
    }
  )
  
  # Save conditions results
  saveRDS(
    cond_Results_DF,
    file = paste0(
      path, 
      "/S1_Condition-", 
      condition, 
      "-Estimates.rds"
    )
  )
  
  # Print message
  print(paste0("Condition ", condition, " Done! ", 
               "(Progress: ", condition, "/", nrow(cond), " = ", 
               round((condition/nrow(cond))*100), "% Complete)"))
  
  if(condition == nrow(cond)) {
    print("~~~~~ Simulation Complete ~~~~~")
  }
  
  # Log computation time 
  OverallPar_time <- rbind(OverallPar_time,
                           cbind(t(as.data.frame(par_time)), cond_num))
  
}


# Add mins to computation time log & save DF 
OverallPar_time <- as.data.frame(OverallPar_time)
OverallPar_time <- cbind(OverallPar_time, 
                         mins = OverallPar_time[, "elapsed"] / 60)

saveRDS(OverallPar_time,
        file = paste0(path, "/S1_Computation-Time.rds"))


