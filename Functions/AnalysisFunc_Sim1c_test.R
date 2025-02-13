# AnalysisFunc_Sim1c_test <- function(
#     PSmodel = "FE",
#     Medmodel = "FE",
#     Outcomemodel = "FE",
#     data, 
#     condition, 
#     condition_num,
#     #-- new arguments --
#     n_MC = 1000,    # number of Monte Carlo draws
#     seed_MC = 12345 # seed for reproducibility
# ) {
#   
#   #--------------------------------------------------------------------------#
#   # 1. Propensity Score Model & IPTW
#   #--------------------------------------------------------------------------#
#   
#   if (PSmodel == "SL") {
#     psmod <- glm(
#       formula = paste0(
#         "t ~ ", 
#         paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + ")
#       ),
#       family = "binomial",
#       data = data
#     )
#     data$ps      <- predict(psmod, type = "response")
#     data$ps_logit <- predict(psmod, type = "link")
#     data$iptw    <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
#   }
#   
#   if (PSmodel == "FE") {
#     psmod <- glm(
#       formula = paste0(
#         "t ~ ",
#         paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + "),
#         " + as.factor(school)"
#       ),
#       family = "binomial",
#       data = data
#     )
#     data$ps      <- predict(psmod, type = "response")
#     data$ps_logit <- predict(psmod, type = "link")
#     data$iptw    <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
#   }
#   
#   if (PSmodel == "RE") {
#     psmod <- lme4::glmer(
#       formula = paste0(
#         "t ~ ",
#         paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + "),
#         " + (1 | school)"
#       ),
#       family = "binomial",
#       data = data
#     )
#     data$ps      <- predict(psmod, type = "response")
#     data$ps_logit <- predict(psmod, type = "link")
#     data$iptw    <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
#   }
#   
#   if (PSmodel == "RE-Mean") {
#     # Add cluster means
#     data <- merge(
#       data, 
#       aggregate(cbind(x1, x2, x3) ~ school, data = data, FUN = mean), 
#       by = "school", 
#       suffixes = c("", "_mean")
#     )
#     psmod <- lme4::glmer(
#       formula = paste0(
#         "t ~ ",
#         paste0(colnames(data)[grep("^\\d+$", colnames(data))], collapse = " + "),
#         " + x1_mean + x2_mean + x3_mean + (1 | school)"
#       ),
#       family = "binomial",
#       data = data
#     )
#     data$ps      <- predict(psmod, type = "response")
#     data$ps_logit <- predict(psmod, type = "link")
#     data$iptw    <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
#   }
#   
#   #--------------------------------------------------------------------------#
#   # 2. Mediation Model
#   #--------------------------------------------------------------------------#
#   
#   if (Medmodel == "SL") {
#     med <- glm(m ~ t + x1 + x2 + x3, data = data, weights = iptw)
#   }
#   
#   if (Medmodel == "FE") {
#     med <- glm(
#       m ~ t + x1 + x2 + x3 + as.factor(school),
#       data = data,
#       weights = iptw
#     )
#   }
#   
#   if (Medmodel == "RE") {
#     data$L2weight <- 1  # for WeMix
#     med <- WeMix::mix(
#       formula = m ~ t + x1 + x2 + x3 + (1 | school),
#       data = data,
#       weights = c("iptw", "L2weight")
#     )
#   }
#   
#   if (Medmodel == "RE-Mean") {
#     # Add cluster means of t
#     data <- merge(
#       data, 
#       setNames(aggregate(x = data$t, by = list(data$school), FUN = mean), 
#                c("school","t_mean")),
#       by = "school"
#     )
#     data$L2weight <- 1
#     med <- WeMix::mix(
#       formula = m ~ t + t_mean + x1 + x2 + x3 + (1 | school),
#       data = data,
#       weights = c("iptw", "L2weight")
#     )
#   }
#   
#   #--------------------------------------------------------------------------#
#   # 3. Outcome Model
#   #--------------------------------------------------------------------------#
#   
#   if (Outcomemodel == "SL") {
#     out <- glm(y ~ m + t + x1 + x2 + x3, data = data, weights = iptw)
#   }
#   
#   if (Outcomemodel == "FE") {
#     out <- glm(
#       y ~ m + t + x1 + x2 + x3 + as.factor(school),
#       data = data,
#       weights = iptw
#     )
#   }
#   
#   if (Outcomemodel == "RE") {
#     data$L2weight <- 1
#     out <- WeMix::mix(
#       formula = y ~ m + t + x1 + x2 + x3 + (1 | school),
#       data = data,
#       weights = c("iptw", "L2weight")
#     )
#   }
#   
#   if (Outcomemodel == "RE-Mean") {
#     # Add cluster means of m
#     data <- merge(
#       data,
#       setNames(aggregate(x = data$m, by = list(data$school), FUN = mean), 
#                c("school","m_mean")),
#       by = "school"
#     )
#     data$L2weight <- 1
#     out <- WeMix::mix(
#       formula = y ~ m + m_mean + t + t_mean + x1 + x2 + x3 + (1 | school),
#       data = data,
#       weights = c("iptw","L2weight")
#     )
#   }
#   
#   #--------------------------------------------------------------------------#
#   # 4. Extract point estimates (a, b, direct effect)
#   #--------------------------------------------------------------------------#
#   
#   # 'med' model: a-path (coefficient for t)
#   a_path_est <- summary(med)$coef["t", "Estimate"]
#   a_path_se  <- summary(med)$coef["t", "Std. Error"]
#   
#   # 'out' model: b-path (coefficient for m) & direct effect (coefficient for t)
#   b_path_est <- summary(out)$coef["m", "Estimate"]
#   b_path_se  <- summary(out)$coef["m", "Std. Error"]
#   
#   direct_est <- summary(out)$coef["t", "Estimate"]
#   direct_se  <- summary(out)$coef["t", "Std. Error"]
#   
#   #--------------------------------------------------------------------------#
#   # 5. Point estimates of NDE & NIE
#   #--------------------------------------------------------------------------#
#   
#   NIE_est <- a_path_est * b_path_est  # "pure" natural indirect effect
#   NDE_est <- direct_est               # from outcome model coefficient of t
#   
#   #--------------------------------------------------------------------------#
#   # 6. Monte Carlo confidence intervals
#   #--------------------------------------------------------------------------#
#   
#   set.seed(seed_MC)  # for reproducibility
#   
#   # draw from Normal(a, se(a)^2) and Normal(b, se(b)^2)
#   a_draws <- rnorm(n_MC, mean = a_path_est, sd = a_path_se)
#   b_draws <- rnorm(n_MC, mean = b_path_est, sd = b_path_se)
#   
#   # Indirect effect = a*b for each draw
#   NIE_draws <- a_draws * b_draws
#   
#   # Direct effect draws
#   c_draws <- rnorm(n_MC, mean = direct_est, sd = direct_se)
#   
#   # CI's from empirical quantiles
#   NIE_CI <- quantile(NIE_draws, c(0.025, 0.975))
#   NDE_CI <- quantile(c_draws,   c(0.025, 0.975))
#   
#   #--------------------------------------------------------------------------#
#   # 7. Prepare output
#   #--------------------------------------------------------------------------#
#   
#   results <- data.frame(
#     analysisCond = paste0(PSmodel, "_", Medmodel, "_", Outcomemodel),
#     PS          = PSmodel,
#     Medmodel    = Medmodel,
#     Outmodel    = Outcomemodel,
#     
#     # point estimates
#     NDE_est     = NDE_est,
#     NIE_est     = NIE_est,
#     
#     # MC intervals
#     NDE_LCL     = NDE_CI[1],
#     NDE_UCL     = NDE_CI[2],
#     NIE_LCL     = NIE_CI[1],
#     NIE_UCL     = NIE_CI[2],
#     
#     # condition info
#     ICC         = condition[condition_num, "icc"],
#     clust_size  = condition[condition_num, "clust_size"],
#     conditionNum= condition_num,
#     
#     # a- and b-path info
#     a_path_est = a_path_est,
#     a_path_se  = a_path_se,
#     b_path_est = b_path_est,
#     b_path_se  = b_path_se,
#     
#     # direct path info
#     direct_est = direct_est,
#     direct_se  = direct_se
#   )
#   
#   return(results)
# }




# MVN approach 
AnalysisFunc_Sim1c_test <- function(#AnalysisFunc_Sim1b_MVN <- function(
    PSmodel = "FE",
    Medmodel = "FE",
    Outcomemodel = "FE",
    data, 
    condition, 
    condition_num,
    n_MC = 1000,         # number of Monte Carlo draws
    seed_MC = 12345      # seed for reproducibility
) {
  #----------------------------------------------------#
  # 1. Propensity Score Model & IPTW
  #----------------------------------------------------#
  
  if (PSmodel == "SL") {
    psmod <- glm(
      formula = paste0(
        "t ~ ", 
        paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + ")
      ),
      family = "binomial",
      data = data
    )
    data$ps       <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    data$iptw     <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
  }
  
  if (PSmodel == "FE") {
    psmod <- glm(
      formula = paste0(
        "t ~ ",
        paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + "),
        " + as.factor(school)"
      ),
      family = "binomial",
      data = data
    )
    data$ps       <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    data$iptw     <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
  }
  
  if (PSmodel == "RE") {
    psmod <- lme4::glmer(
      formula = paste0(
        "t ~ ",
        paste0(colnames(data)[grep("^x\\d$", colnames(data))], collapse = " + "),
        " + (1 | school)"
      ),
      family = "binomial",
      data = data
    )
    data$ps       <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    data$iptw     <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
  }
  
  if (PSmodel == "RE-Mean") {
    # Add cluster means
    data <- merge(
      data, 
      aggregate(cbind(x1, x2, x3) ~ school, data = data, FUN = mean), 
      by = "school", 
      suffixes = c("", "_mean")
    )
    psmod <- lme4::glmer(
      formula = paste0(
        "t ~ ",
        paste0(colnames(data)[grep("^\\d+$", colnames(data))], collapse = " + "),
        " + x1_mean + x2_mean + x3_mean + (1 | school)"
      ),
      family = "binomial",
      data = data
    )
    data$ps       <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    data$iptw     <- with(data, (t / ps) + ((1 - t) / (1 - ps)))
  }
  
  #----------------------------------------------------#
  # 2. Mediation Model
  #----------------------------------------------------#
  
  if (Medmodel == "SL") {
    med <- glm(m ~ t + x1 + x2 + x3, data = data, weights = iptw)
  }
  
  if (Medmodel == "FE") {
    med <- glm(
      m ~ t + x1 + x2 + x3 + as.factor(school),
      data = data,
      weights = iptw
    )
  }
  
  if (Medmodel == "RE") {
    data$L2weight <- 1  # for WeMix
    med <- WeMix::mix(
      formula = m ~ t + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  
  if (Medmodel == "RE-Mean") {
    # Add cluster means of t
    data <- merge(
      data, 
      setNames(aggregate(x = data$t, by = list(data$school), FUN = mean), 
               c("school","t_mean")),
      by = "school"
    )
    data$L2weight <- 1
    med <- WeMix::mix(
      formula = m ~ t + t_mean + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  
  #----------------------------------------------------#
  # 3. Outcome Model
  #----------------------------------------------------#
  
  if (Outcomemodel == "SL") {
    out <- glm(y ~ m + t + x1 + x2 + x3, data = data, weights = iptw)
  }
  
  if (Outcomemodel == "FE") {
    out <- glm(
      y ~ m + t + x1 + x2 + x3 + as.factor(school),
      data = data,
      weights = iptw
    )
  }
  
  if (Outcomemodel == "RE") {
    data$L2weight <- 1
    out <- WeMix::mix(
      formula = y ~ m + t + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  
  if (Outcomemodel == "RE-Mean") {
    # Add cluster means of m
    data <- merge(
      data,
      setNames(aggregate(x = data$m, by = list(data$school), FUN = mean), 
               c("school","m_mean")),
      by = "school"
    )
    data$L2weight <- 1
    out <- WeMix::mix(
      formula = y ~ m + m_mean + t + t_mean + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw","L2weight")
    )
  }
  
  #----------------------------------------------------#
  # 4. Extract point estimates for a, b, and direct
  #----------------------------------------------------#
  
  #  -- Mediator model: "a" = coefficient on "t"
  a_path_est <- summary(med)$coef["t", "Estimate"]
  a_path_se  <- summary(med)$coef["t", "Std. Error"]
  
  #  -- Outcome model: "b" = coefficient on "m"; direct effect = coefficient on "t"
  b_path_est   <- summary(out)$coef["m", "Estimate"]
  b_path_se    <- summary(out)$coef["m", "Std. Error"]
  
  direct_est   <- summary(out)$coef["t", "Estimate"]
  direct_se    <- summary(out)$coef["t", "Std. Error"]
  
  #----------------------------------------------------#
  # 5. Compute point estimates for NIE & NDE
  #----------------------------------------------------#
  
  NIE_est <- a_path_est * b_path_est
  NDE_est <- direct_est
  
  #----------------------------------------------------#
  # 6. Monte Carlo draws via MASS::mvrnorm
  #    with block-diagonal covariance
  #----------------------------------------------------#
  #
  # Because 'a' comes from a separate model, we have no
  # cross-correlation with (b, direct). So we form a 1x1
  # block for 'a' and a 2x2 block for (direct, b).
  # Then combine them with a block-diagonal matrix.
  #
  
  set.seed(seed_MC)
  
  #--- (A) Cov for 'a' from mediator model
  #    We'll treat it as a 1-dimensional normal with Var = a_path_se^2
  cov_a  <- matrix(a_path_se^2, nrow = 1, ncol = 1)
  mean_a <- a_path_est
  
  #--- (B) Cov for (direct, b) from outcome model
  #    We'll extract that 2x2 submatrix from vcov(...) for the outcome model
  #    ensuring we pick the correct parameter names if they differ.
  #    Typically, the rownames might be "t" and "m", or they might have
  #    something like "m" "t" in some order. We'll do this carefully:
  
  out_coef_names <- names(coef(out))
  # We'll identify indices for "t" and "m"
  idx_t <- which(out_coef_names == "t")
  idx_m <- which(out_coef_names == "m")
  # The relevant 2x2 submatrix
  vcov_out <- vcov(out)  # for glm or WeMix, assuming it works
  # we store direct first, then b
  out_cov_sub <- vcov_out[c(idx_t, idx_m), c(idx_t, idx_m)]
  
  # Means in the same order: direct, then b
  mean_out_sub <- c(direct_est, b_path_est)
  
  #--- (C) Block-diagonal combination
  #    We'll use 'Matrix::bdiag' to combine them
  library(Matrix)  # for bdiag
  block_cov <- bdiag(cov_a, out_cov_sub)   # 3x3 block diagonal
  
  # "bdiag" returns a sparse Matrix, so we can convert to a regular matrix
  block_cov <- as.matrix(block_cov)
  
  # Combine means in the same dimension: (a, direct, b)
  big_mean <- c(mean_a, mean_out_sub)
  
  #--- (D) Draw from the 3D distribution
  library(MASS)
  param_draws <- MASS::mvrnorm(n_MC, mu = big_mean, Sigma = block_cov)
  # param_draws columns: [ ,1] = a, [ ,2] = direct, [ ,3] = b
  
  #--- (E) For each draw, compute NDE and NIE
  a_draw     <- param_draws[, 1]
  direct_draw<- param_draws[, 2]
  b_draw     <- param_draws[, 3]
  
  NIE_draw <- a_draw * b_draw
  NDE_draw <- direct_draw
  
  #--- (F) 95% MC confidence intervals
  NIE_CI <- stats::quantile(NIE_draw, c(0.025, 0.975))
  NDE_CI <- stats::quantile(NDE_draw, c(0.025, 0.975))
  
  #----------------------------------------------------#
  # 7. Return results
  #----------------------------------------------------#
  
  results <- data.frame(
    analysisCond = paste0(PSmodel, "_", Medmodel, "_", Outcomemodel),
    PS          = PSmodel,
    Medmodel    = Medmodel,
    Outmodel    = Outcomemodel,
    
    # point estimates
    NDE_est     = NDE_est,
    NIE_est     = NIE_est,
    
    # MC intervals
    NDE_LCL     = NDE_CI[1],
    NDE_UCL     = NDE_CI[2],
    NIE_LCL     = NIE_CI[1],
    NIE_UCL     = NIE_CI[2],
    
    # condition info
    ICC         = condition[condition_num, "icc"],
    clust_size  = condition[condition_num, "clust_size"],
    conditionNum= condition_num,
    
    # a-path info
    a_path_est = a_path_est,
    a_path_se  = a_path_se,
    
    # b-path info
    b_path_est = b_path_est,
    b_path_se  = b_path_se,
    
    # direct effect info
    direct_est = direct_est,
    direct_se  = direct_se
  )
  
  return(results)
}











