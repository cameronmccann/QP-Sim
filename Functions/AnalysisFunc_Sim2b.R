#---------------------------------------------------#
# QP Project 
# Data Analysis for Simulation 2 
#' 
#' `AnalysisFunc_Sim2()` analyzes simulated clustered data for Simulation 2 
#' based on the specified propensity score, mediation, and outcome models. 
#' The function returns:
#' \itemize{
#'   \item The fitted PS, mediator, and outcome models (depending on SL, FE, RE, or RE-Mean).
#'   \item The PNDE, TNDE, PNIE, and TNIE point estimates and their corresponding 
#'         Monte Carlo confidence intervals.
#'   \item Simulation condition identifiers (ICC, cluster size, etc.).
#'   \item Various path estimates with their standard errors.
#' }
#' 
#' @param PSmodel Model to use for propensity score: "SL", "FE", "RE", or "RE-Mean".
#' @param Medmodel Mediation model to use: "SL", "FE", "RE", or "RE-Mean".
#' @param Outcomemodel Outcome model to use: "SL", "FE", "RE", or "RE-Mean".
#' @param data The data set to analyze 
#' @param condition A \code{data.frame} containing a row for each condition with columns 
#'                  \code{num_clust}, \code{clust_size}, \code{num_x}, and \code{icc}.
#' @param condition_num The row number in \code{condition} corresponding to the scenario.
#' @param n_MC Number of Monte Carlo draws to use for confidence intervals (default: 1000).
#' @param seed_MC Random seed for the Monte Carlo draws (default: 123456).
#'
#'#' @examples
#' # Example usage:
#' AnalysisFunc_Sim2b(
#'   PSmodel = "FE", Medmodel = "FE", Outcomemodel = "FE",
#'   data = myData, condition = cond, condition_num = 1
#' )
#' 
AnalysisFunc_Sim2b <- function(PSmodel = "FE",
                              Medmodel = "FE",
                              Outcomemodel = "FE",
                              data = data, 
                              condition = cond, 
                              condition_num = cond_num, 
                              n_MC = 1000, 
                              seed_MC = 123456) {
  # PS Models ---------------------------------------------------------------
  ## Single-Level
  if (PSmodel == "SL") {
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
  }
  ## Fixed-Effect
  if (PSmodel == "FE") {
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
  }
  ## Random-Effect
  if (PSmodel == "RE") {
    ## Estimate PS model
    psmod <- lme4::glmer(
      formula = paste0("t ~ ",
                       paste0(colnames(data)[grep("^x\\d$", colnames(data))],
                              collapse = " + "),
                       " + (1 | school)"),
      family = "binomial",
      data = data
    )
    data$ps <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    ## IPTW
    data <- cbind(data, iptw = with(data,
                                    (t / ps) +
                                      (1 - t) / (1 - ps)))
  }
  ## Random-Effect Mean
  if (PSmodel == "RE-Mean") {
    ### Add mean columns 
    data$x1_mean <- ave(data$x1, data$school, FUN = mean)
    data$x2_mean <- ave(data$x2, data$school, FUN = mean)
    data$x3_mean <- ave(data$x3, data$school, FUN = mean)
    data$x4_mean <- ave(data$x4, data$school, FUN = mean)
    data$x5_mean <- ave(data$x5, data$school, FUN = mean)
    data$x6_mean <- ave(data$x6, data$school, FUN = mean)
    ## Estimate PS model
    psmod <- lme4::glmer(
      formula = paste0("t ~ ", 
                       paste0(colnames(data)[grep("^x\\d+$", colnames(data))], 
                              collapse = " + "), 
                       " + x1_mean + x2_mean + x3_mean + x4_mean + x5_mean + x6_mean + (1 | school)"), 
      family = "binomial", 
      data = data
    )
    data$ps <- predict(psmod, type = "response")
    data$ps_logit <- predict(psmod, type = "link")
    ## IPTW
    data <- cbind(data, iptw = with(data,
                                    (t / ps) +
                                      (1 - t) / (1 - ps)))
  }
  
  # Med models  -------------------------------------------------------------
  ## Single-Level
  if (Medmodel == "SL") {
    med <- glm(formula = m ~ t + x1 + x2 + x3,
               data = data,
               weights = iptw)
  }
  ## Fixed-Effect
  if (Medmodel == "FE") {
    med <- glm(
      formula = m ~ t + x1 + x2 + x3 + as.factor(school),
      data = data,
      weights = iptw
    )
  }
  ## Random-Effect
  if (Medmodel == "RE") {
    ### Add column of just 1s for level-2 weight
    data$L2weight <- 1 #data <- cbind(data, L2weight = rep(1, nrow(data)))
    med <- WeMix::mix(
      formula = m ~ t + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  ## Random-Effect Mean
  if (Medmodel == "RE-Mean") {
    ### Add mean columns 
    data$t_mean <- ave(data$t, data$school, FUN = mean)
    # data <- merge(data, 
    #               setNames(aggregate(x = data$t, 
    #                                  by = list(data$school), 
    #                                  FUN = mean), 
    #                        c("school", "t_mean")), 
    #               by = "school")
    ### Add column of just 1s for level-2 weight
    data$L2weight <- 1 #data <- cbind(data, L2weight = rep(1, nrow(data)))
    med <- WeMix::mix(
      formula = m ~ t + t_mean + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  
  
  # Outcome model  ----------------------------------------------------------
  ## Single-Level
  if (Outcomemodel == "SL") {
    out <- glm(
      formula = y ~ m + t + t:m + x1 + x2 + x3 + x4 + x5 + x6,
      data = data,
      weights = iptw
    )
  }
  ## Fixed-Effect
  if (Outcomemodel == "FE") {
    out <- glm(
      formula = y ~ m + t + t:m + x1 + x2 + x3 + x4 + x5 + x6 + as.factor(school),
      data = data,
      weights = iptw
    )
  }
  ## Random-Effect
  if (Outcomemodel == "RE") {
    ### Add column of just 1s for level-2 weight
    data$L2weight <- 1 #data <- cbind(data, L2weight = rep(1, nrow(data)))
    out <-
      WeMix::mix(
        formula = y ~ m + t + t:m + x1 + x2 + x3 + x4 + x5 + x6 + (1 | school),
        data = data,
        weights = c("iptw", "L2weight")
      )
  }
  ## Random-Effect Mean
  if (Outcomemodel == "RE-Mean") {
    ### Add mean columns 
    data$t_mean <- ave(data$t, data$school, FUN = mean)
    data$m_mean <- ave(data$m, data$school, FUN = mean)
    # data <- merge(data, 
    #               setNames(aggregate(x = data$m, 
    #                                  by = list(data$school), 
    #                                  FUN = mean), 
    #                        c("school", "m_mean")), 
    #               by = "school")
    ### Add column of just 1s for level-2 weight
    data <- cbind(data, L2weight = rep(1, nrow(data)))
    out <-
      WeMix::mix(
        formula = y ~ m + m_mean + t + t_mean + t:m + t_mean:m_mean + x1 + x2 + x3 + x4 + x5 + x6 + (1 | school),
        data = data,
        weights = c("iptw", "L2weight")
      )
  }
  
  # Extract point estimates (a, b, & direct effect) -------------------------
  # # a-path
  # a_path_est = summary(med)$coef["t", "Estimate"] 
  # a_path_se = summary(med)$coef["t", "Std. Error"] 
  # # b-path
  # b_path_est = summary(out)$coef["m", "Estimate"] 
  # b_path_se = summary(out)$coef["m", "Std. Error"] 
  # # direct effect
  # direct_est = summary(out)$coef["t", "Estimate"] 
  # direct_se = summary(out)$coef["t", "Std. Error"] 
  
  # PNDE = Y(1, M(0)) - Y(0, M(0)) = c1 + c3*a0
  # TNDE = Y(1, M(1)) - Y(0, M(1)) = c1 + c3*(a0 + a1)
  # PNIE = Y(0, M(1)) - Y(0, M(0)) = c2*a1
  # TNIE = Y(1, M(1)) - Y(1, M(0)) = a1*(c2 + c3)
  
  a0_est <- as.numeric(summary(med)$coef["(Intercept)", "Estimate"])
  a1_est <- as.numeric(summary(med)$coef["t", "Estimate"]) # a_path_est
  c0_est <- as.numeric(summary(out)$coef["(Intercept)", "Estimate"])
  c1_est <- as.numeric(summary(out)$coef["t", "Estimate"]) # c_path_est
  c2_est <- as.numeric(summary(out)$coef["m", "Estimate"]) # b_path_est
  c3_est <- as.numeric(summary(out)$coef["m:t", "Estimate"])
  
  # Point estimates of NDE & NIE --------------------------------------------
  # NIE_est <- a_path_est * b_path_est
  # NDE_est <- direct_est
  PNDE_est <- c1_est + c3_est * a0_est
  TNDE_est <- c1_est + c3_est * (a0_est + a1_est)
  PNIE_est <- c2_est * a1_est
  TNIE_est <- a1_est * (c2_est + c3_est)
  
  # Monte Carlo confidence intervals ----------------------------------------
  # covariance matrix
  if (Medmodel %in% c("SL", "FE")) {
    med_vcov <- vcov(med)
  } else {
    med_vcov <- med$cov_mat
  }
  if (Outcomemodel %in% c("SL", "FE")) {
    out_vcov <- vcov(out)
  } else {
    out_vcov <- out$cov_mat
  }
  # subset parameters
  med_params <- c("(Intercept)", "t")
  out_params <- c("(Intercept)", "t", "m", "m:t")
  med_vcov_sub <- med_vcov[med_params, med_params, drop = FALSE]
  out_vcov_sub <- out_vcov[out_params, out_params, drop = FALSE]
  
  # c) Build the mean vector in the correct order
  big_mean <- c(
    a0_est,          # (Intercept) in mediator model
    a1_est,      # t in mediator model
    c0_est,          # (Intercept) in outcome model
    c1_est,          # t in outcome model
    c2_est,          # m in outcome model
    c3_est           # m:t in outcome model
  )
  
  # d) Build block-diagonal covariance matrix
  block_cov <- as.matrix(Matrix::bdiag(med_vcov_sub, out_vcov_sub))
  
  # e) Draw from the 6D normal distribution
  set.seed(seed_MC)
  param_draws <- MASS::mvrnorm(n_MC, mu = big_mean, Sigma = block_cov)
  
  # f) Label draws for clarity
  # param_draws columns:
  #    1 -> a0 (Intercept in med)
  #    2 -> a1 (T in med)
  #    3 -> c0 (Intercept in out)
  #    4 -> c1 (T in out)
  #    5 -> c2 (M in out)
  #    6 -> c3 (M:T in out)
  a0_draw <- param_draws[, 1]
  a1_draw <- param_draws[, 2]
  c0_draw <- param_draws[, 3]
  c1_draw <- param_draws[, 4]
  c2_draw <- param_draws[, 5]
  c3_draw <- param_draws[, 6]
  
  # g) Compute draws of each effect
  PNDE_draw <- c1_draw + c3_draw * a0_draw
  TNDE_draw <- c1_draw + c3_draw * (a0_draw + a1_draw)
  PNIE_draw <- c2_draw * a1_draw
  TNIE_draw <- a1_draw * (c2_draw + c3_draw)
  
  # h) Compute MC-based 95% CIs
  PNDE_CI <- stats::quantile(PNDE_draw, c(0.025, 0.975))
  TNDE_CI <- stats::quantile(TNDE_draw, c(0.025, 0.975))
  PNIE_CI <- stats::quantile(PNIE_draw, c(0.025, 0.975))
  TNIE_CI <- stats::quantile(TNIE_draw, c(0.025, 0.975))
  
  # set.seed(seed_MC)
  # # draw from N(est, SE^2)
  # a_draws <- rnorm(n_MC, mean = a_path_est, sd = a_path_se)
  # b_draws <- rnorm(n_MC, mean = b_path_est, sd = b_path_se)
  # # NIE draws
  # NIE_draws <- a_draws * b_draws
  # # NDE draws
  # c_draws <- rnorm(n_MC, mean = direct_est, sd = direct_se)
  # # CIs
  # NIE_CI <- quantile(NIE_draws, c(0.025, 0.975))
  # NDE_CI <- quantile(c_draws, c(0.025, 0.975))
  # 
  # set.seed(seed_MC)
  # # cov for a (from mediator model)
  # cov_a <- matrix(a_path_se^2, nrow = 1, ncol = 1)
  # mean_a <- a_path_est
  # # cov for driect & b (from outcome model)
  # out_coef_names <- names(coef(out))
  # idx_t <- which(out_coef_names == "t")
  # idx_m <- which(out_coef_names == "m")
  # # var covariance matrix
  # if (Outcomemodel %in% c("SL", "FE")) {
  #   vcov_out <- vcov(out)  #
  # }
  # if (Outcomemodel %in% c("RE", "RE-Mean")) {
  #   vcov_out <- out$cov_mat
  # }
  # # store direct & b
  # out_cov_sub <- vcov_out[c(idx_t, idx_m), c(idx_t, idx_m)]
  # # means for direct & b (same order as cov)
  # mean_out_sub <- c(direct_est, b_path_est)
  # 
  # #
  # library(Matrix)
  # block_cov <- Matrix::bdiag(cov_a, out_cov_sub) # 3x3 block diagonal
  # # convert sparse matrix to regular matrix
  # block_cov <- as.matrix(block_cov)
  # # combine means
  # big_mean <- c(mean_a, mean_out_sub)
  # 
  # # draw from 3D distr
  # library(MASS) #param_draws columns: [ ,1] = a, [ ,2] = direct, [ ,3] = b
  # param_draws <- MASS::mvrnorm(n_MC, mu = big_mean, Sigma = block_cov)
  # # compute NDE & NIE
  # a_draw     <- param_draws[, 1]
  # direct_draw<- param_draws[, 2]
  # b_draw     <- param_draws[, 3]
  # 
  # NIE_draw <- a_draw * b_draw
  # NDE_draw <- direct_draw
  # 
  # # 95% MC confidence intervals
  # NIE_CI <- stats::quantile(NIE_draw, c(0.025, 0.975))
  # NDE_CI <- stats::quantile(NDE_draw, c(0.025, 0.975))
  
  # Store results  -----------------------------------------------------------
  
  results <- data.frame(
    # Analysis conditions
    analysisCond = paste0(PSmodel, "_", Medmodel, "_", Outcomemodel), 
    PS = PSmodel,          # Propensity Score model
    medmodel = Medmodel, 
    outModel = Outcomemodel, # Outcome model
    # PNDE, TNIE, TNDE, & PNIE est & CIs
    PNDE_est = PNDE_est, 
    PNDE_LCL = PNDE_CI[1],
    PNDE_UCL = PNDE_CI[2],
    
    TNDE_est = TNDE_est,
    TNDE_LCL = TNDE_CI[1],
    TNDE_UCL = TNDE_CI[2],
    
    PNIE_est = PNIE_est,
    PNIE_LCL = PNIE_CI[1],
    PNIE_UCL = PNIE_CI[2],
    
    TNIE_est = TNIE_est,
    TNIE_LCL = TNIE_CI[1],
    TNIE_UCL = TNIE_CI[2],
    # condition info 
    ICC = condition[condition_num, "icc"],
    clust_size = condition[condition_num, "clust_size"],
    conditionNum = condition_num, 
    # Mediator model intercept & T (in case needed)
    intercept_medModel_est = as.numeric(summary(med)$coef["(Intercept)", "Estimate"]),
    intercept_medModel_se  = as.numeric(summary(med)$coef["(Intercept)", "Std. Error"]),
    t_medModel_est         = as.numeric(summary(med)$coef["t", "Estimate"]),
    t_medModel_se          = as.numeric(summary(med)$coef["t", "Std. Error"]),
    # Outcome model intercept, T, M, T:M
    intercept_outModel_est = as.numeric(summary(out)$coef["(Intercept)", "Estimate"]),
    intercept_outModel_se  = as.numeric(summary(out)$coef["(Intercept)", "Std. Error"]),
    t_outModel_est         = as.numeric(summary(out)$coef["t", "Estimate"]),
    t_outModel_se          = as.numeric(summary(out)$coef["t", "Std. Error"]),
    m_outModel_est         = as.numeric(summary(out)$coef["m", "Estimate"]),
    m_outModel_se          = as.numeric(summary(out)$coef["m", "Std. Error"]),
    tm_outModel_est        = as.numeric(summary(out)$coef["m:t", "Estimate"]),
    tm_outModel_se         = as.numeric(summary(out)$coef["m:t", "Std. Error"])
    
    
    # # a- & b-path info 
    # a_path_est = a_path_est, 
    # a_path_se = a_path_se, 
    # b_path_est = b_path_est, 
    # b_path_se = b_path_se, 
    # # direct path info
    # direct_est = direct_est, 
    # direct_se = direct_se
    # 
    # # Direct Effects 
    # ## Pure Natural Direct Effect
    # PNDE_est = as.numeric(summary(out)$coef["t", "Estimate"] + summary(out)$coef["m:t", "Estimate"] * summary(med)$coef["(Intercept)", "Estimate"]), 
    # ## Total Natural Direct Effect
    # TNDE_est = as.numeric(summary(out)$coef["t", "Estimate"] + summary(out)$coef["m:t", "Estimate"] * (summary(med)$coef["(Intercept)", "Estimate"] + summary(med)$coef["t", "Estimate"])), 
    # 
    # # Indirect Effects 
    # ## Total Natural Indirect Effect
    # TNIE_est = as.numeric(summary(med)$coef["t", "Estimate"] * (summary(out)$coef["m", "Estimate"] + summary(out)$coef["m:t", "Estimate"])), 
    # ## Pure Natural Indirect Effect
    # PNIE_est = as.numeric(summary(med)$coef["t", "Estimate"] * summary(out)$coef["m", "Estimate"]), 
    # 
    # # Simulation condition info 
    # ICC = condition[condition_num, "icc"],                # Intraclass Correlation
    # clust_size = condition[condition_num, "clust_size"],  # Cluster size
    # conditionNum = condition_num,                         # Condition ID 
    # 
    # # Mediator model details
    # intercept_medModel_est = as.numeric(summary(med)$coef["(Intercept)", "Estimate"]),   # Intercept estimate
    # intercept_medModel_se = as.numeric(summary(med)$coef["(Intercept)", "Std. Error"]),  # Intercept SE
    # t_medModel_est = as.numeric(summary(med)$coef["t", "Estimate"]),                     # T-M estimate 
    # t_medModel_se = as.numeric(summary(med)$coef["t", "Std. Error"]),                    # T-M SE
    # 
    # # Outcome model details
    # intercept_outModel_est = as.numeric(summary(out)$coef["(Intercept)", "Estimate"]),   # Intercept estimate
    # intercept_outModel_se = as.numeric(summary(out)$coef["(Intercept)", "Std. Error"]),  # Intercept SE
    # t_outModel_est = as.numeric(summary(out)$coef["t", "Estimate"]),                     # T-Y estimate 
    # t_outModel_se = as.numeric(summary(out)$coef["t", "Std. Error"]),                    # T-Y SE 
    # m_outModel_est = as.numeric(summary(out)$coef["m", "Estimate"]),                     # M-Y estimate 
    # m_outModel_se = as.numeric(summary(out)$coef["m", "Std. Error"]),                    # M-Y SE 
    # tm_outModel_est = as.numeric(summary(out)$coef["m:t", "Estimate"]),                  # TM interaction estimate 
    # tm_outModel_se = as.numeric(summary(out)$coef["m:t", "Std. Error"])                  # TM interaction SE 
    
  )
  
  return(results)
  
}


##################################### END ######################################
