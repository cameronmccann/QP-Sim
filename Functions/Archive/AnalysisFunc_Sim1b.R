#---------------------------------------------------#
# QP Project 
# Data Analysis for Simulation 1 
#' 
#' `AnalysisFunc_Sim1b()` analyzes simulated clustered data for Simulation 1 based 
#' on the specified propensity score, mediation, and outcome models. The function 
#' returns, in addition to the specified models, the estimated direct (Total 
#' Natural Direct Effect; TNDE) & indirect (Pure Natural Indirect Effect; PNIE) 
#' as well as the simulation condition (including the level of ICC & cluster size). 
#' The a- & b-path estimates and standard errors are returned too.  
#' 
#' @param PSmodel Model to use for propensity score (SL, FE, RE, or RE-Mean) 
#' @param Medmodel Mediation model to use (SL, FE, RE, or RE-Mean) 
#' @param Outcomemodel Outcome model to use (SL, FE, RE, or RE-Mean) 
#' @param data The data set to analyze 
#' @param condition Dataframe containing a row for each condition & the following 
#' variables (columns): num_clust, clust_size, num_x, & icc 
#' @param condition_num condition number, which corresponds to the row number in 
#' condition dataframe
#' @returns Returns a dataframe containing seed number, replication number, 
#' the PS model used (PS), the mediator and outcome model used (outModel), 
#' the direct effect estimate (NDE_est) and indirect effect estimate (NIE_est). 
#' @examples
#' AnalysisFunc_Sim1(PSmodel = "FE", Medmodel = "FE", Outcomemodel = "FE", data = data)
#' 
AnalysisFunc_Sim1b <- function(
    PSmodel = "FE",
    Medmodel = "FE",
    Outcomemodel = "FE",
    data = data, 
    condition = cond, 
    condition_num = cond_num, 
    n_MC = 1000, 
    seed_MC = 123456
    ) {
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
    data <- merge(data, 
                  aggregate(cbind(x1, x2, x3) ~ school, data = data, FUN = mean), 
                  by = "school",
                  suffixes = c("", "_mean"))
    ## Estimate PS model
    psmod <- lme4::glmer(
      formula = paste0("t ~ ", 
                       paste0(colnames(data)[grep("^x\\d+$", colnames(data))], 
                              collapse = " + "), 
                       " + x1_mean + x2_mean + x3_mean + (1 | school)"), 
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
    data <- cbind(data, L2weight = rep(1, nrow(data)))
    med <- WeMix::mix(
      formula = m ~ t + x1 + x2 + x3 + (1 | school),
      data = data,
      weights = c("iptw", "L2weight")
    )
  }
  ## Random-Effect Mean
  if (Medmodel == "RE-Mean") {
    ### Add mean columns 
    data <- merge(data, 
                  setNames(aggregate(x = data$t, 
                                     by = list(data$school), 
                                     FUN = mean), 
                           c("school", "t_mean")), 
                  by = "school")
    ### Add column of just 1s for level-2 weight
    data <- cbind(data, L2weight = rep(1, nrow(data)))
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
      formula = y ~ m + t + x1 + x2 + x3,
      data = data,
      weights = iptw
    )
  }
  ## Fixed-Effect
  if (Outcomemodel == "FE") {
    out <- glm(
      formula = y ~ m + t + x1 + x2 + x3 + as.factor(school),
      data = data,
      weights = iptw
    )
  }
  ## Random-Effect
  if (Outcomemodel == "RE") {
    ### Add column of just 1s for level-2 weight
    data <- cbind(data, L2weight = rep(1, nrow(data)))
    out <-
      WeMix::mix(
        formula = y ~ m + t + x1 + x2 + x3 + (1 | school),
        data = data,
        weights = c("iptw", "L2weight")
      )
  }
  ## Random-Effect Mean
  if (Outcomemodel == "RE-Mean") {
    ### Add mean columns 
    data <- merge(data, 
                  setNames(aggregate(x = data$m, 
                                     by = list(data$school), 
                                     FUN = mean), 
                           c("school", "m_mean")), 
                  by = "school")
    ### Add column of just 1s for level-2 weight
    data <- cbind(data, L2weight = rep(1, nrow(data)))
    out <-
      WeMix::mix(
        formula = y ~ m + m_mean + t + t_mean + x1 + x2 + x3 + (1 | school),
        data = data,
        weights = c("iptw", "L2weight")
      )
  }
  

  # Extract point estimates (a, b, & direct effect) -------------------------
  # a-path
  a_path_est = summary(med)$coef["t", "Estimate"] 
  a_path_se = summary(med)$coef["t", "Std. Error"] 
  # b-path
  b_path_est = summary(out)$coef["m", "Estimate"] 
  b_path_se = summary(out)$coef["m", "Std. Error"] 
  # direct effect
  direct_est = summary(out)$coef["t", "Estimate"] 
  direct_se = summary(out)$coef["t", "Std. Error"] 
  

  # Point estimates of NDE & NIE --------------------------------------------
  NIE_est <- a_path_est * b_path_est
  NDE_est <- direct_est
  
  # Monte Carlo confidence intervals ----------------------------------------
  set.seed(seed_MC)
  # draw from N(est, SE^2)
  a_draws <- rnorm(n_MC, mean = a_path_est, sd = a_path_se)
  b_draws <- rnorm(n_MC, mean = b_path_est, sd = b_path_se)
  # NIE draws 
  NIE_draws <- a_draws * b_draws
  # NDE draws
  c_draws <- rnorm(n_MC, mean = direct_est, sd = direct_se)
  # CIs
  NIE_CI <- quantile(NIE_draws, c(0.025, 0.975))
  NDE_CI <- quantile(c_draws, c(0.025, 0.975))
  
  # Store results  -----------------------------------------------------------
  results <- c(
    analysisCond = paste0(PSmodel, "_", Medmodel, "_", Outcomemodel),
    PS = PSmodel,
    medmodel = Medmodel, 
    outModel = Outcomemodel,
    # Point estimates 
    NDE_est = NDE_est, 
    NIE_est = NIE_est, 
    # MC intervals 
    NDE_LCL = NDE_CI[1], 
    NDE_UCL = NDE_CI[2], 
    NIE_LCL = NIE_CI[1], 
    NIE_UCL = NIE_CI[2], 
    # condition info 
    ICC = condition[condition_num, "icc"],
    clust_size = condition[condition_num, "clust_size"],
    conditionNum = condition_num, 
    # a- & b-path info 
    a_path_est = a_path_est, 
    a_path_se = a_path_se, 
    b_path_est = b_path_est, 
    b_path_se = b_path_se, 
    # direct path info
    direct_est = direct_est, 
    direct_se = direct_se
    
    # NDE_est = summary(out)$coef["t", "Estimate"],
    # NIE_est = summary(med)$coef["t", "Estimate"] *
    #   summary(out)$coef["m", "Estimate"],
    # ICC = condition[condition_num, "icc"],
    # clust_size = condition[condition_num, "clust_size"],
    # conditionNum = condition_num, 
    # 
    # # Extra info 
    # a_path_est = summary(med)$coef["t", "Estimate"], 
    # a_path_se = summary(med)$coef["t", "Std. Error"], 
    # 
    # b_path_est = summary(out)$coef["m", "Estimate"], 
    # b_path_se = summary(out)$coef["m", "Std. Error"], 
    # 
    # direct_est = summary(out)$coef["t", "Estimate"], 
    # direct_se = summary(out)$coef["t", "Std. Error"] 
  )
  
  return(results)
}


##################################### END ######################################



