#---------------------------------------------------#
# QP Project 
# Data Analysis for Simulation 2 
#' 
#' `AnalysisFunc_Sim2()` analyzes simulated clustered data for Simulation 2 
#' based on the specified propensity score, mediation, and outcome models. The 
#' function returns, in addition to the specified models, the estimated 
#' direct (Pure Natural Direct Effect; PNDE) & indirect (Total Natural Indirect 
#' Effect; TNIE) as well as the simulation condition (including the level of ICC 
#' & cluster size). The a- & b-path estimates and standard errors are returned too.  
  #' 
#' @param PSmodel Model to use for propensity score (SL, FE, or RE) 
#' @param Medmodel Mediation model to use (SL, FE, RE, or RE-Mean) 
#' @param Outcomemodel Outcome model to use (SL, FE, RE, or RE-Mean) 
#' @param data The data set to analyze 
#' @param condition DF containing a row for each condition & the following 
#' variables (columns): num_clust, clust_size, num_x, & icc 
#' @param condition_num condition number, which corresponds to the row number 
#' in condition dataframe
#' @returns Returns a dataframe containing seed number, replication number, 
#' the PS model used (PS), the mediator and outcome model used (outModel), 
#' the direct effect estimate (NDE_est) and indirect effect estimate (NIE_est). 
#' @examples
#' AnalysisFunc_Sim2(PSmodel = "FE", Medmodel = "FE", Outcomemodel = "FE", data = data)
#' 
AnalysisFunc_Sim2 <- function(PSmodel = "FE",
                               Medmodel = "FE",
                               Outcomemodel = "FE",
                               data = data, 
                               condition = cond, 
                               condition_num = cond_num) {
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
    data <- cbind(data, L2weight = rep(1, nrow(data)))
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
        formula = y ~ m + m_mean + t + t_mean + t:m + t_mean:m_mean + x1 + x2 + x3 + x4 + x5 + x6 + (1 | school),
        data = data,
        weights = c("iptw", "L2weight")
      )
  }
  
  
  # Store results  -----------------------------------------------------------
  
  results <- data.frame(
    # Analysis conditions
    analysisCond = paste0(PSmodel, "_", Outcomemodel), # Combined PS and outcome model name
    PS = PSmodel,          # Propensity Score model
    outModel = Outcomemodel, # Outcome model
    
    # Direct Effects 
    ## Pure Natural Direct Effect
    PNDE_est = as.numeric(summary(out)$coef["t", "Estimate"] + summary(out)$coef["m:t", "Estimate"] * summary(med)$coef["(Intercept)", "Estimate"]), 
    ## Total Natural Direct Effect
    TNDE_est = as.numeric(summary(out)$coef["t", "Estimate"] + summary(out)$coef["m:t", "Estimate"] * (summary(med)$coef["(Intercept)", "Estimate"] + summary(med)$coef["t", "Estimate"])), 
    
    # Indirect Effects 
    ## Total Natural Indirect Effect
    TNIE_est = as.numeric(summary(med)$coef["t", "Estimate"] * (summary(out)$coef["m", "Estimate"] + summary(out)$coef["m:t", "Estimate"])), 
    ## Pure Natural Indirect Effect
    PNIE_est = as.numeric(summary(med)$coef["t", "Estimate"] * summary(out)$coef["m", "Estimate"]), 
    
    # Simulation condition info 
    ICC = condition[condition_num, "icc"],                # Intraclass Correlation
    clust_size = condition[condition_num, "clust_size"],  # Cluster size
    conditionNum = condition_num,                         # Condition ID 
    
    # Mediator model details
    intercept_medModel_est = as.numeric(summary(med)$coef["(Intercept)", "Estimate"]),   # Intercept estimate
    intercept_medModel_se = as.numeric(summary(med)$coef["(Intercept)", "Std. Error"]),  # Intercept SE
    t_medModel_est = as.numeric(summary(med)$coef["t", "Estimate"]),                     # T-M estimate 
    t_medModel_se = as.numeric(summary(med)$coef["t", "Std. Error"]),                    # T-M SE
    
    # Outcome model details
    intercept_outModel_est = as.numeric(summary(out)$coef["(Intercept)", "Estimate"]),   # Intercept estimate
    intercept_outModel_se = as.numeric(summary(out)$coef["(Intercept)", "Std. Error"]),  # Intercept SE
    t_outModel_est = as.numeric(summary(out)$coef["t", "Estimate"]),                     # T-Y estimate 
    t_outModel_se = as.numeric(summary(out)$coef["t", "Std. Error"]),                    # T-Y SE 
    m_outModel_est = as.numeric(summary(out)$coef["m", "Estimate"]),                     # M-Y estimate 
    m_outModel_se = as.numeric(summary(out)$coef["m", "Std. Error"]),                    # M-Y SE 
    tm_outModel_est = as.numeric(summary(out)$coef["m:t", "Estimate"]),                  # TM interaction estimate 
    tm_outModel_se = as.numeric(summary(out)$coef["m:t", "Std. Error"])                  # TM interaction SE 
    
  )
  
  return(results)
  
}


##################################### END ######################################



