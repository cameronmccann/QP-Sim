#---------------------------------------------------#
# QP Project - Empirical Application 

#' Monte Carlo Confidence Interval Estimation for Mediation Effects
#'
#' This function computes 95% Monte Carlo confidence intervals for mediation effects.
#' It can be used for two different types of mediation estimates:
#'
#' 1. When `effect_type = "PNIE"` (Pure Natural Indirect Effect) the function assumes that
#'    the outcome model does not include a treatment-mediator interaction. In this case:
#'    - PNIE_est = a_path_est * b_path_est
#'    - TNDE_est = outcome model treatment coefficient (direct effect)
#'
#' 2. When `effect_type = "TNIE"` (Total Natural Indirect Effect) the function assumes that
#'    the outcome model includes a treatment-mediator interaction. In this case:
#'    - TNIE_est = a_path_est * b_path_est + interaction_coef
#'    - PNDE_est = outcome model treatment coefficient (direct effect from the interaction model)
#'
#' In both cases the function draws Monte Carlo samples (either independently or jointly)
#' from the estimated sampling distributions (normal approximations) of the relevant parameters
#' and then calculates the 2.5th and 97.5th percentiles as the confidence interval.
#'
#' @param mediator_fit A fitted mediator model object.
#' @param outcome_fit A fitted outcome model object.
#' @param a_coef_name Character. The name of the treatment coefficient in the mediator model (default "sportPartic_w1").
#' @param b_coef_name Character. The name of the mediator (b–path) coefficient in the outcome model (default "selfEst_w3_sc").
#' @param direct_coef_name Character. The name of the treatment coefficient in the outcome model (default "sportPartic_w1").
#' @param n_MC Integer. Number of Monte Carlo simulation iterations (default 1000).
#' @param seed_MC Integer. Seed for reproducibility (default 123456).
#' @param outcome_model_type Character. How to extract the outcome model covariance:
#'        "SL" or "FE" (use vcov()) versus "RE" or "RE-Mean" (use outcome_fit$cov_mat).
#' @param effect_type Character. Either "PNIE" (for the non-interaction models) or "TNIE"
#'        (for the interaction models). See details.
#' @param interaction_coef_name Character. When effect_type = "TNIE", the name of the interaction
#'        coefficient in the outcome model (default "selfEst_w3_sc:sportPartic_w1").
#' @param use_joint Logical. If TRUE, simulation is done jointly using the block–diagonal covariance
#'        matrix; if FALSE (default), parameters are drawn independently.
#'
#' @return A list containing:
#'   - If effect_type = "PNIE":
#'       \item{a_path_est}{Mediator model treatment estimate.}
#'       \item{a_path_se}{Standard error for a_path.}
#'       \item{b_path_est}{Outcome model mediator estimate.}
#'       \item{b_path_se}{Standard error for b_path.}
#'       \item{direct_est}{Outcome model treatment estimate (TNDE).}
#'       \item{direct_se}{Standard error for direct effect.}
#'       \item{PNIE_est}{Point estimate for PNIE (a_path_est * b_path_est).}
#'       \item{TNDE_est}{Point estimate for TNDE (direct_est).}
#'       \item{indirect_CI}{95% Monte Carlo CI for PNIE.}
#'       \item{direct_CI}{95% Monte Carlo CI for TNDE.}
#'
#'   - If effect_type = "TNIE":
#'       \item{a_path_est}{Mediator model treatment estimate.}
#'       \item{a_path_se}{Standard error for a_path.}
#'       \item{b_path_est}{Outcome model mediator estimate.}
#'       \item{b_path_se}{Standard error for b_path.}
#'       \item{interaction_coef}{Outcome model interaction coefficient.}
#'       \item{interaction_coef_se}{Standard error for the interaction coefficient.}
#'       \item{direct_est}{Outcome model treatment estimate (PNDE).}
#'       \item{direct_se}{Standard error for direct effect.}
#'       \item{TNIE_est}{Point estimate for TNIE (a_path_est * b_path_est + interaction_coef).}
#'       \item{PNDE_est}{Point estimate for PNDE (direct_est).}
#'       \item{indirect_CI}{95% Monte Carlo CI for TNIE.}
#'       \item{direct_CI}{95% Monte Carlo CI for PNDE.}
#'
#' @examples
#' # For a fixed-effect mediator/outcome model without interaction:
#' result_PNIE <- monteCarloCI(mediator_fit = med_fefe,
#'                             outcome_fit = out_fefe,
#'                             a_coef_name = "sportPartic_w1",
#'                             b_coef_name = "selfEst_w3_sc",
#'                             direct_coef_name = "sportPartic_w1",
#'                             n_MC = 1000,
#'                             seed_MC = 123456,
#'                             outcome_model_type = "FE",
#'                             effect_type = "PNIE",
#'                             use_joint = FALSE)
#'
#' # For an outcome model with interaction:
#' result_TNIE <- monteCarloCI(mediator_fit = med_fefe_interac,
#'                             outcome_fit = out_fefe_interac,
#'                             a_coef_name = "sportPartic_w1",
#'                             b_coef_name = "selfEst_w3_sc",
#'                             direct_coef_name = "sportPartic_w1",
#'                             n_MC = 1000,
#'                             seed_MC = 123456,
#'                             outcome_model_type = "FE",
#'                             effect_type = "TNIE",
#'                             interaction_coef_name = "selfEst_w3_sc:sportPartic_w1",
#'                             use_joint = FALSE)
#'
monteCarloCIb <- function(mediator_fit, outcome_fit, 
                         a_coef_name = "sportPartic_w1", 
                         b_coef_name = "selfEst_w3_sc", 
                         direct_coef_name = "sportPartic_w1", 
                         n_MC = 1000, 
                         seed_MC = 123456, 
                         PS_model_type = "SL", 
                         mediator_model_type = "SL", 
                         outcome_model_type = "SL", 
                         effect_type = c("PNIE", "TNIE"),
                         interaction_coef_name = "selfEst_w3_sc:sportPartic_w1",
                         use_joint = FALSE, 
                         output_type = c("list", "dataframe")) {
  
  effect_type <- match.arg(effect_type)
  output_type <- match.arg(output_type)
  
  
  # PNDE = Y(1, M(0)) - Y(0, M(0)) = c1 + c3*a0
  # TNDE = Y(1, M(1)) - Y(0, M(1)) = c1 + c3*(a0 + a1)
  # PNIE = Y(0, M(1)) - Y(0, M(0)) = c2*a1
  # TNIE = Y(1, M(1)) - Y(1, M(0)) = a1*(c2 + c3)
  
  if (effect_type == "TNIE") {
    a0_est <- as.numeric(summary(mediator_fit)$coef["(Intercept)", "Estimate"])
    a1_est <- as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Estimate"]) # a_path_est
    c0_est <- as.numeric(summary(outcome_fit)$coef["(Intercept)", "Estimate"])
    c1_est <- as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Estimate"]) # c_path_est
    c2_est <- as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Estimate"]) # b_path_est
    c3_est <- as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
    
    # Point estimates of NDE & NIE --------------------------------------------
    # NIE_est <- a_path_est * b_path_est
    # NDE_est <- direct_est
    PNDE_est <- c1_est + c3_est * a0_est
    TNDE_est <- c1_est + c3_est * (a0_est + a1_est)
    PNIE_est <- c2_est * a1_est
    TNIE_est <- a1_est * (c2_est + c3_est)
    
    # Monte Carlo confidence intervals ----------------------------------------
    # covariance matrix
    if (mediator_model_type %in% c("SL", "FE")) {
      med_vcov <- vcov(mediator_fit)
    } else {
      med_vcov <- mediator_fit$cov_mat
    }
    if (outcome_model_type %in% c("SL", "FE")) {
      out_vcov <- vcov(outcome_fit)
    } else {
      out_vcov <- outcome_fit$cov_mat
    }
    # subset parameters
    med_params <- c("(Intercept)", "sportPartic_w1")
    out_params <- c("(Intercept)", "sportPartic_w1", "selfEst_w3_sc", "selfEst_w3_sc:sportPartic_w1")
    med_vcov_sub <- med_vcov[med_params, med_params, drop = FALSE]
    out_vcov_sub <- out_vcov[out_params, out_params, drop = FALSE]
    
    # c) Build the mean vector in the correct order
    big_mean <- c(
      a0_est,          # (Intercept) in mediator model
      a1_est,          # t in mediator model
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
    
    # Store results  -----------------------------------------------------------
    
    results <- data.frame(
      # Analysis conditions
      analysisCond = paste0(PS_model_type, "_", mediator_model_type, "_", outcome_model_type), 
      PS = PS_model_type,          # Propensity Score model
      medmodel = mediator_model_type, 
      outModel = outcome_model_type, # Outcome model
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
      
      
      # Mediator model intercept & T (in case needed)
      intercept_medModel_est = as.numeric(summary(mediator_fit)$coef["(Intercept)", "Estimate"]),
      intercept_medModel_se  = as.numeric(summary(mediator_fit)$coef["(Intercept)", "Std. Error"]),
      t_medModel_est         = as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Estimate"]),
      t_medModel_se          = as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Std. Error"]),
      # Outcome model intercept, T, M, T:M
      intercept_outModel_est = as.numeric(summary(outcome_fit)$coef["(Intercept)", "Estimate"]),
      intercept_outModel_se  = as.numeric(summary(outcome_fit)$coef["(Intercept)", "Std. Error"]),
      t_outModel_est         = as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Estimate"]),
      t_outModel_se          = as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Std. Error"]),
      m_outModel_est         = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Estimate"]),
      m_outModel_se          = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Std. Error"]),
      tm_outModel_est        = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"]),
      tm_outModel_se         = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Std. Error"])
    )
    
  } else {
    
    a0_est <- as.numeric(summary(mediator_fit)$coef["(Intercept)", "Estimate"])
    a1_est <- as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Estimate"]) # a_path_est
    c0_est <- as.numeric(summary(outcome_fit)$coef["(Intercept)", "Estimate"])
    c1_est <- as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Estimate"]) # c_path_est
    c2_est <- as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Estimate"]) # b_path_est
    c3_est <- as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
    
    # Point estimates of NDE & NIE --------------------------------------------
    # NIE_est <- a_path_est * b_path_est
    # NDE_est <- direct_est
    # PNDE_est <- c1_est + c3_est * a0_est
    TNDE_est <- c1_est + c3_est * (a0_est + a1_est)
    PNIE_est <- c2_est * a1_est
    # TNIE_est <- a1_est * (c2_est + c3_est)
    
    # Monte Carlo confidence intervals ----------------------------------------
    # covariance matrix
    if (mediator_model_type %in% c("SL", "FE")) {
      med_vcov <- vcov(mediator_fit)
    } else {
      med_vcov <- mediator_fit$cov_mat
    }
    if (outcome_model_type %in% c("SL", "FE")) {
      out_vcov <- vcov(outcome_fit)
    } else {
      out_vcov <- outcome_fit$cov_mat
    }
    # subset parameters
    med_params <- c("(Intercept)", "sportPartic_w1")
    out_params <- c("(Intercept)", "sportPartic_w1", "selfEst_w3_sc", "selfEst_w3_sc:sportPartic_w1")
    med_vcov_sub <- med_vcov[med_params, med_params, drop = FALSE]
    out_vcov_sub <- out_vcov[out_params, out_params, drop = FALSE]
    
    # c) Build the mean vector in the correct order
    big_mean <- c(
      a0_est,          # (Intercept) in mediator model
      a1_est,          # t in mediator model
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
    # PNDE_draw <- c1_draw + c3_draw * a0_draw
    TNDE_draw <- c1_draw + c3_draw * (a0_draw + a1_draw)
    PNIE_draw <- c2_draw * a1_draw
    # TNIE_draw <- a1_draw * (c2_draw + c3_draw)
    
    # h) Compute MC-based 95% CIs
    # PNDE_CI <- stats::quantile(PNDE_draw, c(0.025, 0.975))
    TNDE_CI <- stats::quantile(TNDE_draw, c(0.025, 0.975))
    PNIE_CI <- stats::quantile(PNIE_draw, c(0.025, 0.975))
    # TNIE_CI <- stats::quantile(TNIE_draw, c(0.025, 0.975))
    
    # Store results  -----------------------------------------------------------
    
    results <- data.frame(
      # Analysis conditions
      analysisCond = paste0(PS_model_type, "_", mediator_model_type, "_", outcome_model_type), 
      PS = PS_model_type,          # Propensity Score model
      medmodel = mediator_model_type, 
      outModel = outcome_model_type, # Outcome model
      # PNDE, TNIE, TNDE, & PNIE est & CIs
      # PNDE_est = PNDE_est, 
      # PNDE_LCL = PNDE_CI[1],
      # PNDE_UCL = PNDE_CI[2],
      
      TNDE_est = TNDE_est,
      TNDE_LCL = TNDE_CI[1],
      TNDE_UCL = TNDE_CI[2],
      
      PNIE_est = PNIE_est,
      PNIE_LCL = PNIE_CI[1],
      PNIE_UCL = PNIE_CI[2],
      
      # TNIE_est = TNIE_est,
      # TNIE_LCL = TNIE_CI[1],
      # TNIE_UCL = TNIE_CI[2],
      
      
      # Mediator model intercept & T (in case needed)
      intercept_medModel_est = as.numeric(summary(mediator_fit)$coef["(Intercept)", "Estimate"]),
      intercept_medModel_se  = as.numeric(summary(mediator_fit)$coef["(Intercept)", "Std. Error"]),
      t_medModel_est         = as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Estimate"]),
      t_medModel_se          = as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Std. Error"]),
      # Outcome model intercept, T, M, T:M
      intercept_outModel_est = as.numeric(summary(outcome_fit)$coef["(Intercept)", "Estimate"]),
      intercept_outModel_se  = as.numeric(summary(outcome_fit)$coef["(Intercept)", "Std. Error"]),
      t_outModel_est         = as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Estimate"]),
      t_outModel_se          = as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Std. Error"]),
      m_outModel_est         = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Estimate"]),
      m_outModel_se          = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Std. Error"]),
      tm_outModel_est        = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"]),
      tm_outModel_se         = as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Std. Error"])
    )
    
  }

}
