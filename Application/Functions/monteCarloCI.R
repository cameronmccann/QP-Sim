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
monteCarloCI <- function(mediator_fit, outcome_fit, 
                         a_coef_name = "sportPartic_w1", 
                         b_coef_name = "selfEst_w3_sc", 
                         direct_coef_name = "sportPartic_w1", 
                         n_MC = 1000, 
                         seed_MC = 123456, 
                         PS_model_type = "SL", 
                         outcome_model_type = "SL", 
                         effect_type = c("PNIE", "TNIE"),
                         interaction_coef_name = "selfEst_w3_sc:sportPartic_w1",
                         use_joint = FALSE, 
                         output_type = c("list", "dataframe")) {
  
  effect_type <- match.arg(effect_type)
  output_type <- match.arg(output_type)
  
  # Extract mediator a-path estimate and SE
  mediator_summary <- summary(mediator_fit)$coef
  if (!(a_coef_name %in% rownames(mediator_summary))) {
    stop(paste("Coefficient", a_coef_name, "not found in mediator model."))
  }
  a_path_est <- mediator_summary[a_coef_name, "Estimate"]
  a_path_se  <- mediator_summary[a_coef_name, "Std. Error"]
  
  # Extract outcome model coefficients for treatment and mediator
  outcome_summary <- summary(outcome_fit)$coef
  if (!(direct_coef_name %in% rownames(outcome_summary))) {
    stop(paste("Coefficient", direct_coef_name, "not found in outcome model (direct effect)."))
  }
  if (!(b_coef_name %in% rownames(outcome_summary))) {
    stop(paste("Coefficient", b_coef_name, "not found in outcome model (mediator effect)."))
  }
  direct_est <- outcome_summary[direct_coef_name, "Estimate"]
  direct_se  <- outcome_summary[direct_coef_name, "Std. Error"]
  b_path_est <- outcome_summary[b_coef_name, "Estimate"]
  b_path_se  <- outcome_summary[b_coef_name, "Std. Error"]
  
  # For the "TNIE" type, extract the interaction coefficient and SE.
  if (effect_type == "TNIE") {
    if (!(interaction_coef_name %in% rownames(outcome_summary))) {
      stop(paste("Interaction coefficient", interaction_coef_name, "not found in outcome model."))
    }
    interaction_coef    <- outcome_summary[interaction_coef_name, "Estimate"]
    interaction_coef_se <- outcome_summary[interaction_coef_name, "Std. Error"]
  }
  
  # Compute point estimates
  if (effect_type == "PNIE") {
    indirect_est <- a_path_est * b_path_est  # PNIE estimate
    # For non-interaction models, the direct effect is TNDE.
    direct_effect_est <- direct_est         # TNDE estimate
  } else {  # effect_type == "TNIE"
    indirect_est <- a_path_est * b_path_est + interaction_coef  # TNIE estimate
    direct_effect_est <- direct_est                             # PNDE estimate
  }
  
  set.seed(seed_MC)
  
  #
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
  med_params <- c("(Intercept)", "sportPartic_w1") # allow sportPartic_w1 to be changed?
  out_params <- c("(Intercept)", "sportPartic_w1", "selfEst_w3_sc", "selfEst_w3_sc:sportPartic_w1")
  med_vcov_sub <- med_vcov[med_params, med_params, drop = FALSE]
  out_vcov_sub <- out_vcov[out_params, out_params, drop = FALSE]

  # Build mean vector in correct order
  big_mean <- c(
    as.numeric(summary(mediator_fit)$coef["(Intercept)", "Estimate"]),
    as.numeric(summary(mediator_fit)$coef["sportPartic_w1", "Estimate"]), # a_path_est
    as.numeric(summary(outcome_fit)$coef["(Intercept)", "Estimate"]),
    as.numeric(summary(outcome_fit)$coef["sportPartic_w1", "Estimate"]),  # c_path_est
    as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc", "Estimate"]),  # b_path_est
    as.numeric(summary(outcome_fit)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
  )

  # Build block-diagonal covariance matrix
  block_cov <- as.matrix(Matrix::bdiag(med_vcov_sub, out_vcov_sub))
  # Draw from dist
  set.seed(seed_MC)
  params_draws <- MASS::mvrnorm(n_MC, mu = c(), Sigma = block_cov)

  
  #
  if (!use_joint) {
    # Draw independent samples from each normal distribution
    a_draws      <- rnorm(n_MC, mean = a_path_est, sd = a_path_se)
    b_draws      <- rnorm(n_MC, mean = b_path_est, sd = b_path_se)
    direct_draws <- rnorm(n_MC, mean = direct_est, sd = direct_se)
    if (effect_type == "TNIE") {
      interaction_draws <- rnorm(n_MC, mean = interaction_coef, sd = interaction_coef_se)
      indirect_draws <- a_draws * b_draws + interaction_draws
    } else {
      indirect_draws <- a_draws * b_draws
    }
    # Confidence intervals for the indirect and direct effects
    indirect_CI <- quantile(indirect_draws, c(0.025, 0.975))
    direct_CI   <- quantile(direct_draws, c(0.025, 0.975))
    
  } else {
    # Joint simulation using covariance information
    if (outcome_model_type %in% c("SL", "FE")) {
      vcov_out <- vcov(outcome_fit)
    } else if (outcome_model_type %in% c("RE", "RE-Mean")) {
      vcov_out <- outcome_fit$cov_mat
    } else {
      stop("Invalid outcome_model_type. Choose from 'SL', 'FE', 'RE', or 'RE-Mean'.")
    }
    
    # Identify indices for the outcome parameters
    outcome_coef_names <- rownames(outcome_summary)
    idx_direct <- which(outcome_coef_names == direct_coef_name)
    idx_b      <- which(outcome_coef_names == b_coef_name)
    
    if (length(idx_direct) == 0 || length(idx_b) == 0) {
      stop("Specified coefficient names not found in outcome model covariance matrix.")
    }
    # Build outcome covariance submatrix
    if (effect_type == "TNIE") {
      idx_interaction <- which(outcome_coef_names == interaction_coef_name)
      if (length(idx_interaction) == 0) {
        stop("Interaction coefficient not found in outcome model covariance matrix.")
      }
      out_cov_sub  <- vcov_out[c(idx_direct, idx_b, idx_interaction), c(idx_direct, idx_b, idx_interaction)]
      mean_out_sub <- c(direct_est, b_path_est, interaction_coef)
    } else {
      out_cov_sub  <- vcov_out[c(idx_direct, idx_b), c(idx_direct, idx_b)]
      mean_out_sub <- c(direct_est, b_path_est)
    }
    
    # Create block diagonal covariance with a-path variance
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The 'Matrix' package is required for joint simulation. Please install it.")
    }
    cov_a <- matrix(a_path_se^2, nrow = 1)
    block_cov <- Matrix::bdiag(cov_a, out_cov_sub)
    block_cov <- as.matrix(block_cov)
    big_mean <- if (effect_type == "TNIE") {
      c(a_path_est, mean_out_sub)  # length 3+1=4? Actually: a, direct, b, interaction (3 parameters from outcome if TNIE)
      # Here: a_path + direct_est + b_path_est + interaction_coef
    } else {
      c(a_path_est, mean_out_sub)  # a and two outcome parameters
    }
    
    # Draw joint samples
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("The 'MASS' package is required for joint simulation. Please install it.")
    }
    param_draws <- MASS::mvrnorm(n = n_MC, mu = big_mean, Sigma = block_cov)
    # For effect_type == "PNIE": columns: 1 = a, 2 = direct, 3 = b.
    # For effect_type == "TNIE": columns: 1 = a, 2 = direct, 3 = b, 4 = interaction.
    a_draw_joint      <- param_draws[, 1]
    direct_draw_joint <- param_draws[, 2]
    b_draw_joint      <- param_draws[, if (effect_type == "PNIE") 3 else 3]
    if (effect_type == "TNIE") {
      interaction_draw_joint <- param_draws[, 4]
      indirect_draws <- a_draw_joint * b_draw_joint + interaction_draw_joint
    } else {
      indirect_draws <- a_draw_joint * b_draw_joint
    }
    indirect_CI <- quantile(indirect_draws, c(0.025, 0.975))
    direct_CI   <- quantile(direct_draw_joint, c(0.025, 0.975))
  }
  
  if (effect_type == "PNIE") {
    if (output_type == "list") {
      return(list(
        PS_model_type = PS_model_type, 
        outcome_model_type = outcome_model_type, 
        effect_type = effect_type,
        a_path_est = a_path_est,
        a_path_se  = a_path_se,
        b_path_est = b_path_est,
        b_path_se  = b_path_se,
        direct_est = direct_est,
        direct_se  = direct_se,
        NIE_est   = indirect_est,
        NDE_est   = direct_est,
        NIE_CI = indirect_CI,
        NDE_CI   = direct_CI, 
        mediator_model = mediator_fit, 
        outcome_model = outcome_fit
      ))
    } else {
      return(data.frame(
        PS_model_type = PS_model_type, 
        outcome_model_type = outcome_model_type, 
        effect_type = effect_type,
        NIE_est   = indirect_est,
        NIE_LL = indirect_CI[1], 
        NIE_UL = indirect_CI[2], 
        NDE_est = direct_est,
        NDE_LL = direct_CI[1], 
        NDE_UL = direct_CI[2],
        a_path_est = a_path_est,
        a_path_se  = a_path_se,
        b_path_est = b_path_est,
        b_path_se  = b_path_se,
        direct_est = direct_est,
        direct_se  = direct_se
      ))
    }
    
  } else {  # effect_type == "TNIE"
    if (output_type == "list") {
      return(list(
        PS_model_type = PS_model_type, 
        outcome_model_type = outcome_model_type, 
        effect_type = effect_type,
        a_path_est = a_path_est,
        a_path_se  = a_path_se,
        b_path_est = b_path_est,
        b_path_se  = b_path_se,
        interaction_coef = interaction_coef,
        interaction_coef_se = interaction_coef_se,
        direct_est = direct_est,
        direct_se  = direct_se,
        NIE_est   = indirect_est,
        NDE_est   = direct_est,
        NIE_CI = indirect_CI,
        NDE_CI   = direct_CI, 
        mediator_model = mediator_fit, 
        outcome_model = outcome_fit
      ))
    } else {
      return(data.frame(
        PS_model_type = PS_model_type, 
        outcome_model_type = outcome_model_type, 
        effect_type = effect_type,
        NIE_est   = indirect_est,
        NIE_LL = indirect_CI[1], 
        NIE_UL = indirect_CI[2], 
        NDE_est = direct_est,
        NDE_LL = direct_CI[1], 
        NDE_UL = direct_CI[2],
        a_path_est = a_path_est,
        a_path_se  = a_path_se,
        b_path_est = b_path_est,
        b_path_se  = b_path_se,
        interaction_coef = interaction_coef,
        interaction_coef_se = interaction_coef_se,
        direct_est = direct_est,
        direct_se  = direct_se
      ))
    }
  }
}
