#---------------------------------------------------#
# QP Project - Empirical Application 
# Last edited: 2025-06-12
#---------------------------------------------------#

#' Bootstrap Confidence Interval Estimation for Mediation Analysis
#'
#' `bootstrapCIb()` performs percentile bootstrap resampling to estimate confidence 
#' intervals for direct and indirect effects in a mediation analysis using clustered 
#' data with propensity score weighting. The function supports multiple model specifications 
#' including single-level ("SL"), fixed-effects ("FE"), random effects ("RE"), and 
#' random effects with cluster means ("RECM"). 
#' 
#' @param iterations Number of bootstrap iterations to perform (default is 50).
#' @param iptw Inverse probability of treatment weights to be used in model fitting.
#' @param data The dataset containing the variables needed for the mediation analysis.
#' @param model The type of model to fit for the mediation and outcome analyses. Options are:
#'   - `"SL"` for single-level (Super Learner),
#'   - `"FE"` for Fixed Effects,
#'   - `"RE"` for Random Effects,
#'   - `"RECM"` for Random Effects with Cluster Means.
#' @param cores Number of CPU cores to use for parallel processing (default is 2).
#' @param core_seeds Optional vector of seeds for each core; if NULL, random seeds will be generated.
#' @param effect_type Type of indirect effect to estimate: `"PNIE"` for Pure Natural Indirect Effect or 
#'   `"TNIE"` for Total Natural Indirect Effect (default is `"TNIE"`). Deprecated: parameter is forced to TNIE since output provides all effects.
#'
#' @returns A list containing:
#'   - `effect_type`: Type of effect estimated ("PNIE" or "TNIE").
#'   - `analysisCond`: Condition specifying the analysis type (e.g., "iptw_FE_SL").
#'   - `PS`: Propensity score weights applied in the analysis.
#'   - `PNDE_est`, `PNDE_LCL`, `PNDE_UCL`: Estimated Pure Natural Direct Effect and its confidence intervals.
#'   - `TNDE_est`, `TNDE_LCL`, `TNDE_UCL`: Estimated Total Natural Direct Effect and its confidence intervals.
#'   - `PNIE_est`, `PNIE_LCL`, `PNIE_UCL`: Estimated Pure Natural Indirect Effect and its confidence intervals.
#'   - `TNIE_est`, `TNIE_LCL`, `TNIE_UCL`: Estimated Total Natural Indirect Effect and its confidence intervals.
#'   - `mediator_converged_count`: Count of bootstrap iterations where the mediator model converged.
#'   - `outcome_converged_count`: Count of bootstrap iterations where the outcome model converged.
#'   - `both_converged_count`: Count of iterations where both models converged.
#'
#' @example
#' # Example usage for single-level:
#' result_sl <- bootstrapCIb(iterations = 100,
#'                          iptw = iptw_re, 
#'                          data = data, 
#'                          model = "SL", 
#'                          cores = 4, 
#'                          core_seeds = c(4561, 4562, 4563, 4564))
#'
#' # Example usage for fixed-effects:
#' result_fe <- bootstrapCIb(iterations = 100,
#'                          iptw = iptw_re, 
#'                          data = data, 
#'                          model = "FE", 
#'                          cores = 4, 
#'                          core_seeds = c(4561, 4562, 4563, 4564))
#'
#' # Example usage for random effects:
#' result_re <- bootstrapCIb(iterations = 100,
#'                          iptw = iptw_re, 
#'                          data = data, 
#'                          model = "RE", 
#'                          cores = 4, 
#'                          core_seeds = c(4561, 4562, 4563, 4564))
#'
#' # Example usage for random effects with cluster means:
#' result_recm <- bootstrapCIb(iterations = 100,
#'                            iptw = iptw_re, 
#'                            data = data, 
#'                            model = "RECM", 
#'                            cores = 4, 
#'                            core_seeds = c(4561, 4562, 4563, 4564))
#'
bootstrapCIb <- function(iterations = 50, iptw, data, model = "SL", cores = 2, core_seeds = NULL, effect_type = "PNIE") {
  
  # always set to TNIE to have interaction in models
  effect_type <- "TNIE"
  
  # Load required libraries for parallel processing and model fitting
  library(parallel)
  library(WeMix)
  
  # Convert the iptw argument to a string for use in model fitting
  iptw_str <- deparse(substitute(iptw))
  
  # Generate unique seeds for each core if not provided
  if (is.null(core_seeds)) {
    core_seeds <- sample.int(1e6, cores)  # Random seed generation
  } else if (length(core_seeds) != cores) {
    stop("The number of core_seeds must match the number of cores")
  }
  
  # Function to perform a single bootstrap iteration
  bootstrap_iteration <- function(i, core_seed, iptw_str, data, model, effect_type) {
    set.seed(core_seed + i)  # Set a unique seed for reproducibility
    
    # Resample clusters with replacement to create a bootstrap dataset
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    if (model == "SL") {
      # Single-level mediator model
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
            parentalEdu_w1_sc + familyStruct_w1,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) NULL)
      
      # Outcome model for SL
      outcome_formula <- depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
        selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1
      # outcome_formula <- if (effect_type == "TNIE") {
      #   depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
      #     selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
      #     white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1
      #   # depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #   #   parentalEdu_w1_sc + familyStruct_w1
      # } else {
      #   depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #     parentalEdu_w1_sc + familyStruct_w1
      # }
      
      outcome <- tryCatch({
        glm(
          formula = outcome_formula,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) NULL)
      
    } else if (model == "FE") {
      # Fixed-effects mediator model
      mediator <- tryCatch({
        glm(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
            parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2),
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) NULL)
      
      # Outcome model for FE
      outcome_formula <- depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
        selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)
      # outcome_formula <- if (effect_type == "TNIE") {
      #   depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #     parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)
      # } else {
      #   depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #     parentalEdu_w1_sc + familyStruct_w1 + as.factor(CLUSTER2)
      # }
      
      outcome <- tryCatch({
        glm(
          formula = outcome_formula,
          data = data_boot,
          weights = data_boot[[iptw_str]]
        )
      }, error = function(e) NULL)
      
    } else if (model == "RE") {
      # Random effects mediator model
      data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
      mediator <- tryCatch({
        WeMix::mix(
          formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
            parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
          data = data_boot,
          weights = c(iptw_str, "L2weight")
        )
      }, error = function(e) NULL)
      
      # Outcome model for RE
      outcome_formula <- depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + 
        selfEst_w3_sc:sportPartic_w1 + age_w1_sc + sex_w1 + 
        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # outcome_formula <- if (effect_type == "TNIE") {
      #   depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #     parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # } else {
      #   depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
      #     parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # }
      
      outcome <- tryCatch({
        WeMix::mix(
          formula = outcome_formula,
          data = data_boot,
          weights = c(iptw_str, "L2weight")
        )
      }, error = function(e) NULL)
      
    } else if (model == "RECM") {
      # Random effects with cluster means: add level-2 weights and compute cluster means
      data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
      
      # Calculate the cluster means for treatment and mediator variables
      cluster_means_sportPartic_w1 <- tapply(data_boot$sportPartic_w1, data_boot$CLUSTER2, mean, na.rm = TRUE)
      cluster_means_selfEst_w3_sc <- tapply(data_boot$selfEst_w3_sc, data_boot$CLUSTER2, mean, na.rm = TRUE)
      
      # Add these cluster means back to the dataset
      data_boot$cluster_mean_sportPartic_w1 <- cluster_means_sportPartic_w1[as.character(data_boot$CLUSTER2)]
      data_boot$cluster_mean_selfEst_w3_sc <- cluster_means_selfEst_w3_sc[as.character(data_boot$CLUSTER2)]
      
      # Mediator model for RECM
      mediator <- tryCatch({
        WeMix::mix(
          formula = selfEst_w3 ~ sportPartic_w1 + cluster_mean_sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 +
            parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),
          data = data_boot,
          weights = c(iptw_str, "L2weight")
        )
      }, error = function(e) NULL)
      
      # Outcome model for RECM
      outcome_formula <- depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + selfEst_w3_sc:sportPartic_w1 + 
        cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc +
        age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # outcome_formula <- if (effect_type == "TNIE") {
      #   depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc +
      #     age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # } else {
      #   depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + cluster_mean_sportPartic_w1 + cluster_mean_selfEst_w3_sc +
      #     age_w1_sc + sex_w1 + white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
      # }
      
      outcome <- tryCatch({
        WeMix::mix(
          formula = outcome_formula,
          data = data_boot,
          weights = c(iptw_str, "L2weight")
        )
      }, error = function(e) NULL)
      
    } else {
      stop("Invalid model type. Choose one of: 'SL', 'FE', 'RE', 'RECM'.")
    }
    
    # Check convergence of both models
    mediator_converged <- !is.null(mediator)
    outcome_converged <- !is.null(outcome)
    
    # Calculate indirect and direct effects if both models converged
    if (mediator_converged && outcome_converged) {
      # if (effect_type == "TNIE") {
        a0_est <- as.numeric(summary(mediator)$coef["(Intercept)", "Estimate"])
        a1_est <- as.numeric(summary(mediator)$coef["sportPartic_w1", "Estimate"]) # a_path_est
        c0_est <- as.numeric(summary(outcome)$coef["(Intercept)", "Estimate"])
        c1_est <- as.numeric(summary(outcome)$coef["sportPartic_w1", "Estimate"]) # c_path_est
        c2_est <- as.numeric(summary(outcome)$coef["selfEst_w3_sc", "Estimate"]) # b_path_est
        c3_est <- as.numeric(summary(outcome)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
        # NIE_est <- a_path_est * b_path_est
        # NDE_est <- direct_est
        PNDE_est <- c1_est + c3_est * a0_est
        TNDE_est <- c1_est + c3_est * (a0_est + a1_est)
        PNIE_est <- c2_est * a1_est
        TNIE_est <- a1_est * (c2_est + c3_est)
      #   
      #   # # TNIE (Total Natural Indirect Effect)
      #   # indirect_effect <- summary(mediator)$coef["sportPartic_w1", "Estimate"] *
      #   #   (summary(outcome)$coef["selfEst_w3_sc", "Estimate"] +
      #   #      summary(outcome)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"])
      #   # # PNDE (Pure Natural Direct Effect)
      #   # direct_effect <- summary(outcome)$coef["sportPartic_w1", "Estimate"] +
      #   #   summary(outcome)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"] *
      #   #   summary(mediator)$coef["(Intercept)", "Estimate"]
      # } else {
      #   # PNIE (Pure Natural Indirect Effect)
      #   indirect_effect <- summary(mediator)$coef["sportPartic_w1", "Estimate"] *
      #     summary(outcome)$coef["selfEst_w3_sc", "Estimate"]
      #   # TNDE (Total Natural Direct Effect)
      #   # direct_effect <- summary(outcome)$coef["sportPartic_w1", "Estimate"] + summary(outcome)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"] * (summary(mediator_fit)$coef["(Intercept)", "Estimate"] + summary(mediator_fit)$coef["sportPartic_w1", "Estimate"])
      #   direct_effect <- summary(outcome)$coef["sportPartic_w1", "Estimate"]
      # }
    } else {
      indirect_effect <- NA
      direct_effect <- NA
    }
    
    # Return a list of effects and convergence status
    # list(indirect_effect = indirect_effect, direct_effect = direct_effect,
    #      mediator_converged = mediator_converged, outcome_converged = outcome_converged)
    list(PNDE_est = PNDE_est, 
         # PNDE_LCL = PNDE_CI[1],
         # PNDE_UCL = PNDE_CI[2],
         TNDE_est = TNDE_est,
         # TNDE_LCL = TNDE_CI[1],
         # TNDE_UCL = TNDE_CI[2],
         PNIE_est = PNIE_est,
         # PNIE_LCL = PNIE_CI[1],
         # PNIE_UCL = PNIE_CI[2],
         TNIE_est = TNIE_est,
         # TNIE_LCL = TNIE_CI[1],
         # TNIE_UCL = TNIE_CI[2],
         mediator_converged = mediator_converged, 
         outcome_converged = outcome_converged)
  }
  
  # Execute bootstrap iterations in parallel
  results <- mclapply(1:iterations, function(i) {
    core_index <- (i - 1) %% cores + 1  # Assign iterations to cores
    bootstrap_iteration(i, core_seeds[core_index], iptw_str, data, model, effect_type)
  }, mc.cores = cores)
  
  # Extract results from bootstrap iterations
  mediator_converged <- sapply(results, function(x) x$mediator_converged)
  outcome_converged <- sapply(results, function(x) x$outcome_converged)
  PNDE_est <- sapply(results, function(x) x$PNDE_est)
  TNDE_est <- sapply(results, function(x) x$TNDE_est)
  PNIE_est <- sapply(results, function(x) x$PNIE_est)
  TNIE_est <- sapply(results, function(x) x$TNIE_est)
  
  # Calculate percentile bootstrap confidence intervals for effects
  PNDE_CI <- stats::quantile(PNDE_est, probs = c(0.025, 0.975), na.rm = TRUE)
  TNDE_CI <- stats::quantile(TNDE_est, probs = c(0.025, 0.975), na.rm = TRUE)
  PNIE_CI <- stats::quantile(PNIE_est, probs = c(0.025, 0.975), na.rm = TRUE)
  TNIE_CI <- stats::quantile(TNIE_est, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return a comprehensive list of results
  list(
    effect_type = effect_type, 
    analysisCond = paste0(toupper(sub("^iptw_", "", iptw_str)), "_", model, "_", model), 
    PS = paste0(toupper(sub("^iptw_", "", iptw_str))),  
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
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}
