#---------------------------------------------------#
# QP Project - Empirical Application 
#' 
#' # Bootstrap Confidence Interval Estimation for Random Effects Mediation Analysis
#' 
#' `bootstrap_ci_re_paral()` performs percentile bootstrap resampling to estimate confidence 
#' intervals for direct and indirect effects in a mediation analysis using clustered 
#' data with random effects. The function allows for parallel execution of bootstrap 
#' iterations. It returns the estimated confidence intervals for the direct 
#' effect (Total Natural Direct Effect; TNDE) and indirect effect (Pure Natural 
#' Indirect Effect; PNIE). The function supports random effects modeling for both 
#' mediator and outcome analyses.
#' 
#' @param iterations Number of bootstrap iterations to perform (default is 50).
#' @param iptw Inverse probability of treatment weights to be used in model fitting.
#' @param data The dataset containing the variables needed for the mediation analysis.
#' @param cores Number of CPU cores to use for parallel processing (default is 2).
#' @param core_seeds Optional vector of seeds for each core; if NULL, random seeds will be generated.
#' @returns A list containing:
#'   - `indirect_ci`: Confidence interval for the indirect effect.
#'   - `direct_ci`: Confidence interval for the direct effect.
#'   - `indirect_effects`: Vector of estimated indirect effects from each bootstrap iteration.
#'   - `direct_effects`: Vector of estimated direct effects from each bootstrap iteration.
#'   - `mediator_converged_count`: Count of bootstrap iterations where the mediator model converged.
#'   - `outcome_converged_count`: Count of bootstrap iterations where the outcome model converged.
#'   - `both_converged_count`: Count of iterations where both models converged.
#' 
#' @example
#' # Example usage of bootstrap_ci_re_paral
#' result_re <- bootstrap_ci_re_paral(iterations = 10,
#'                                     iptw = iptw_re, 
#'                                     data = data, 
#'                                     cores = 6, 
#'                                     core_seeds = c(4561, 4562, 4563, 4564, 4565, 4566))
#' 
bootstrap_ci_re_paral <- function(iterations = 50, iptw, data, cores = 2, core_seeds = NULL) {
  # Load required libraries for parallel processing and mixed models
  library(parallel)
  library(WeMix)
  
  # Convert iptw argument to string for model fitting
  iptw_str <- deparse(substitute(iptw))
  
  # Generate unique seeds for each core if not provided
  if (is.null(core_seeds)) {
    core_seeds <- sample.int(1e6, cores)  # Random seed generation
  } else if (length(core_seeds) != cores) {
    stop("The number of core_seeds must match the number of cores")  # Validate core_seeds length
  }
  
  # Function to perform a single bootstrap iteration
  bootstrap_iteration <- function(i, core_seed, iptw_str, data) {
    set.seed(core_seed + i)  # Set a unique seed for reproducibility
    
    # Resample clusters with replacement to create a bootstrap dataset
    cluster_boot <- sample(unique(data$CLUSTER2), replace = TRUE)
    data_boot <- data[data$CLUSTER2 %in% cluster_boot, ]
    
    # Add a column of ones for level-2 weights
    data_boot <- cbind(data_boot, L2weight = rep(1, nrow(data_boot)))
    
    # Fit the random effects models for the bootstrap sample
    mediator_rere <- tryCatch({
      WeMix::mix(
        formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),  
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)  # Handle errors gracefully
    
    outcome_rere <- tryCatch({
      WeMix::mix(
        formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + 
          parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
        data = data_boot, 
        weights = c(iptw_str, "L2weight")
      )
    }, error = function(e) NULL)  # Handle errors gracefully
    
    # Check convergence of both models
    mediator_converged <- !is.null(mediator_rere)
    outcome_converged <- !is.null(outcome_rere)
    
    # Calculate indirect and direct effects if both models converged
    if (mediator_converged && outcome_converged) {
      indirect_effect <- summary(mediator_rere)$coef["sportPartic_w1", "Estimate"] * 
        summary(outcome_rere)$coef["selfEst_w3_sc", "Estimate"]
      direct_effect <- summary(outcome_rere)$coef["sportPartic_w1", "Estimate"]
    } else {
      indirect_effect <- NA  # Assign NA if not converged
      direct_effect <- NA
    }
    
    # Return a list of effects and convergence status
    list(indirect_effect = indirect_effect, direct_effect = direct_effect, 
         mediator_converged = mediator_converged, outcome_converged = outcome_converged)
  }
  
  # Run bootstrap iterations in parallel
  results <- mclapply(1:iterations, function(i) {
    core_index <- (i - 1) %% cores + 1  # Assign iterations to cores
    bootstrap_iteration(i, core_seeds[core_index], iptw_str, data)
  }, mc.cores = cores)
  
  # Extract results from the list of bootstrap iterations
  indirect_effects <- sapply(results, function(x) x$indirect_effect)
  direct_effects <- sapply(results, function(x) x$direct_effect)
  mediator_converged <- sapply(results, function(x) x$mediator_converged)
  outcome_converged <- sapply(results, function(x) x$outcome_converged)
  
  # Calculate percentile bootstrap confidence intervals for effects
  indirect_ci_rere <- quantile(indirect_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  direct_ci_rere <- quantile(direct_effects, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return a comprehensive list of results
  list(
    indirect_ci = indirect_ci_rere,
    direct_ci = direct_ci_rere,
    indirect_effects = indirect_effects,
    direct_effects = direct_effects,
    mediator_converged_count = sum(mediator_converged),
    outcome_converged_count = sum(outcome_converged),
    both_converged_count = sum(mediator_converged & outcome_converged)
  )
}

##################################### END ######################################
