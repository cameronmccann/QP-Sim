

#' Percentile-Bootstrap Confidence Intervals for Mediation Effects
#'
#' This function bootstraps the mediator and outcome models, re-computes
#' PNIE/TNIE and PNDE/TNDE on each resample, and then takes the 2.5th/97.5th
#' percentiles as the confidence bounds.  The return value is a data.frame
#' with exactly the same columns as your monteCarloCIb() output.
#'
#' @param data         Original dataset (must contain all vars used in the models).
#' @param iptw         Name of the weight vector in `data` (string), or NULL if none.
#' @param mediator_formula  A one-sided formula for the mediator model (e.g. ~ sportPartic_w1 + ...)
#' @param outcome_formula   Formula for the outcome model *without* treatment×mediator interaction.
#' @param outcome_formula_int  Formula for the outcome model *with* treatment×mediator interaction.
#' @param a_coef_name        Name of the treatment coef in the mediator model (string).
#' @param b_coef_name        Name of the mediator coef in the outcome model (string).
#' @param direct_coef_name   Name of the treatment coef in the outcome model (string).
#' @param interaction_coef_name  Name of the T×M interaction term in the outcome model (string).
#' @param cluster_var    If non-NULL (string), do cluster bootstrap by sampling levels of this var.
#' @param model          One of "SL" or "FE" (single-level or fixed-effects GLM).
#'                       You can insert your own RE/RECM calls here.
#' @param effect_type    "PNIE" or "TNIE"
#' @param R              Number of bootstrap replicates (default 1000).
#' @param cores          Number of cores for `parallel::mclapply()` (default 1).
#' @param seed           Seed for reproducibility.
#'
#' @return A single-row data.frame with columns
#'   analysisCond, PS, medmodel, outModel,
#'   PNIE_est, PNIE_LCL, PNIE_UCL,
#'   TNDE_est, TNDE_LCL, TNDE_UCL,
#'   PNDE_est, PNDE_LCL, PNDE_UCL,
#'   TNIE_est, TNIE_LCL, TNIE_UCL,
#'   intercept_medModel_est, intercept_medModel_se,
#'   t_medModel_est,         t_medModel_se,
#'   intercept_outModel_est, intercept_outModel_se,
#'   t_outModel_est,         t_outModel_se,
#'   m_outModel_est,         m_outModel_se,
#'   tm_outModel_est,        tm_outModel_se
#'
bootstrapCI <- function(data,
                        iptw = NULL,
                        mediator_formula,
                        outcome_formula,
                        outcome_formula_int,
                        a_coef_name = "sportPartic_w1",
                        b_coef_name = "selfEst_w3_sc",
                        direct_coef_name = "sportPartic_w1",
                        interaction_coef_name = "selfEst_w3_sc:sportPartic_w1",
                        cluster_var = NULL,
                        model = c("SL","FE"),
                        effect_type = c("PNIE","TNIE"),
                        R = 1000,
                        cores = 1,
                        seed = 123456) {
  require(parallel)
  effect_type <- match.arg(effect_type)
  model       <- match.arg(model)
  set.seed(seed)
  
  # ----- 1. Statistic function: fit models on one bootstrap draw and return all paths -----
  stat_fun <- function(ids) {
    # ids is either row-indices or cluster-levels depending on cluster_var
    if (is.null(cluster_var)) {
      d0 <- data[ids, ]
    } else {
      clusters <- unique(data[[cluster_var]])
      sel_cl   <- clusters[ids]
      d0 <- data[data[[cluster_var]] %in% sel_cl, ]
    }
    
    # Fit mediator model
    if (model == "SL") {
      med <- glm(update(mediator_formula, . ~ .), data = d0,
                 weights = if(is.null(iptw)) NULL else d0[[iptw]])
    } else {  # FE
      med <- glm(update(mediator_formula, . ~ . + as.factor(get(cluster_var))),
                 data = d0,
                 weights = if(is.null(iptw)) NULL else d0[[iptw]])
    }
    
    # Fit outcome model (with or without interaction)
    form_out <- if (effect_type=="TNIE") outcome_formula_int else outcome_formula
    if (model == "SL") {
      out <- glm(form_out, data = d0,
                 weights = if(is.null(iptw)) NULL else d0[[iptw]])
    } else {
      out <- glm(update(form_out, . ~ . + as.factor(get(cluster_var))),
                 data = d0,
                 weights = if(is.null(iptw)) NULL else d0[[iptw]])
    }
    
    # Extract coefs
    cf_med <- coef(med)
    cf_out <- coef(out)
    
    a0 <- cf_med["(Intercept)"]
    a1 <- cf_med[a_coef_name]
    
    b1 <- cf_out[b_coef_name]
    c1 <- cf_out[direct_coef_name]
    c3 <- if (effect_type=="TNIE") cf_out[interaction_coef_name] else 0
    
    # Effects
    PNIE <- a1 * b1
    TNDE <- c1 + c3 * (a0 + a1)
    PNDE <- c1 + c3 * a0
    TNIE <- a1 * (b1 + c3)
    
    c(PNIE = PNIE,
      TNDE = TNDE,
      PNDE = PNDE,
      TNIE = TNIE,
      a0 = a0,
      a1 = a1,
      b1 = b1,
      c1 = c1,
      c3 = c3,
      out_int = cf_out["(Intercept)"]
    )
  }
  
  # ----- 2. Set up indices to bootstrap: -----
  if (is.null(cluster_var)) {
    idx_list <- replicate(R, sample(nrow(data), replace = TRUE), simplify = FALSE)
  } else {
    clusters <- unique(data[[cluster_var]])
    idx_list <- replicate(R,
                          sample(seq_along(clusters), replace = TRUE),
                          simplify = FALSE)
  }
  
  # ----- 3. Run bootstrap in parallel: -----
  draws <- mclapply(idx_list, stat_fun, mc.cores = cores)
  draws_mat <- do.call(rbind, draws)
  colnames(draws_mat) <- names(draws[[1]])
  
  # ----- 4. Point estimates on original data: -----
  orig <- stat_fun(
    if (is.null(cluster_var)) seq_len(nrow(data))
    else seq_along(unique(data[[cluster_var]]))
  )
  
  # ----- 5. Percentile CIs: -----
  qs <- apply(draws_mat, 2, function(x) quantile(x, c(.025,.975), na.rm=TRUE))
  
  # ----- 6. Assemble the final data.frame: -----
  res <- data.frame(
    analysisCond = model,
    PS           = model,
    medmodel     = model,
    outModel     = model,
    
    PNIE_est = orig["PNIE"],
    PNIE_LCL = qs["2.5%",    "PNIE"],
    PNIE_UCL = qs["97.5%",   "PNIE"],
    
    TNDE_est = orig["TNDE"],
    TNDE_LCL = qs["2.5%",    "TNDE"],
    TNDE_UCL = qs["97.5%",   "TNDE"],
    
    PNDE_est = orig["PNDE"],
    PNDE_LCL = qs["2.5%",    "PNDE"],
    PNDE_UCL = qs["97.5%",   "PNDE"],
    
    TNIE_est = orig["TNIE"],
    TNIE_LCL = qs["2.5%",    "TNIE"],
    TNIE_UCL = qs["97.5%",   "TNIE"],
    
    intercept_medModel_est = orig["a0"],
    intercept_medModel_se  = sd(draws_mat[,"a0"], na.rm=TRUE),
    t_medModel_est         = orig["a1"],
    t_medModel_se          = sd(draws_mat[,"a1"], na.rm=TRUE),
    
    intercept_outModel_est = orig["out_int"],
    intercept_outModel_se  = sd(draws_mat[,"out_int"], na.rm=TRUE),
    t_outModel_est         = orig["c1"],
    t_outModel_se          = sd(draws_mat[,"c1"], na.rm=TRUE),
    
    m_outModel_est         = orig["b1"],
    m_outModel_se          = sd(draws_mat[,"b1"], na.rm=TRUE),
    tm_outModel_est        = orig["c3"],
    tm_outModel_se         = sd(draws_mat[,"c3"], na.rm=TRUE)
  )
  
  return(res)
}
