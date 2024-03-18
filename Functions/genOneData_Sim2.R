#---------------------------------------------------#
# QP Project 
# Data Generation for Simulation 2 
#' 
#' `genOneData_Sim2()` generates clustered data consisting of a level-1 treatment, 
#' 6 level-1 confounders, a level-1 mediator, and a level-1 outcome as well as 
#' a level-2 confounder, with control over the number of clusters, cluster size, 
#' and ICC. Specifically, the dataframe returned from the function consists of the 
#' following variables: observation ID (id); cluster ID (school); 
#' 3 level-1 confounders of T, M, Y relations (x1-x3); 3 level-1 T-Y confounders (x4-x6); 
#' a level-2 T-Y confounder (z); (t_ast); the true propensity score of an observation (ps_true); 
#' level-1 treatment assignment (t); level-1 mediator value (m); and the level-1 outcome measure (y). 
#' 
#' @param num_clust Number of clusters 
#' @param clust_size Cluster size for each cluster (i.e., number of observations per cluster) 
#' @param num_x Number of x (level-1) confounders 
#' @param iccx,icct,iccm,iccy The intraclass correlation (ICC) for x confounders, treatment, mediator, and outcome 
#' @returns Returns a dataframe of generated data 
#' @examples
#' genOneData_Sim2(num_clust = 30, clust_size = 30)
#' 
genOneData_Sim2 <-
  function(num_clust,
           clust_size,
           num_x = 6,
           iccx = 0.2,
           icct = 0.2,
           iccm = 0.2,
           iccy = 0.2) {
    
    # set parameters
    #x_t
    a_x <- 0.2 # on trt 
    #x_m
    b_x <- 0.8; # on med 
    #x_y
    c_x <- 0.8; # on outcome 
    
    #z_t 
    a_z <- 0.5 # on trt 
    #z_y
    c_z <- 1.4 # on outcome 
    
    treat_m <- 1 # trt on med   
    treat_y <- 1.3 # trt on outcome 
    med_y <- 1 # med on outcome     
    sd_logit <- 1.81 #
    
    N <- num_clust * clust_size
    
    
    ## Generate variables
    # X confounders & ID variables 
    xb <- mvtnorm::rmvnorm(n = num_clust, 
                           mean = rep(0, num_x), 
                           sigma = diag(iccx, nrow = num_x))[rep(1:num_clust, each = clust_size), ]
    xe <- mvtnorm::rmvnorm(n = N, 
                           mean = rep(0, num_x), 
                           sigma = diag(1-iccx, nrow = num_x))
    x <- xb + xe
    
    data <- data.frame(
      id = 1:N, 
      school = rep(1:num_clust, each = clust_size), 
      x
    )
    
    colnames(data)[-c(1:2)] <- paste0(rep("x", ncol(data) - 2), # Rename x variables 
                                      1:(ncol(data) - 2))
    
    # Covariate z
    z <- rep(stats::rnorm(num_clust, sd = 1), each = clust_size)
    data$z <- z
    
    # Treatment
    tb <- sd_logit * stats::rnorm(num_clust, sd = sqrt(icct))
    t <-
      rep(tb, each = clust_size) + stats::rlogis(N) * sqrt(1 - icct)
    data$t_ast <-
      sd_logit * a_x * rowSums(data[, grepl(pattern = "^x\\d", colnames(data))]) + # code assumes the x variables have the same relationship strength
      t + sd_logit * a_z * z
    data$ps_true <- stats::pnorm(data$t_ast)  # True PS
    data$t <- 1 * (data$t_ast > 0) # Treatment assignment
    
    # Mediator
    mb <-
      stats::rnorm(num_clust, sd = sqrt(iccm)) # mediator explained by clustering
    mr <-
      rep(mb, each = clust_size) + stats::rnorm(N, sd = sqrt(1 - iccm)) # residual (or variance not explained by clustering)
    data$m <- treat_m * data$t +
      b_x * rowSums(data[, grepl(pattern = "^x[1-3]", colnames(data))]) + # x1-x3 
      mr
    
    # Outcome
    yb <-
      stats::rnorm(num_clust, sd = sqrt(iccy)) # outcome explained by clustering (or random-effect)
    yr <-
      rep(yb, each = clust_size) + stats::rnorm(N, sd = sqrt(1 - iccy)) # residual (or error/epsilon term)
    data$y <- med_y * data$m +
      treat_y * data$t +
      c_x * rowSums(data[, grepl(pattern = "^x\\d", colnames(data))]) + #x1-x6 
      c_z * data$z + yr
    
    return(data)
    
  }



##################################### END ######################################



