#---------------------------------------------------#
# QP Project 
# Data Generation for Simulation 1 
#' 
#' `genOneData_Sim1()` generates clustered data consisting of a level-1 treatment, 
#' 3 level-1 confounders, a level-1 mediator, and a level-1 outcome as well as 
#' a level-2 confounder, with control over the number of clusters, cluster size, 
#' and ICC. Specifically, the dataframe returned from this function consists of 
#' the following variables: observation ID (id); cluster ID (school); 
#' 3 level-1 confounders of T, M, Y relations (x1-x3); a level-2 confounder of 
#' T, M, Y relations (z); (t_ast); the true propensity score of an observation (ps_true); 
#' a level-1 treatment assignment (t); level-1 mediator value (m); and the level-1 outcome measure (y). 
#' 
#' @param num_clust Number of clusters
#' @param clust_size Cluster size for each cluster (i.e., number of observations per cluster)
#' @param num_x Number of x (level-1) confounders
#' @param iccx,icct,iccm,iccy The intraclass correlation (ICC) for covariate x, treatment, mediator, and outcome
#' @returns Returns a dataframe of generated data
#' @examples
#' genOneData_Sim1(num_clust = 30, clust_size = 30)
#' 
genOneData_Sim1 <-
  function(num_clust,
           clust_size,
           num_x = 3,
           iccx = 0.2,
           icct = 0.2,
           iccm = 0.2,
           iccy = 0.2) {
    
    ## set parameters
    # x 
    a_x <- 0.25 # on trt R^2=0.02798124  #((a_x^2)*num_x * (1 - 0.2)) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3))
    b_x <- 0.3 # on med R^2=0.02993624  #((b_x^2)*num_x * (1 - 0.2)) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
    c_x <- 0.35 # on outcome R^2=0.03030198  #((c_x^2)*num_x * (1 - 0.2)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
    # z 
    a_z <- 1.03 # on trt R^2=0.197902  #(a_z^2) / ((a_x^2)*num_x + a_z^2 + (pi^2/3)/4 + (pi^2/3)) 
    b_z <- 1.21 # on med R^2=0.202915  #(b_z^2) / ((b_x^2)*num_x + b_z^2 + treat_m^2 + (pi^2/3)/4 + (pi^2/3)) 
    c_z <- 1.4 # on outcome R^2=0.2020132  #(c_z^2) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
    # effects 
    treat_m <- 1.17 # trt on med   R^2=0.1517767  #(treat_m^2 * (1 - 0.2)) / (b_z^2 + treat_m^2 + (b_x^2)*num_x + (pi^2/3)/4 + (pi^2/3)) 
    treat_y <- 1.35 # trt on outcome  R^2=0.1502731  #(treat_y^2 * (1 - 0.2)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3)) 
    med_y <- 1.2 # med on outcome     R^2=0.1187343  #(med_y^2 * (1 - 0.2)) / (c_z^2 + treat_y^2 + (c_x^2)*num_x + med_y^2 + (pi^2/3)/4 + (pi^2/3))  
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
      b_x * rowSums(data[, grepl(pattern = "^x\\d", colnames(data))]) +
      b_z * data$z + mr

    # Outcome
    yb <-
      stats::rnorm(num_clust, sd = sqrt(iccy)) # outcome explained by clustering (or random-effect)
    yr <-
      rep(yb, each = clust_size) + stats::rnorm(N, sd = sqrt(1 - iccy)) # residual (or error/epsilon term)
    data$y <- med_y * data$m +
      treat_y * data$t +
      c_x * rowSums(data[, grepl(pattern = "^x\\d", colnames(data))]) +
      c_z * data$z + yr
    
    
    return(data)
    
  }


##################################### END ######################################


