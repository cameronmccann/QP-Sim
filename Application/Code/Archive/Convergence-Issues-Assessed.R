


# Test Dr Liu Solution (with healthInsur_gap) -----------------------------


## Data imputation & scaling  ----------------------------------------------

# Need to scale: 
# Scale, rename, & drop non-scaled variable
data <- data %>%
  mutate(
    age_w1_sc = as.vector(scale(age_w1)),
    parentalEdu_w1_sc = as.vector(scale(parentalEdu_w1)),
    feelings_w1_sc = as.vector(scale(feelings_w1))
  ) %>%
  select(-age_w1, -parentalEdu_w1, -feelings_w1)  # Remove original columns if desired

# Select variables in PS model 
t <- data %>% 
  select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_w3, healthInsur_w4, healthInsur_w1)) %>% 
# select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_gap))
# healthInsur_w3, healthInsur_w4, 
# healthInsur_w1))
  select(c("sex_w1", "white_w1", "black_w1", 
           "sportPartic_w1", "selfEst_w3", "depress_w4", 
           "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_gap", 
           "age_w1_sc", "feelings_w1_sc"))

# Missing pattern 
md.pattern(t, rotate.names = TRUE)

# Imputation 
imp <- mice(t, m = 1, seed = 875421)

# Create complete data 
t_imp <- complete(imp, "long", include = TRUE)

# head(t); head(t_imp[t_imp$.imp == 1, ])

# replace missing values in original data set
var <- colnames(t)
for (i in 1:length(var)) {
  
  data[, var[i]] <- t_imp[t_imp$.imp == 1, var[i]]
  
}

# Checking data after imputing values 
# mice::md.pattern(data[, var])
# summary(data$sportPartic_w1)
# head(data[, var])

# Add a scaled version of the mediator to use in outcome model 
data <- data %>%
  mutate(selfEst_w3_sc = as.vector(scale(selfEst_w3)))


## Correlations among PS variables -----------------------------------------

round(cor(data[, c("sex_w1", "white_w1", "black_w1", 
                   "sportPartic_w1", "selfEst_w3", "depress_w4", 
                   "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_gap", 
                   "age_w1_sc", "feelings_w1_sc")]), 
      digits = 2)
#                   sex_w1 white_w1 black_w1 sportPartic_w1 selfEst_w3 depress_w4 familyStruct_w1 parentalEdu_w1_sc healthInsur_gap age_w1_sc feelings_w1_sc
# sex_w1              1.00    -0.03     0.05          -0.10      -0.08       0.12            0.01             -0.06           -0.04     -0.05           0.10
# white_w1           -0.03     1.00    -0.69           0.01      -0.03      -0.09            0.20             -0.06            0.08      0.02          -0.11
# black_w1            0.05    -0.69     1.00           0.00       0.06       0.06           -0.22              0.17           -0.01     -0.01           0.06
# sportPartic_w1     -0.10     0.01     0.00           1.00       0.05      -0.09            0.03              0.14            0.02     -0.08          -0.12
# selfEst_w3         -0.08    -0.03     0.06           0.05       1.00      -0.31            0.05              0.05            0.00      0.01          -0.27
# depress_w4          0.12    -0.09     0.06          -0.09      -0.31       1.00           -0.11             -0.10           -0.01      0.01           0.33
# familyStruct_w1     0.01     0.20    -0.22           0.03       0.05      -0.11            1.00             -0.01            0.10     -0.02          -0.09
# parentalEdu_w1_sc  -0.06    -0.06     0.17           0.14       0.05      -0.10           -0.01              1.00            0.14      0.00          -0.15
# healthInsur_gap    -0.04     0.08    -0.01           0.02       0.00      -0.01            0.10              0.14            1.00     -0.03          -0.02
# age_w1_sc          -0.05     0.02    -0.01          -0.08       0.01       0.01           -0.02              0.00           -0.03      1.00           0.10
# feelings_w1_sc      0.10    -0.11     0.06          -0.12      -0.27       0.33           -0.09             -0.15           -0.02      0.10           1.00
corrplot::corrplot(cor(data[, c("sex_w1", "white_w1", "black_w1", 
                                "sportPartic_w1", "selfEst_w3", "depress_w4", 
                                "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_gap", 
                                "age_w1_sc", "feelings_w1_sc")]),
                   method = c("ellipse"))

PerformanceAnalytics::chart.Correlation(data[, c("sex_w1", "white_w1", "black_w1", 
                                                 "sportPartic_w1", "selfEst_w3", "depress_w4", 
                                                 "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_gap", 
                                                 "age_w1_sc", "feelings_w1_sc")])


## RE PS model stepwise from intercept only model --------------------------------------------------------------

# Run RE PS model & add predictors (with new Other variable)
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + healthInsur_gap + parentalEdu_w1_sc + familyStruct_w1 + 
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) # Model converged 

# Run SL & FE models 
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                white_w1 + black_w1 + healthInsur_gap + parentalEdu_w1_sc + familyStruct_w1",
                family = "binomial", 
                data = data)
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                white_w1 + black_w1 + healthInsur_gap + parentalEdu_w1_sc + familyStruct_w1 + 
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)


### Proceed with Mediator & Outcome models ----------------------------------

#### PS & IPTW ---------------------------------------------------------------
# SL PS model 
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
                family = "binomial", 
                data = data)
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))
# FE PS model 
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))
# RE PS model 
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap +
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) 
data$ps_re <- predict(psmod_re, type = "response")
data$ps_re_logit <- predict(psmod_re, type = "link")
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))


#### Drop Nonoverlap ---------------------------------------------------------

# SL 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_sl < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_sl > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.01)), 
  "): ", sum(I(data$iptw_sl < quantile(data$iptw_sl, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.99)), 
  "): ", sum(I(data$iptw_sl > quantile(data$iptw_sl, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_sl, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < percentiles[1], percentiles[1], 
                          ifelse(iptw_sl > percentiles[2], percentiles[2], 
                                 iptw_sl)))

# FE 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_fe < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_fe > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.01)), 
  "): ", sum(I(data$iptw_fe < quantile(data$iptw_fe, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.99)), 
  "): ", sum(I(data$iptw_fe > quantile(data$iptw_fe, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_fe, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < percentiles[1], percentiles[1], 
                          ifelse(iptw_fe > percentiles[2], percentiles[2], 
                                 iptw_fe)))

# RE
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_re < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_re > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_re, probs = c(0.01)), 
  "): ", sum(I(data$iptw_re < quantile(data$iptw_re, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_re, probs = c(0.99)), 
  "): ", sum(I(data$iptw_re > quantile(data$iptw_re, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_re, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_re = ifelse(iptw_re < percentiles[1], percentiles[1], 
                          ifelse(iptw_re > percentiles[2], percentiles[2], 
                                 iptw_re)))


#### Estimate Effects --------------------------------------------------------

##### Mediator models ---------------------------------------------------------

# Mediator models (SL) 
# Single-Level PS & SL med/outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect PS & SL med/outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_fe
  )
# Random-Effect PS & SL med/outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_re
  )

# Mediator models (FE) 
# Single-Level PS & FE med/outcome
med_slfe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
med_fefe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Random-Effect PS & FE med/outcome
med_refe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_re)

# Mediator models (RE) 
### Add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))
# Single-Level PS & RE med/outcome
med_slre <-
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


##### Outcome models ----------------------------------------------------------

# Outcome models (SL) 
# Single-Level
out_slsl <-
  glm(
    formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_fe
  )
# Random-Effect
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap",
    data = data,
    weights = iptw_re
  )

# Outcome models (FE) 
# Single-Level PS & FE med/outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)
# Random-Effect PS & FE med/outcome
out_refe <-glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + as.factor(CLUSTER2)", 
      data = data,
      weights = iptw_re
)

# Outcome models (RE) 
# If not done so with mediator model, add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))
# Single-Level PS & RE med/outcome
out_slre <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_gap + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


##### Display Estimates -------------------------------------------------------
# Add NDE & NIE estimates 
results_DF <- data.frame(
  cond = c("slre", "fere", "rere"),
  NDE = c(
    summary(out_slre)$coef["sportPartic_w1", "Estimate"], 
    summary(out_fere)$coef["sportPartic_w1", "Estimate"], 
    summary(out_rere)$coef["sportPartic_w1", "Estimate"]
  ),
  
  NIE = c(
    summary(med_slre)$coef["sportPartic_w1", "Estimate"] * summary(out_slre)$coef["selfEst_w3_sc", "Estimate"], 
    summary(med_fere)$coef["sportPartic_w1", "Estimate"] * summary(out_fere)$coef["selfEst_w3_sc", "Estimate"], 
    summary(med_rere)$coef["sportPartic_w1", "Estimate"] * summary(out_rere)$coef["selfEst_w3_sc", "Estimate"]
  )
)
# Display 
results_DF






# Test Dr Liu Solution (with healthInsur_w1) ------------------------------


## Data imputation & scaling  ----------------------------------------------

# Need to scale: 
# Scale, rename, & drop non-scaled variable
data <- data %>%
  mutate(
    age_w1_sc = as.vector(scale(age_w1)),
    parentalEdu_w1_sc = as.vector(scale(parentalEdu_w1)),
    feelings_w1_sc = as.vector(scale(feelings_w1))
  ) %>%
  select(-age_w1, -parentalEdu_w1, -feelings_w1)  # Remove original columns if desired

# Select variables in PS model 
t <- data %>% 
  select(!c(AID, CLUSTER2, sport_w1, n)) %>% 
  # select(!c(AID, CLUSTER2, sport_w1, n, healthInsur_gap))
  # healthInsur_w3, healthInsur_w4, 
  # healthInsur_w1))
  select(c("sex_w1", "white_w1", "black_w1", 
           "sportPartic_w1", "selfEst_w3", "depress_w4", 
           "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_w3", "healthInsur_w4", "healthInsur_w1", 
           "age_w1_sc", "feelings_w1_sc"))

# Missing pattern 
md.pattern(t, rotate.names = TRUE)

# Imputation 
imp <- mice(t, m = 1, seed = 875421)

# Create complete data 
t_imp <- complete(imp, "long", include = TRUE)

# head(t); head(t_imp[t_imp$.imp == 1, ])

# replace missing values in original data set
var <- colnames(t)
for (i in 1:length(var)) {
  
  data[, var[i]] <- t_imp[t_imp$.imp == 1, var[i]]
  
}

# Checking data after imputing values 
# mice::md.pattern(data[, var])
# summary(data$sportPartic_w1)
# head(data[, var])

# Add a scaled version of the mediator to use in outcome model 
data <- data %>%
  mutate(selfEst_w3_sc = as.vector(scale(selfEst_w3)))


## Correlations among PS variables -----------------------------------------

round(cor(data[, c("sex_w1", "white_w1", "black_w1", 
                   "sportPartic_w1", "selfEst_w3", "depress_w4", 
                   "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_w3", "healthInsur_w4", "healthInsur_w1", 
                   "age_w1_sc", "feelings_w1_sc")]), 
      digits = 2)
#                   sex_w1 white_w1 black_w1 sportPartic_w1 selfEst_w3 depress_w4 familyStruct_w1 parentalEdu_w1_sc healthInsur_w3 healthInsur_w4 healthInsur_w1 age_w1_sc feelings_w1_sc
# sex_w1              1.00    -0.02     0.05          -0.10      -0.07       0.12            0.02             -0.03           0.06           0.10           0.00     -0.05           0.10
# white_w1           -0.02     1.00    -0.69           0.01      -0.03      -0.09            0.20             -0.05           0.07           0.02          -0.12      0.02          -0.12
# black_w1            0.05    -0.69     1.00           0.00       0.06       0.06           -0.22              0.16          -0.02          -0.02           0.03     -0.01           0.06
# sportPartic_w1     -0.10     0.01     0.00           1.00       0.05      -0.09            0.04              0.12           0.10           0.09          -0.06     -0.08          -0.11
# selfEst_w3         -0.07    -0.03     0.06           0.05       1.00      -0.31            0.05              0.04           0.11           0.07          -0.05      0.01          -0.27
# depress_w4          0.12    -0.09     0.06          -0.09      -0.31       1.00           -0.11             -0.07          -0.12          -0.10           0.06      0.00           0.33
# familyStruct_w1     0.02     0.20    -0.22           0.04       0.05      -0.11            1.00              0.00           0.09           0.04          -0.16     -0.02          -0.10
# parentalEdu_w1_sc  -0.03    -0.05     0.16           0.12       0.04      -0.07            0.00              1.00           0.13           0.11          -0.20      0.00          -0.10
# healthInsur_w3      0.06     0.07    -0.02           0.10       0.11      -0.12            0.09              0.13           1.00           0.21          -0.18     -0.02          -0.09
# healthInsur_w4      0.10     0.02    -0.02           0.09       0.07      -0.10            0.04              0.11           0.21           1.00          -0.11      0.02          -0.07
# healthInsur_w1      0.00    -0.12     0.03          -0.06      -0.05       0.06           -0.16             -0.20          -0.18          -0.11           1.00      0.04           0.07
# age_w1_sc          -0.05     0.02    -0.01          -0.08       0.01       0.00           -0.02              0.00          -0.02           0.02           0.04      1.00           0.10
# feelings_w1_sc      0.10    -0.12     0.06          -0.11      -0.27       0.33           -0.10             -0.10          -0.09          -0.07           0.07      0.10           1.00
corrplot::corrplot(cor(data[, c("sex_w1", "white_w1", "black_w1", 
                                "sportPartic_w1", "selfEst_w3", "depress_w4", 
                                "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_w3", "healthInsur_w4", "healthInsur_w1", 
                                "age_w1_sc", "feelings_w1_sc")]),
                   method = c("ellipse"))

PerformanceAnalytics::chart.Correlation(data[, c("sex_w1", "white_w1", "black_w1", 
                                                 "sportPartic_w1", "selfEst_w3", "depress_w4", 
                                                 "familyStruct_w1", "parentalEdu_w1_sc", "healthInsur_w3", "healthInsur_w4", "healthInsur_w1", 
                                                 "age_w1_sc", "feelings_w1_sc")])


## RE PS model stepwise from intercept only model --------------------------------------------------------------

# Run RE PS model & add predictors (with new Other variable)
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) # Model converged 

# Run SL & FE models 
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1",
                family = "binomial", 
                data = data)
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)

### Proceed with Mediator & Outcome models ----------------------------------

#### PS & IPTW ---------------------------------------------------------------
# SL PS model 
psmod_sl <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1",
                family = "binomial", 
                data = data)
data$ps_sl <- predict(psmod_sl, type = "response") 
data$ps_sl_logit <- predict(psmod_sl, type = "link")
data <- cbind(data, iptw_sl = with(data, (sportPartic_w1 / ps_sl) + (1 - sportPartic_w1) / (1 - ps_sl)))
# FE PS model 
psmod_fe <- glm(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                as.factor(CLUSTER2)",
                family = "binomial", 
                data = data)
data$ps_fe <- predict(psmod_fe, type = "response")
data$ps_fe_logit <- predict(psmod_fe, type = "link")
data <- cbind(data, iptw_fe = with(data, (sportPartic_w1 / ps_fe) + (1 - sportPartic_w1) / (1 - ps_fe)))
# RE PS model 
psmod_re <- lme4::glmer(formula = "sportPartic_w1 ~ feelings_w1_sc + sex_w1 + age_w1_sc + 
                        white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w1 +
                        (1 | CLUSTER2)",
                        family = "binomial", 
                        data = data) 
data$ps_re <- predict(psmod_re, type = "response")
data$ps_re_logit <- predict(psmod_re, type = "link")
data <- cbind(data, iptw_re = with(data, (sportPartic_w1 / ps_re) + (1 - sportPartic_w1) / (1 - ps_re)))


#### Drop Nonoverlap ---------------------------------------------------------

# SL 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_sl < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_sl > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.01)), 
  "): ", sum(I(data$iptw_sl < quantile(data$iptw_sl, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_sl, probs = c(0.99)), 
  "): ", sum(I(data$iptw_sl > quantile(data$iptw_sl, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_sl, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_sl = ifelse(iptw_sl < percentiles[1], percentiles[1], 
                          ifelse(iptw_sl > percentiles[2], percentiles[2], 
                                 iptw_sl)))

# FE 
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_fe < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_fe > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.01)), 
  "): ", sum(I(data$iptw_fe < quantile(data$iptw_fe, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_fe, probs = c(0.99)), 
  "): ", sum(I(data$iptw_fe > quantile(data$iptw_fe, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_fe, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_fe = ifelse(iptw_fe < percentiles[1], percentiles[1], 
                          ifelse(iptw_fe > percentiles[2], percentiles[2], 
                                 iptw_fe)))

# RE
paste0("Number of PSs < 0.0111: ", sum(I(data$ps_re < 0.0111)), 
       "; Number of PSs > 0.999: ", sum(I(data$ps_re > 0.999)))
paste0(
  "Number of cases < 1st percentile of IPTW (", quantile(data$iptw_re, probs = c(0.01)), 
  "): ", sum(I(data$iptw_re < quantile(data$iptw_re, probs = c(0.01)))), 
  "; Number of cases > 99th percentile of IPTW (", quantile(data$iptw_re, probs = c(0.99)), 
  "): ", sum(I(data$iptw_re > quantile(data$iptw_re, probs = c(0.99))))
)
# Change < 1st percentile to 1st the percentile value & > 99th to the 99th percentile value 
percentiles <- quantile(data$iptw_re, probs = c(0.01, 0.99))
data <- data %>% 
  mutate(iptw_re = ifelse(iptw_re < percentiles[1], percentiles[1], 
                          ifelse(iptw_re > percentiles[2], percentiles[2], 
                                 iptw_re)))


#### Estimate Effects --------------------------------------------------------

##### Mediator models ---------------------------------------------------------

# Mediator models (SL) 
# Single-Level PS & SL med/outcome
med_slsl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect PS & SL med/outcome
med_fesl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_fe
  )
# Random-Effect PS & SL med/outcome
med_resl <-
  glm(
    formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3",
    data = data,
    weights = iptw_re
  )

# Mediator models (FE) 
# Single-Level PS & FE med/outcome
med_slfe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
med_fefe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)

# Random-Effect PS & FE med/outcome
med_refe <- glm(formula = "selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_re)

# Mediator models (RE) 
### Add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))
# Single-Level PS & RE med/outcome
med_slre <-
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
med_fere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
med_rere <- 
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w3 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


##### Outcome models ----------------------------------------------------------

# Outcome models (SL) 
# Single-Level
out_slsl <-
  glm(
    formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_sl
  )
# Fixed-Effect
out_fesl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_fe
  )
# Random-Effect
out_resl <-
  glm(
    "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4",
    data = data,
    weights = iptw_re
  )

# Outcome models (FE) 
# Single-Level PS & FE med/outcome
out_slfe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_sl)
# Fixed-Effect PS & FE med/outcome 
out_fefe <- glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data, 
      weights = iptw_fe)
# Random-Effect PS & FE med/outcome
out_refe <-glm(formula = "depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + as.factor(CLUSTER2)", 
      data = data,
      weights = iptw_re
)

# Outcome models (RE) 
# If not done so with mediator model, add column of just 1s for level-2 weight
data <- cbind(data, L2weight = rep(1, nrow(data)))
# Single-Level PS & RE med/outcome
out_slre <-
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_sl", "L2weight")
  )
# Fixed-Effect PS & RE med/outcome 
out_fere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_fe", "L2weight")
  )
# Random-Effect PS & RE med/outcome
out_rere <- 
  WeMix::mix(
    formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
      white_w1 + black_w1 + 
      parentalEdu_w1_sc + familyStruct_w1 + healthInsur_w4 + (1 | CLUSTER2),
    data = data,
    weights = c("iptw_re", "L2weight")
  )


##### Display Estimates -------------------------------------------------------
# Add NDE & NIE estimates 
results_DF <- data.frame(
  cond = c("slre", "fere", "rere"),
  NDE = c(
    summary(out_slre)$coef["sportPartic_w1", "Estimate"], 
    summary(out_fere)$coef["sportPartic_w1", "Estimate"], 
    summary(out_rere)$coef["sportPartic_w1", "Estimate"]
  ),
  
  NIE = c(
    summary(med_slre)$coef["sportPartic_w1", "Estimate"] * summary(out_slre)$coef["selfEst_w3_sc", "Estimate"], 
    summary(med_fere)$coef["sportPartic_w1", "Estimate"] * summary(out_fere)$coef["selfEst_w3_sc", "Estimate"], 
    summary(med_rere)$coef["sportPartic_w1", "Estimate"] * summary(out_rere)$coef["selfEst_w3_sc", "Estimate"]
  )
)
# Display 
results_DF




















































