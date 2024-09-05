# Display the final dataframe with PS model & Mediator/Outcome labels
print(results_DF)
#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 1     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050 -0.3007998 0.09570006 -0.5084805  0.004777044 -0.3077056 0.10218566 -0.5085235  0.006911447  Fixed-Effect                     Fixed-Effect
# 2     fesl -0.2558156 -0.2563010 -0.09886261 -0.0261259555 -0.2921923 0.08181015 -0.5118420 -0.010035295 -0.2961822 0.08569538 -0.5117962 -0.011418954  Fixed-Effect                     Single-Level
# 3     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374 -0.2797301 0.09890207 -0.5069342 -0.000239732 -0.2829443 0.10403690 -0.5077968  0.002558890 Random-Effect                     Fixed-Effect
# 4     resl -0.2959685 -0.2961819 -0.06592642  0.0012385894 -0.2334787 0.10900586 -0.5444086 -0.064872212 -0.2386974 0.11433888 -0.5457297 -0.064180945 Random-Effect                     Single-Level
# 5     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906 -0.2637974 0.10503455 -0.4788932  0.010653283 -0.2650747 0.10698408 -0.4789265  0.014197515  Single-Level                     Fixed-Effect
# 6     slsl -0.3317336 -0.3318426 -0.01230346  0.0247090071 -0.1747149 0.16068737 -0.5693915 -0.084157900 -0.1768113 0.16116459 -0.5715905 -0.083434503  Single-Level                     Single-Level
# 7     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.3304645  0.090037440  Fixed-Effect                    Random-Effect
# 8  fere_cm -0.2425750 -0.2431934 -0.10177153 -0.0098265525 -0.3116530 0.08420532 -0.5122988 -0.004824141 -0.3299738 0.08873388 -0.3299738  0.088733884  Fixed-Effect Random-Effect with Cluster Means
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.2862409  0.098414294 Random-Effect                    Random-Effect
# 10 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.2981107  0.093321648 Random-Effect Random-Effect with Cluster Means
# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.2486733  0.129917535  Single-Level                    Random-Effect
# 12 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.2800797  0.104324713  Single-Level Random-Effect with Cluster Means




#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model

# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.2486733  0.129917535  Single-Level                    Random-Effect
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.2862409  0.098414294 Random-Effect                    Random-Effect
# 7     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.3304645  0.090037440  Fixed-Effect                    Random-Effect


# > quantile(slre_ci_TNIE_DF$direct, probs = c(0.025, 0.975))
# 2.5%       97.5% 
#   -0.51217804 -0.04402773 
# > quantile(rere_ci_TNIE_DF$direct, probs = c(0.025, 0.975)) # RERE PNDE 
# 2.5%       97.5% 
#   -0.52461418 -0.03481273 
# > quantile(fere_ci_TNIE_DF$direct, probs = c(0.025, 0.975)) # FERE PNDE 
# 2.5%        97.5% 
#   -0.512662836 -0.001785989 



# are PNDE CIs for RE & RE-Mean med/outcome models the same ? 
#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.2486733  0.129917535  Single-Level                    Random-Effect
# 12 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108 -0.2777162 0.10351657 -0.4686721 -0.003627143 -0.2800797 0.10432471 -0.2800797  0.104324713  Single-Level Random-Effect with Cluster Means
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.2862409  0.098414294 Random-Effect                    Random-Effect
# 10 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984 -0.2902577 0.09149804 -0.5107439 -0.015202755 -0.2981107 0.09332165 -0.2981107  0.093321648 Random-Effect Random-Effect with Cluster Means
# Answer: They are not identical but similar


# Check differences in PS truncation's effects on estimates  --------------
# Change PS truncation back to prior values & obtain estimates. Then compare these estimates to both old estimates & current (recently ran) estimates 

## New Run 
# RE PS Model did not converge 
# PS Truncation 
  # SL
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"
    # [1] "Number of cases < 1st percentile of IPTW (1.25694007958331): 32; Number of cases > 99th percentile of IPTW (4.323208231669): 32"
  # FE 
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 11"
    # [1] "Number of cases < 1st percentile of IPTW (1.02488369152563): 32; Number of cases > 99th percentile of IPTW (6.72733775692784): 32"
  # RE 
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"
    # [1] "Number of cases < 1st percentile of IPTW (1.10316815375557): 32; Number of cases > 99th percentile of IPTW (5.02691967023504): 32"
# Estimates 
#       cond       TNDE       PNDE        PNIE          TNIE
# 1     slsl -0.3317336 -0.3318426 -0.01230346  0.0247090071
# 2     fesl -0.2558156 -0.2563010 -0.09886261 -0.0261259555
# 3     resl -0.2959685 -0.2961819 -0.06592642  0.0012385894
# 4     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906
# 5     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050
# 6     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374
# 7     slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456
# 8     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596
# 10 slre_cm -0.2203850 -0.2207466 -0.07965423 -0.0397668108
# 11 fere_cm -0.2425750 -0.2431934 -0.10177153 -0.0098265525
# 12 rere_cm -0.2470984 -0.2476017 -0.08812400 -0.0065533984

## Old Run 
# RE PS Model did not converge 
# PS Truncation 
  # SL 
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"
    # [1] "Number of cases < 1st percentile of IPTW (1.26087144090849): 35; Number of cases > 99th percentile of IPTW (5.21931640120409): 4"
  # FE 
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 11"
    # [1] "Number of cases < 1st percentile of IPTW (1.02588244458161): 35; Number of cases > 99th percentile of IPTW (24.1939199626553): 4"
  # RE 
    # [1] "Number of PSs < 0.01: 0; Number of PSs > 0.99: 0"
    # [1] "Number of cases < 1st percentile of IPTW (1.10667547931046): 35; Number of cases > 99th percentile of IPTW (10.5038136824577): 4"
# Estimates 
#       cond       TNDE       PNDE        PNIE         TNIE
# 1     slsl -0.3273004 -0.3274150 -0.01279479  0.023831306
# 2     fesl -0.2022917 -0.2038619 -0.03144804  0.057169062
# 3     resl -0.2705579 -0.2712985 -0.03156782  0.048868309
# 4     slfe -0.2211653 -0.2216441 -0.08317793 -0.034546268
# 5     fefe -0.2230344 -0.2247051 -0.03952605  0.076516517
# 6     refe -0.2272425 -0.2285758 -0.04929975  0.054256810
# 7     slre -0.2628031 -0.2630943 -0.05248441 -0.009126956
# 8     fere -0.2196004 -0.2211833 -0.03915470  0.064140810
# 9     rere -0.2460845 -0.2471651 -0.03843323  0.055476533
# 10 slre_cm -0.2162055 -0.2165688 -0.07977256 -0.040393256
# 11 fere_cm -0.2180661 -0.2195907 -0.03995004  0.058505878
# 12 rere_cm -0.2269542 -0.2281142 -0.04609519  0.043167499





# Re-obtain estimates & CIs by hand & compare results ---------------------
# Re-obtain estimate & CIs for a couple conditions by hand (not the functions & compare results)

## Estimates ---------------------------------------------------------------

### FEFE --------------------------------------------------------------------
summary(med_fefe)$coef["sportPartic_w1", "Estimate"] * summary(out_fefe)$coef["selfEst_w3_sc", "Estimate"] # PNIE
# summary(out_fefe_interac)$coef["sportPartic_w1", "Estimate"] + summary(out_fefe_interac)$coef["selfEst_w3_sc:sportPartic_w1", "Estimate"] # 
summary(out_fefe)$coef["sportPartic_w1", "Estimate"] # TNDE
summary(out_fefe_interac)$coef["sportPartic_w1", "Estimate"] # PNDE 
# by hand _ oriignal estimates 
# PNIE = -0.03952605 != -0.09660659
# TNDE = -0.2230344 != -0.2439341
# PNDE = -0.2247051 != -0.2445695

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 1     fefe -0.2439341 -0.2445695 -0.09660659  0.0126441050 -0.3007998 0.09570006 -0.5084805  0.004777044 -0.3077056 0.10218566 -0.5085235  0.006911447  Fixed-Effect                     Fixed-Effect

# FEFE condition have different estimates 

### FERE --------------------------------------------------------------------
summary(med_fere)$coef["sportPartic_w1", "Estimate"] * summary(out_fere)$coef["selfEst_w3_sc", "Estimate"]
summary(out_fere)$coef["sportPartic_w1", "Estimate"]
summary(out_fere_interac)$coef["sportPartic_w1", "Estimate"] 
# by hand _ original estimates 
# PNIE = -0.0391547 != -0.09700125
# TNDE = -0.2196004 != -0.2497230
# PNDE = -0.2211833 != -0.2503269

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 7     fere -0.2497230 -0.2503269 -0.09700125 -0.0005990211 -0.3119372 0.08752047 -0.5126724 -0.003286732 -0.3304645 0.09003744 -0.3304645  0.090037440  Fixed-Effect                    Random-Effect

# FERE condition have different estimates 

### RERE --------------------------------------------------------------------
summary(med_rere)$coef["sportPartic_w1", "Estimate"] * summary(out_rere)$coef["selfEst_w3_sc", "Estimate"]
summary(out_rere)$coef["sportPartic_w1", "Estimate"]
summary(out_rere_interac)$coef["sportPartic_w1", "Estimate"] 
# by hand _ original estimates 
# PNIE = -0.07559626 == -0.07559626
# TNDE = -0.2682174 == -0.2682174
# PNDE = -0.2686312 == -0.2686312

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 9     rere -0.2682174 -0.2686312 -0.07559626  0.0099211596 -0.2724058 0.09636612 -0.5241173 -0.035136923 -0.2862409 0.09841429 -0.2862409  0.098414294 Random-Effect                    Random-Effect

# RERE condition have identical estimates 

### REFE --------------------------------------------------------------------
summary(med_refe)$coef["sportPartic_w1", "Estimate"] * summary(out_refe)$coef["selfEst_w3_sc", "Estimate"]
summary(out_refe)$coef["sportPartic_w1", "Estimate"]
summary(out_refe_interac)$coef["sportPartic_w1", "Estimate"] 
# by hand _ original estimates 
# PNIE = -0.08816393 == -0.08816393
# TNDE = -0.2464966 == -0.2464966
# PNDE = -0.2471019 == -0.2471019

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 3     refe -0.2464966 -0.2471019 -0.08816393  0.0085269374 -0.2797301 0.09890207 -0.5069342 -0.000239732 -0.2829443 0.10403690 -0.5077968  0.002558890 Random-Effect                     Fixed-Effect

# REFE condition have identical estimates 

### SLRE --------------------------------------------------------------------
summary(med_slre)$coef["sportPartic_w1", "Estimate"] * summary(out_slre)$coef["selfEst_w3_sc", "Estimate"]
summary(out_slre)$coef["sportPartic_w1", "Estimate"]
summary(out_slre_interac)$coef["sportPartic_w1", "Estimate"] 
# by hand _ original estimates 
# PNIE = -0.05248441 != -0.05197434
# TNDE = -0.2628031 != -0.2672821
# PNDE = -0.2630943 != -0.2675694

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 11    slre -0.2672821 -0.2675694 -0.05197434 -0.0080550456 -0.2450241 0.12480705 -0.5120571 -0.043267219 -0.2486733 0.12991754 -0.2486733  0.129917535  Single-Level                    Random-Effect

# While not the exactly the same, like RERE, SLRE are close unlike the other two conditions. 

### SLFE --------------------------------------------------------------------
summary(med_slfe)$coef["sportPartic_w1", "Estimate"] * summary(out_slfe)$coef["selfEst_w3_sc", "Estimate"]
summary(out_slfe)$coef["sportPartic_w1", "Estimate"]
summary(out_slfe_interac)$coef["sportPartic_w1", "Estimate"] 
# by hand _ original estimates 
# PNIE = -0.08332282 != -0.08332282
# TNDE = -0.2253132 != -0.2253132
# PNDE = -0.2257916 != -0.2257916

#       cond       TNDE       PNDE        PNIE          TNIE    PNIE_LL    PNIE_UL    TNDE_LL      TNDE_UL    TNIE_LL    TNIE_UL    PNDE_LL      PNDE_UL            PS                            Model
# 5     slfe -0.2253132 -0.2257916 -0.08332282 -0.0341267906 -0.2637974 0.10503455 -0.4788932  0.010653283 -0.2650747 0.10698408 -0.4789265  0.014197515  Single-Level                     Fixed-Effect

# SLFE condition have identical estimates 

### Report ------------------------------------------------------------------

# FEFE & FERE have different 
# RERE & REFE & SLFE have identical 
# SLRE have similar 

## CI ---------------------------------------------------------------

### FEFE --------------------------------------------------------------------



# Double check convergence for RE & RE-Mean mediator & outcome mod --------
# Double check RE & RE-Mean mediator & outcome models converged for estimates 



# Double check model formulas/statements in functions for CIs -------------

# TNIE & PNDE Estimate formula 
formula = depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + 
  white_w1 + black_w1 + parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),


## bootstrap_ci_re_paral_2() -----------------------------------------------

# Function 
# Fit the random effects models for the bootstrap sample
mediator_rere <- tryCatch({
  WeMix::mix(
    formula = selfEst_w3 ~ sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + # 
      parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2),  
    data = data_boot, 
    weights = c(iptw_str, "L2weight")
  )
}, error = function(e) NULL)  # Handle errors gracefully

outcome_formula <- if (effect_type == "TNIE") {
  depress_w4 ~ selfEst_w3_sc * sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + # TNIE & PNDE formula 
    parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
} else {
  depress_w4 ~ selfEst_w3_sc + sportPartic_w1 + age_w1_sc + sex_w1 + white_w1 + black_w1 + # PNIE & TNDE formula 
    parentalEdu_w1_sc + familyStruct_w1 + (1 | CLUSTER2)
}


## bootstrap_ci_paral_2() --------------------------------------------------


## bootstrap_ci_re_mean_paral() --------------------------------------------








# Check on nonconvergence for RE PS model  --------------------------------
# Check on nonconvergence for RE PS model 








