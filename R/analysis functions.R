#####################################################
# Program: analysis functions.R
# Programmer: Chase
# Date created: 06.27.2024
#
# Purpose: The program contains all the R functions
# for analyzing the data. The goal of doing this
# is to ensure that these functions only need to be
# editing one time and will apply throughout all
# simulations easily.
#
# Modifications:
# - Removed the MAR only analyses. Instead, 
#   incorporating those into the percentages.
#####################################################

# Test
# tar_data <- tar_read(cohort_data_0_0.7_0.245_.0.2_7_0.7_0.7_Parameters_Abortion07_Preeclampsia07.xlsx_1_10000)


run_analysis <- function(data, rr_abortion, rr_preec, marginal_p_miss_severity,
                         beta12, marginal_p_miss_miscarriage, gamma1,
                         pnc_wk){
  
  hold <- prep_data_for_analysis(data, pnc_wk) %>% 
    conduct_analysis() %>% 
    mutate(
      rr_abortion = rr_abortion,
      rr_preec = rr_preec,
      marginal_p_miss_severity = marginal_p_miss_severity,
      beta12 = beta12,
      marginal_p_miss_miscarriage = marginal_p_miss_miscarriage,
      gamma1 = gamma1,
      pnc_wk = pnc_wk
    )
  
  return(hold)
  
}





#########################################
# FUNCTION: prep_data_for_analysis()
# PURPOSE: Prepare the data for analysis
# functions.

# Checked: 08.09.2024
#########################################

prep_data_for_analysis <- function(data, pnc_wk){
  
  # First, need to do some data cleaning
  hold <- data %>% 
    mutate(
      
      # Make indicator variables for observed delivery
      obs_delivery_mar_mnar = ifelse(ltfu_mar_mnar == 'not' & pregout_t_mar_mnar >= 18,
                                     1,
                                     0),
      
      # Make indicator for observed pregnancy outcome
      obs_outcome_mar_mnar = ifelse(ltfu_mar_mnar == 'not',
                                    1,
                                    0),
      
      
      # Establish the time-to-event outcomes for survival analyses
      ## 0 = censor
      ## 1 = preeclampsia, regardless of preg outcome
      ## 2 = fetal death wo preeclampsia
      ## 3 = live birth wo preeclampsia
      # Calculate times to event
      
      ## Outcome indicator - multivariate, important for AJ estimator
      final_pregout_mar_mnar = case_when(pregout_mar_mnar == 'unknown' ~ 0,
                                         preeclampsia_mar_mnar == 1 ~ 1,
                                         pregout_mar_mnar == 'fetaldeath' ~ 2,
                                         pregout_mar_mnar == 'livebirth' ~ 3),
      
      ## Time to event - Establish the time-to-event for the outcome indicator
      final_pregout_marmnar_tte = ifelse(pregout_mar_mnar == 'unknown',
                                         t_ltfu_mar_mnar - pnc_wk,
                                         pregout_t_mar_mnar - pnc_wk),
      
      ## Add small constant to follow-up time for immediate censors
      final_pregout_marmnar_tte = ifelse(final_pregout_marmnar_tte == 0,
                                         0.0001,
                                         final_pregout_marmnar_tte),
      
      # Sensitivity analyses
      
      ## Outcome (pre-delivery preeclampsia) immediately upon censoring
      ## among all pregnancies LTFU
      preeclampsia_marmnar_immediate_outc_all = ifelse(pregout_mar_mnar == 'unknown',
                                                       1,
                                                       preeclampsia_mar_mnar),
      
      ## Outcome immediately upon censoring only among treated. Assume
      ## that untreated do not have outcome
      preeclampsia_marmnar_immediate_outc_trt = case_when(pregout_mar_mnar == 'unknown' &
                                                            trt == 1 ~ 1,
                                                          pregout_mar_mnar == 'unknown' &
                                                            trt == 0 ~ 0,
                                                          pregout_mar_mnar != 'unknown' ~ preeclampsia_mar_mnar),
      
      ## Outcome immediately upon censoring only among untreated. Assume
      ## that treated do not have outcome
      preeclampsia_marmnar_immediate_outc_untrt = case_when(pregout_mar_mnar == 'unknown' &
                                                              trt == 0 ~ 1,
                                                            pregout_mar_mnar == 'unknown' &
                                                              trt == 1 ~ 0,
                                                            pregout_mar_mnar != 'unknown' ~ preeclampsia_mar_mnar),
      
      ## No outcome and go the full follow-up without outcome. 
      preeclampsia_marmnar_no_outc = ifelse(pregout_mar_mnar == 'unknown',
                                            0,
                                            preeclampsia_mar_mnar)
      
    )
  
  return(hold)
  
}

#########################################
# FUNCTION: conduct_analysis()
# PURPOSE: Conduct the analyses for the
# dataset.
#########################################

conduct_analysis <- function(data){
  
  # Get the severity distribution for each simulation
  sev_dist <- get_severity_dist(data)
  
  # Get the risks from the potential outcomes
  potential_risks <- calculate_pot_risks(data)
  
  # Get the risks among observed deliveries

    observed_deliveries_mar_mnar <- data %>% 
      filter(obs_delivery_mar_mnar == 1) %>% 
      calculate_risks(sev_dist, preeclampsia_pre_miss)
  
  # Get the risks among observed pregnancy outcomes
    
    observed_outcomes_mar_mnar <- data %>% 
      filter(obs_outcome_mar_mnar == 1) %>% 
      calculate_risks(sev_dist, preeclampsia_pre_miss)
  
  # Time-to-event analyses
    
    tte_mar_mnar <- calculate_aj_risks(data, 
                                       sev_dist)
    
  # Sensitivity analyses
    
    ## Assume that all missing outcomes have an outcome
    
      sens_anal_mar_mnar_all_outc <- calculate_risks(data, sev_dist,
                                                     preeclampsia_marmnar_immediate_outc_all)
      
    ## Assume that treated missing have an outcome
      sens_anal_mar_mnar_trt_outc <- calculate_risks(data, sev_dist,
                                                     preeclampsia_marmnar_immediate_outc_trt)
      
    ## Assume that untreated missing have an outcome
      sens_anal_mar_mnar_untrt_out <- calculate_risks(data, sev_dist,
                                                      preeclampsia_marmnar_immediate_outc_untrt)
      
    ## Assume that all missing outcomes do not have an outcome
      
      sens_anal_mar_mnar_no_outc <- calculate_risks(data, sev_dist,
                                                    preeclampsia_marmnar_no_outc)

  
  return(tibble(
    sev_dist = list(sev_dist),
    potential_risks = list(potential_risks),
    observed_deliveries = list(observed_deliveries_mar_mnar),
    observed_outcomes = list(observed_outcomes_mar_mnar),
    tte = list(tte_mar_mnar),
    sens_anal_all_outc = list(sens_anal_mar_mnar_all_outc),
    sens_anal_trt_outc = list(sens_anal_mar_mnar_trt_outc),
    sens_anal_untrt_out = list(sens_anal_mar_mnar_untrt_out),
    sens_anal_no_outc = list(sens_anal_mar_mnar_no_outc)
  ))

}

 




##############################################
# FUNCTION: get_severity_dist()
# PURPOSE: The purpose of this function is to 
# get the distribution of patients by 
# disease severity at baseline.
# Target population: ALL patients.
##############################################

get_severity_dist <- function(data){
  
  sev_dist <- data %>%
    group_by(severity) %>%
    summarize(prop = n() / nrow(data),
              .groups = 'drop')
  
  return(sev_dist)
  
}





##############################################
# FUNCTION: calculate_pot_risks()
# PURPOSE: The purpose of this function is to 
# calculate the risk of preeclampsia using
# the potential outcomes. 
##############################################

calculate_pot_risks <- function(dataset){
  
  # Calculate the risks using the potential outcomes
  pot_out_risks <- dataset %>% 
    summarize(
      risk0 = sum(preeclampsia0) / n(),
      risk1 = sum(preeclampsia1) / n(),
      rd = risk1 - risk0,
      rr = risk1 / risk0
    )

  return(pot_out_risks)
  
}




##############################################
# FUNCTION: calculate_risks()
# PURPOSE: The purpose of this function to 
# calculate the risk of preeclampsia in a 
# closed cohort without right censoring. Only
# one obs per person.

# This implements non-parametric direct
# standardization by the distribution of severity
# at baseline to estimate the population
# average treatment effect (ATE).
##############################################

# Revised so that I can call the outcome variable for calculating risks in the function.

calculate_risks <- function(dataset, sev_dist, risk_var){
  
  risk_var <- enquo(risk_var)
  
  # Calculate the risks within strata of treatment and severity
  strat_risks <- dataset %>% 
    group_by(trt, severity) %>% 
    summarize(
      risk = sum(!!risk_var) / n(),
      .groups = 'drop'
    )
  
  # Merge the severity distribution and multiply through to calculate the 
  # standardized risks within each trt strata.
  merge_sev_prop <- left_join(strat_risks, sev_dist, by = c("severity" = "severity")) %>% 
    mutate(std_risk = risk * prop) %>% 
    group_by(trt) %>% 
    summarize(risk = sum(std_risk),
              .groups = 'drop')
  
  strata_risks_wide <- merge_sev_prop %>% 
    pivot_wider(names_from = trt, 
                values_from = risk, 
                names_prefix = "risk")
  
  risks <- strata_risks_wide %>% 
    ungroup() %>% 
    mutate(rd = risk1 - risk0,
           rr = risk1 / risk0)
    
  
  return(risks)
}











##############################################
# FUNCTION: aj_estimator()
# PURPOSE: The purpose of this function is to 
# implement the AJ estimator to estimate 
# risks.

# This implements non-parametric direct
# standardization by the distribution of severity
# at baseline to estimate the population
# average treatment effect (ATE).
##############################################

aj_estimator <- function(data_subset) { #outcome_ind, time_var
  
  # Run the AJ model
  aj <- survfit(Surv(time, outcome) ~ trt, data = data_subset)

  mod <- summary(aj)
  
  summod <- data.frame(t = mod$time,
                       r = mod$pstate[, 2], 
                       trt = c(rep(0, length(mod[["strata"]][mod[["strata"]] == "trt=0"])), 
                               rep(1, length(mod[["strata"]][mod[["strata"]] == "trt=1"])))
  ) %>%
    # filter(t < 43 + 0.5) %>%  # Deal with the jittering of outcomes
    group_by(trt) %>% 
    summarize(risk = last(r),
              .groups = 'drop') %>% 
    pivot_wider(names_from = trt,
                values_from = risk,
                names_glue = "{.value}{trt}") %>% 
    select(risk0, risk1)
  
  return(summod)
}

           



##############################################
# FUNCTION: calculate_aj_risks()
# PURPOSE: The purpose of this function is to 
# calculate the risk of preeclampsia in a 
# cohort with right censoring using an AJ
# estimator.

# This implements non-parametric direct
# standardization by the distribution of severity
# at baseline to estimate the population
# average treatment effect (ATE).
##############################################

## Revised final_pregout_mar_mnar final_pregout_mar_mnar_tte

calculate_aj_risks <- function(dataset, sev_dist) { #, outcome_ind, time_var
  
  # outcome_ind <- enquo(outcome_ind)
  # time_var <- enquo(time_var)

  # Convert the dataset to a data.table
  mar <- as.data.table(dataset)
  
  # Calculate the count for each group
  mar[, count := .N, by = .(trt, final_pregout_marmnar_tte)]

  # Jitter event times and calculate the adjusted time and outcome
  mar[, jitter := runif(.N, min = -.01, max = .01)]
  mar[, time := ifelse(count > 1, final_pregout_marmnar_tte + jitter, final_pregout_marmnar_tte)]
  mar[, outcome := as.factor(final_pregout_mar_mnar)]
  

  results <- mar[, aj_estimator(.SD), by = severity] %>%
    as.tibble() %>%
    # Convert severity to its original data types if necessary
    mutate(severity = as.numeric(severity)) %>%
    # Left merge the severity distribution onto the risks
    left_join(sev_dist, by = c("severity" = "severity")) %>%
    mutate(risk0_std = risk0 * prop,
           risk1_std = risk1 * prop) %>%
    # Summarize the risks within strata of PNC week
    summarize(risk0 = sum(risk0_std),
              risk1 = sum(risk1_std)) %>%
    mutate(rd = risk1 - risk0,
           rr = risk1 / risk0)
  
  return(results)
}
