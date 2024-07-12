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
#####################################################

# library(tidyverse)
# library(survival)
# library(rlang)
# 
# ## Upload the dataset for testing
# #### Missing: Beta1 = 0.01, Gamma0 = 0.1, Gamma1 = 0.001
# data <- readRDS("ab08preec08_beta001_gamma01_001.rds") %>% 
#   mutate(scenario = "test1")
# 



# Calculate the risks among cohorts with no right-censoring.

# pot_risks <- calculate_pot_risks(data1, pnc_dist)
# 
# risks_obsdel_mar <- calculate_risks(subset(data1, obs_delivery_mar == 1), pnc_dist, sev_dist)
# risks_obsdel_mar_mnar <- calculate_risks(subset(data1, obs_delivery_mar_mnar == 1), pnc_dist, sev_dist)
# 
# risks_obsout_mar <- calculate_risks(subset(data1, obs_outcome_mar == 1), pnc_dist, sev_dist)
# risks_obsout_mar_mnar <- calculate_risks(subset(data1, obs_outcome_mar_mnar == 1), pnc_dist, sev_dist)

# # Time-to-event analyses
# tte_mar <- calculate_aj_risks(dataset, pnc_dist, sev_dist, "final_pregout_mar_tte")
# tte_mnar <- calculate_aj_risks(dataset, pnc_dist, sev_dist, "final_pregout_marmnar_tte")


run_analysis <- function(data){
  
  hold <- lapply(data, prep_data_for_analysis)
  
  hold2 <- conduct_analysis(hold)
  
  return(hold2)
  
}





#########################################
# FUNCTION: prep_data_for_analysis()
# PURPOSE: Prepare the data for analysis
# functions.
#########################################

prep_data_for_analysis <- function(data){
  
  # First, need to do some data cleaning
  hold <- data %>% 
    mutate(
      
      # Make indicator variables for observed delivery
      obs_delivery_mar = ifelse(ltfu_mar == 'not' & pregout_t_mar >= 20,
                                1,
                                0),
      obs_delivery_mar_mnar = ifelse(ltfu_mar_mnar == 'not' & pregout_t_mar_mnar >= 20,
                                     1,
                                     0),
      
      # Make indicator for observed pregnancy outcome
      obs_outcome_mar = ifelse(ltfu_mar == 'not',
                               1,
                               0),
      obs_outcome_mar_mnar = ifelse(ltfu_mar_mnar == 'not',
                                    1,
                                    0),
      
      
      # Establish the time-to-event outcomes for survival analyses
      ## 0 = censor
      ## 1 = preeclampsia, regardless of preg outcome
      ## 2 = fetal death wo preeclampsia
      ## 3 = live birth wo preeclampsia
      # Calculate times to event
      
      # MAR 
      ## Outcome indicator
      final_pregout_mar = case_when(pregout_mar == 'unknown' ~ 0,
                                    preeclampsia_mar == 1 ~ 1,
                                    pregout_mar == 'fetaldeath' ~ 2,
                                    pregout_mar == 'livebirth' ~ 3),
      
      ## Time to event 
      final_pregout_mar_tte = ifelse(pregout_mar == 'unknown',
                                     t_ltfu_mar,
                                     pregout_t_mar),
      
      ## Add small constant to follow-up time for immediate censors - CHECK THIS.
      final_pregout_mar_tte = ifelse(final_pregout_mar_tte == 0,
                                     0.1,
                                     final_pregout_mar_tte),
      
      # MAR+MNAR
      ## Outcome indicator
      final_pregout_mar_mnar = case_when(pregout_mar_mnar == 'unknown' ~ 0,
                                         preeclampsia_mar == 1 ~ 1,
                                         pregout_mar == 'fetaldeath' ~ 2,
                                         pregout_mar == 'livebirth' ~ 3),
      
      ## Time to event 
      final_pregout_marmnar_tte = ifelse(pregout_mar == 'unknown',
                                         t_ltfu_mar_mnar + 1,
                                         pregout_t_mar_mnar + 1),
      
      ## Add small constant to follow-up time for immediate censors - CHECK THIS.
      final_pregout_marmnar_tte = ifelse(final_pregout_marmnar_tte == 0,
                                         0.1,
                                         final_pregout_marmnar_tte),
      
      # Sensitivity analysis
      
      ## Outcome immediately upon censoring
      preeclampsia_mar_immediate_outc = ifelse(pregout_mar == 'unknown',
                                               1,
                                               preeclampsia_mar),
      preeclampsia_marmnar_immediate_outc = ifelse(pregout_mar == 'unknown',
                                                   1,
                                                   preeclampsia_mar_mnar),
      
      ## No outcome and go the full follow-up without outcome. 
      preeclampsia_mar_no_outc = ifelse(pregout_mar == 'unknown',
                                        0,
                                        preeclampsia_mar),
      preeclampsia_marmnar_no_outc = ifelse(pregout_mar_mnar == 'unknown',
                                           0,
                                           preeclampsia_mar_mnar)
      
    )
  
  return(hold)
  
  # hold2 <- split(hold,hold$sim_id)
  # 
  # return(hold2)
  
}

#########################################
# FUNCTION: conduct_analysis()
# PURPOSE: Conduct the analyses for the
# dataset.
#########################################

conduct_analysis <- function(data){
  
  # Get the pnc distribution for each simulation
  pnc_dist <- lapply(data, get_pnc_dist)
  
  # Get the severity distribution for each simulation
  sev_dist <- lapply(data, get_severity_dist)
  
  # Get the risks from the potential outcomes
  potential_risks <- mapply(calculate_pot_risks, data, pnc_dist, SIMPLIFY = FALSE)
  
  # Get the risks among observed deliveries
  
    ## MAR:
    observed_deliveries_mar <- mapply(function(data, pnc_dist, sev_dist) {
      #data_subset <- subset(data, obs_delivery_mar == 1)
      data %>% 
        filter(obs_delivery_mar == 1) %>% 
        calculate_risks(pnc_dist, sev_dist, preeclampsia_pre_miss)
    }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
    
    ## MNAR:
    observed_deliveries_mar_mnar <- mapply(function(data, pnc_dist, sev_dist) {
      #data_subset <- subset(hold2, obs_delivery_mar_mnar == 1)
      data %>% 
        filter(obs_delivery_mar == 1) %>% 
        calculate_risks(pnc_dist, sev_dist, preeclampsia_pre_miss)
    }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
  
  # Get the risks among observed pregnancy outcomes
    
    ## MAR
    observed_outcomes_mar <- mapply(function(data, pnc_dist, sev_dist) {
      #data_subset <- subset(data, obs_outcome_mar == 1)
      data %>% 
        filter(obs_delivery_mar == 1) %>% 
        calculate_risks(pnc_dist, sev_dist, preeclampsia_pre_miss)
    }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
    
    ## MNAR
    observed_outcomes_mar_mnar <- mapply(function(data, pnc_dist, sev_dist) {
      #data_subset <- subset(hold2, obs_outcome_mar_mnar == 1)
      data %>% 
        filter(obs_delivery_mar_mnar == 1) %>% 
        calculate_risks(pnc_dist, sev_dist, preeclampsia_pre_miss)
    }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
  
  # Time-to-event analyses
    
    ## MAR
    # tte_mar <- mapply(calculate_aj_risks, data, pnc_dist, sev_dist, 
    #                   MoreArgs = list(outcome_ind = final_pregout_mar, time_var = final_pregout_mar_tte), 
    #                   SIMPLIFY = FALSE) 
    tte_mar <- mapply(function(data_subset, pnc_dist, sev_dist) {
      calculate_aj_risks(data_subset, pnc_dist, sev_dist, final_pregout_mar, final_pregout_mar_tte)
    }, 
    data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
    
    ## MNAR
    # tte_mar_mnar <- mapply(calculate_aj_risks, data, pnc_dist, sev_dist, 
    #                        MoreArgs = list(time_var = "final_pregout_marmnar_tte"), 
    #                        SIMPLIFY = FALSE)
    tte_mar_mnar <- mapply(function(data_subset, pnc_dist, sev_dist) {
      calculate_aj_risks(data_subset, pnc_dist, sev_dist, final_pregout_mar_mnar, final_pregout_marmnar_tte)
    }, 
    data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
    
  # Sensitivity analyses
    
    ## Assume that all missing outcomes have an outcome
    
      ### MAR:
      sens_anal_mar_all_outc <- mapply(function(data, pnc_dist, sev_dist) {
        data %>% 
          calculate_risks(pnc_dist, sev_dist, preeclampsia_mar_immediate_outc)
      }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
      
      ### MNAR:
      sens_anal_mar_mnar_all_outc <- mapply(function(data, pnc_dist, sev_dist) {
        data %>% 
          calculate_risks(pnc_dist, sev_dist, preeclampsia_marmnar_immediate_outc)
      }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
      
    ## Assume that all missing outcomes do not have an outcome
      
      ## MAR:
      sens_anal_mar_no_outc <- mapply(function(data, pnc_dist, sev_dist) {
        data %>% 
          calculate_risks(pnc_dist, sev_dist, preeclampsia_mar_no_outc)
      }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
      
      ## MNAR:
      sens_anal_mar_mnar_no_outc <- mapply(function(data, pnc_dist, sev_dist) {
        data %>% 
          calculate_risks(pnc_dist, sev_dist, preeclampsia_marmnar_no_outc)
      }, data, pnc_dist, sev_dist, SIMPLIFY = FALSE)
  
  # return(data)
  return(tibble(
    pnc_dist = pnc_dist,
    sev_dist = sev_dist,
    potential_risks = potential_risks,
    observed_deliveries_mar = observed_deliveries_mar,
    observed_deliveries_mar_mnar = observed_deliveries_mar_mnar,
    observed_outcomes_mar = observed_outcomes_mar,
    observed_outcomes_mar_mnar = observed_outcomes_mar_mnar,
    tte_mar = tte_mar,
    tte_mar_mnar = tte_mar_mnar,
    sens_anal_mar_all_outc = sens_anal_mar_all_outc,
    sens_anal_mar_mnar_all_outc = sens_anal_mar_mnar_all_outc,
    sens_anal_mar_no_outc = sens_anal_mar_no_outc,
    sens_anal_mar_mnar_no_outc = sens_anal_mar_mnar_no_outc
  ))

}

 



##############################################
# FUNCTION: get_pnc_dist()
# PURPOSE: The purpose of this function is to 
# get the distribution of patients by the
# gestational week of their first pnc 
# encounter.
##############################################

get_pnc_dist <- function(data){
  
  pnc_dist <- (data %>%
                 group_by(pnc_wk) %>%
                 summarize(prop = n() / nrow(data)))$prop
  
  return(pnc_dist)
  
}



##############################################
# FUNCTION: get_severity_dist()
# PURPOSE: The purpose of this function is to 
# get the distribution of patients by 
# disease severity at baseline.
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

calculate_pot_risks <- function(dataset, pnc_dist){
  
  # Calculate the risks using the potential outcomes
  pot_out_risks_strat <- dataset %>% 
    group_by(pnc_wk)%>% 
    summarize(
      risk0 = sum(preeclampsia0) / n(),
      risk1 = sum(preeclampsia1) / n(),
      rd = risk1 - risk0,
      rr = risk1 / risk0
    )
  
  # Calculate the risks standardized to the PNC dist at baseline
  risk0 <- sum(pot_out_risks_strat$risk0 * pnc_dist)
  risk1 <- sum(pot_out_risks_strat$risk1 * pnc_dist)
  
  pot_out_risks_overall <- tibble(
    pnc_wk = 'all',
    risk0 = risk0,
    risk1 = risk1,
    rd = risk1 - risk0,
    rr = risk1 / risk0
  )
  
  pot_out_risks <- rbind(pot_out_risks_strat, pot_out_risks_overall)
  
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

calculate_risks <- function(dataset, pnc_dist, sev_dist, risk_var){
  
  risk_var <- enquo(risk_var)
  
  # Ensure pnc_wk is character type
  dataset <- dataset %>% 
    mutate(pnc_wk = as.character(pnc_wk))
  
  # Calculate the risks using the potential outcomes
  strat_risks <- dataset %>% 
    group_by(trt, pnc_wk, severity) %>% 
    summarize(
      risk = sum(!!risk_var) / n(),
      .groups = 'drop'
    )
  
  # Merge the severity distribution among everyone at baseline on
  # and multiply through in order to calculate the standardized risks within 
  # each trt and pnc_wk strata.
  merge_sev_prop <- left_join(strat_risks, sev_dist, by = c("severity" = "severity")) %>% 
    rowwise() %>% 
    mutate(std_risk = risk * prop) %>% 
    group_by(trt, pnc_wk) %>% 
    summarize(risk = sum(std_risk),
              .groups = 'drop')
  
  strata_risks_wide <- merge_sev_prop %>% 
    pivot_wider(names_from = trt, values_from = risk, names_prefix = "risk")
  
  # Calculate the risks standardized to the PNC dist at baseline
  risk0 <- sum(strata_risks_wide$risk0 * pnc_dist)
  risk1 <- sum(strata_risks_wide$risk1 * pnc_dist)
  
  risks_overall <- tibble(
    pnc_wk = 'all',
    risk0 = risk0,
    risk1 = risk1
  )
  
  risks <- bind_rows(strata_risks_wide, risks_overall) %>% 
    rowwise() %>% 
    mutate(rd = risk1 - risk0,
           rr = risk1 / risk0) %>% 
    ungroup()
  
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

## REvised

aj_estimator <- function(data_subset, outcome_ind, time_var) {
  outcome_ind <- enquo(outcome_ind)
  time_var <- enquo(time_var)
  
  # Jitter ties
  # data_subset <- data_subset %>% 
  #   group_by(!!time_var) %>% 
  #   add_tally(name = "count") %>% 
  #   ungroup() %>% 
  #   mutate(
  #     # Jitter event times
  #     jitter = runif(n(), min = -.01, max = .01),
  #     time = ifelse(count > 1, !!time_var + jitter, !!time_var)
  #   )
  
  # Convert outcome to factor
  data_subset <- data_subset %>% 
    mutate(outcome = as.factor(!!outcome_ind))
  
  # Run the AJ model
  aj <- survfit(Surv(time, outcome) ~ trt, data = data_subset)
  
  mod <- summary(aj)
  
  summod <- data.frame(t = mod$time,
                       r = mod$pstate[, 2], 
                       se = mod$std.err[, 2],
                       trt = c(rep(0, length(mod[["strata"]][mod[["strata"]] == "trt=0"])), 
                               rep(1, length(mod[["strata"]][mod[["strata"]] == "trt=1"])))
  ) %>%
    filter(t < 43 + 0.5) %>%  # Deal with the jittering of outcomes
    group_by(trt) %>% 
    summarize(risk = last(r),
              .groups = 'drop') %>% 
    pivot_wider(names_from = trt,
                values_from = risk,
                names_glue = "{.value}{trt}") %>% 
    select(risk0, risk1)
  
  return(summod)
}

           



# aj_estimator <- function(data_subset, outcome_ind, time_var) {
#   outcome_ind <- enquo(outcome_ind)
#   time_var <- enquo(time_var)
#   
#   # Jitter ties
#   data_jittered <- data_subset %>% 
#     group_by(!!time_var) %>% 
#     add_tally(name = "count") %>%  # Use a different name for the tally counts
#     ungroup() %>% 
#     mutate(
#       # Jitter event times
#       jitter = runif(nrow(data_subset), min = -.01, max = .01),
#       time = ifelse(count > 1, !!sym(quo_name(time_var)) + jitter, !!sym(quo_name(time_var)))
#     )
#   
#   # Run the AJ model
#   aj <- survfit(Surv(time, factor(!!sym(quo_name(outcome_ind)))) ~ trt, data = data_jittered)
#   
#   mod <- summary(aj)
#   
#   summod <- data.frame(t = mod$time,
#                        r = mod$pstate[, 2], 
#                        se = mod$std.err[, 2],
#                        trt = c(rep(0, length(mod[["strata"]][mod[["strata"]] == "trt=0"])), 
#                                rep(1, length(mod[["strata"]][mod[["strata"]] == "trt=1"])))
#   ) %>%
#     filter(t < 43 + 0.5) %>%  # Deal with the jittering of outcomes
#     group_by(trt) %>% 
#     summarize(risk = last(r),
#               .groups = 'drop') %>% 
#     pivot_wider(names_from = trt,
#                 values_from = risk,
#                 names_glue = "{.value}{trt}") %>% 
#     select(risk0, risk1)
#   
#   return(summod)
# }




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


## Revised

calculate_aj_risks <- function(dataset, pnc_dist, sev_dist, outcome_ind, time_var) {
  outcome_ind <- enquo(outcome_ind)
  time_var <- enquo(time_var)
  
  dataset <- dataset %>% 
    mutate(pnc_wk = as.character(pnc_wk))
  
  # Jitter ties
  mar <- dataset %>% 
    group_by(!!time_var) %>% 
    add_tally(name = "count") %>% 
    ungroup() %>% 
    mutate(
      # Jitter event times
      jitter = runif(n(), min = -.01, max = .01), 
      time = ifelse(count > 1, !!time_var + jitter, !!time_var)
    )
  
  # Split the dataset by severity and pnc_wk
  split_data <- split(mar, list(mar$severity, mar$pnc_wk))
  
  # Apply the aj_estimator function to each subset
  results_list <- lapply(split_data, function(subset) aj_estimator(subset, outcome_ind = !!outcome_ind, 
                                                                   time_var = !!time_var))
  
  # Combine the results into a single data frame
  results <- bind_rows(results_list, .id = "group")
  
  # Extract severity and pnc_wk from the group column
  results <- results %>%
    separate(group, into = c("severity", "pnc_wk"), sep = "\\.")
  
  # Convert severity and pnc_wk to their original data types if necessary
  results <- results %>%
    mutate(severity = as.numeric(severity),
           pnc_wk = as.character(pnc_wk))
  
  # Left merge the severity distribution onto the risks
  results_sev_merge <- left_join(results, sev_dist, by = c("severity" = "severity")) %>% 
    rowwise() %>% 
    mutate(risk0_std = risk0 * prop,
           risk1_std = risk1 * prop) %>% 
    group_by(pnc_wk) %>% 
    summarize(risk0 = sum(risk0_std),
              risk1 = sum(risk1_std))
  
  # Calculate the overall risks
  risk0 <- sum(results_sev_merge$risk0 * pnc_dist)
  risk1 <- sum(results_sev_merge$risk1 * pnc_dist)
  overall_risks <- tibble(
    pnc_wk = 'overall',
    risk0 = risk0,
    risk1 = risk1
  )
  
  # Final risks dataset
  risks <- bind_rows(results_sev_merge, overall_risks) %>% 
    rowwise() %>% 
    mutate(rd = risk1 - risk0,
           rr = risk1 / risk0)
  
  return(risks)
}




# 
# # Define the calculate_aj_risks function
# calculate_aj_risks <- function(dataset, pnc_dist, sev_dist, outcome_ind, time_var) {
#   outcome_ind <- enquo(outcome_ind)
#   time_var <- enquo(time_var)
#   
#   # Jitter ties
#   mar <- dataset %>% 
#     group_by(!!time_var) %>% 
#     add_tally(name = "count") %>% 
#     ungroup() %>% 
#     mutate(
#       # Jitter event times
#       jitter = runif(nrow(dataset), min = -.01, max = .01), 
#       time = ifelse(count > 1, !!sym(quo_name(time_var)) + jitter, !!sym(quo_name(time_var)))
#     )
#   
#   # Split the dataset by severity and pnc_wk
#   split_data <- split(mar, list(mar$severity, mar$pnc_wk))
#   
#   # Apply the aj_mar function to each subset
#   results_list <- lapply(split_data, function(subset) aj_estimator(subset, !!sym(quo_name(outcome_ind)), 
#                                                                    !!sym(quo_name(time_var))))
#   
#   # Combine the results into a single data frame
#   results <- bind_rows(results_list, .id = "group")
#   
#   # Extract severity and pnc_wk from the group column
#   results <- results %>%
#     separate(group, into = c("severity", "pnc_wk"), sep = "\\.")
#   
#   # Convert severity and pnc_wk to their original data types if necessary
#   results <- results %>%
#     mutate(severity = as.numeric(severity),
#            pnc_wk = as.numeric(pnc_wk))
#   
#   # Left merge the severity distribution onto the risks
#   results_sev_merge <- left_join(results, sev_dist, by = c("severity" = "severity")) %>% 
#     rowwise() %>% 
#     mutate(risk0_std = risk0 * prop,
#            risk1_std = risk1 * prop) %>% 
#     group_by(pnc_wk) %>% 
#     summarize(risk0 = sum(risk0_std),
#               risk1 = sum(risk1_std))
#   
#   # Calculate the overall risks
#   risk0 <- sum(results_sev_merge$risk0 * pnc_dist)
#   risk1 <- sum(results_sev_merge$risk1 * pnc_dist)
#   overall_risks <- tibble(
#     pnc_wk = 'overall',
#     risk0 = risk0,
#     risk1 = risk1
#   )
#   
#   # Final risks dataset
#   risks <- rbind(results_sev_merge, overall_risks) %>% 
#     rowwise() %>% 
#     mutate(rd = risk1 - risk0,
#            rr = risk1 / risk0)
#   
#   return(risks)
# }
















#OLD CODE


# calculate_risks <- function(dataset, pnc_dist, sev_dist){
#   
#   # Calculate the risks using the potential outcomes
#   strat_risks <- dataset %>% 
#     group_by(trt, pnc_wk, severity)%>% 
#     summarize(
#       risk = sum(preeclampsia_pre_miss) / n(),
#       .groups = 'drop'
#     )
#   
#   # Merge the severity distribution among everyone at baseline on
#   # and multiply through in order to calculate the standardized risks within 
#   # each trt and pnc_wk strata.
#   merge_sev_prop <- left_join(strat_risks, sev_dist, by = c("severity" = "severity")) %>% 
#     rowwise() %>% 
#     mutate(std_risk = risk*prop) %>% 
#     group_by(trt, pnc_wk) %>% 
#     summarize(risk = sum(std_risk),
#               .groups = 'drop')
#   
#   strata_risks_wide <- merge_sev_prop %>% 
#     pivot_wider(names_from = trt, values_from = risk, names_prefix = "risk")
#   
#   # Calculate the risks standardized to the PNC dist at baseline
#   risk0 <- sum(strata_risks_wide$risk0 * pnc_dist)
#   risk1 <- sum(strata_risks_wide$risk1 * pnc_dist)
#   
#   risks_overall <- tibble(
#     pnc_wk = 'all',
#     risk0 = risk0,
#     risk1 = risk1
#   )
#   
#   risks <- rbind(strata_risks_wide, risks_overall) %>% 
#     rowwise() %>% 
#     mutate(rd = risk1 - risk0,
#            rr = risk1 / risk0) %>% 
#     ungroup()
#   
#   return(risks)
#   
# }




