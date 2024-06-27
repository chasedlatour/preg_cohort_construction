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

## Upload the dataset for testing
#### Missing: Beta1 = 0.01, Gamma0 = 0.1, Gamma1 = 0.001
data <- readRDS("ab08preec08_beta001_gamma01_001.rds") %>% 
  mutate(scenario = "test1")

library(survival)



# Get the distribution by timing of first PNC encounter
pnc_dist <- (data %>% 
  group_by(pnc_wk) %>% 
  summarize(prop = n() / nrow(data)))$prop

sev_dist <- data %>% 
  group_by(severity) %>% 
  summarize(prop = n() / nrow(data),
            .groups = 'drop')

# First, need to do some data cleaning
data1 <- data %>% 
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
    
    ## Time to event - ADDED 1 TO INCLUDE IMMEDIATE CENSORS - CHECK THIS.
    final_pregout_mar_tte = ifelse(pregout_mar == 'unknown',
                                   t_ltfu_mar + 1,
                                   pregout_t_mar + 1),
    
    # MAR+MNAR
    ## Outcome indicator
    final_pregout_mar_mnar = case_when(pregout_mar_mnar == 'unknown' ~ 0,
                                       preeclampsia_mar == 1 ~ 1,
                                       pregout_mar == 'fetaldeath' ~ 2,
                                       pregout_mar == 'livebirth' ~ 3),
    
    ## Time to event - ADDED 1 TO INCLUDE IMMEDIATE CENSORS.
    final_pregout_marmnar_tte = ifelse(pregout_mar == 'unknown',
                                       t_ltfu_mar_mnar + 1,
                                       pregout_t_mar_mnar + 1)
    
  )

# Calculate the risks among cohorts with no right-censoring.

pot_risks <- calculate_pot_risks (data1, pnc_dist)

risks_obsdel_mar <- calculate_risks(subset(data1, obs_delivery_mar == 1), pnc_dist, sev_dist)
risks_obsdel_mar_mnar <- calculate_risks(subset(data1, obs_delivery_mar_mnar == 1), pnc_dist, sev_dist)

risks_obsout_mar <- calculate_risks(subset(data1, obs_outcome_mar == 1), pnc_dist, sev_dist)
risks_obsout_mar_mnar <- calculate_risks(subset(data1, obs_outcome_mar_mnar == 1), pnc_dist, sev_dist)


  


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

calculate_risks <- function(dataset, pnc_dist, sev_dist){
  
  # Calculate the risks using the potential outcomes
  strat_risks <- dataset %>% 
    group_by(trt, pnc_wk, severity)%>% 
    summarize(
      risk = sum(preeclampsia_pre_miss) / n(),
      .groups = 'drop'
    )
  
  # Merge the severity distribution among everyone at baseline on
  # and multiply through in order to calculate the standardized risks within 
  # each trt and pnc_wk strata.
  merge_sev_prop <- left_join(strat_risks, sev_dist, by = c("severity" = "severity")) %>% 
    rowwise() %>% 
    mutate(std_risk = risk*prop) %>% 
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
  
  risks <- rbind(strata_risks_wide, risks_overall) %>% 
    rowwise() %>% 
    mutate(rd = risk1 - risk0,
           rr = risk1 / risk0) %>% 
    ungroup()
  
  return(risks)
  
}







##############################################
# FUNCTION: calculate_aj_risks()
# PURPOSE: The purpose of this function to 
# calculate the risk of preeclampsia in a 
# closed cohort without right censoring. Only
# one obs per person.

# This implements non-parametric direct
# standardization by the distribution of severity
# at baseline to estimate the population
# average treatment effect (ATE).
##############################################


calculate_aj_risks() <- function(dataset, pnc_dist, sev_dist){
  
  
  #COME BACK -- NEED TO DO THIS WITHIN STRATA OF PNC VISIT TIMING
  
  ### MAR Outcomes
  
  #jitter ties
  dataset2 <- dataset %>% 
    group_by(final_pregout_mar_tte) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(
      # Jitter event times
      jitter = runif(nrow(dataset), min = -.01, max = .01), # N in the samp
      time = ifelse(n>1, final_pregout_mar_tte + jitter, final_pregout_mar_tte)
    )
  
  # Run the AJ model
  aj <- survfit(Surv(time, factor(final_pregout_mar)) ~ trt, data = dataset2)
  
  mod <- summary(aj)
  
  summod <- data.frame(t = mod$time,
                       r = mod$pstate[,2], 
                       se = mod$std.err[,2],
                       trt = c(rep(0, length(mod[["strata"]][mod[["strata"]] == "trt=0"])), 
                               rep(1, length(mod[["strata"]][mod[["strata"]] == "trt=1"])))
  ) %>%
    filter(t < 43 + 0.5) %>%  # Deal with the jittering of outcomes
    group_by(trt) %>% 
    summarize(risk = last(r),
              #se = last(se),
              .groups = 'drop') %>% 
    pivot_wider(names_from = trt,
                values_from = risk,
                #values_from = c(r, se), 
                names_glue = "{.value}{trt}") %>% 
    rowwise() %>% 
    mutate(rr = risk1/risk0,
           #selnrr = sqrt((1/r_0)^2 * se_0^2 + (1/r_1)^2*se_1^2),
           #rr_lcl = exp(log(rr) - z*selnrr),
           #rr_ucl = exp(log(rr) + z*selnrr),
           rd = risk1 - risk0 #,
           #se_rd = sqrt(se_1^2 + se_0^2),
           #rd_lcl = rd - 1.96*se_rd,
           #rd_ucl = rd + 1.96*se_rd,
           #Estimator = "Aalen-Johanssen"
           ) %>%
    select(risk0, risk1, rr, #rr_lcl, rr_ucl,
           rd #, rd_lcl, rd_ucl
    )
    
  
  
  ### MNAR Outcomes
  
}



calculate_aj_risks() <- function(dataset, pnc_dist, sev_dist){
  
  
  #COME BACK -- NEED TO DO THIS WITHIN STRATA OF PNC VISIT TIMING
  
  ### MAR Outcomes
  
  #jitter ties
  mar <- dataset %>% 
    group_by(final_pregout_mar_tte) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(
      # Jitter event times
      jitter = runif(nrow(dataset), min = -.01, max = .01), # N in the samp
      time = ifelse(n>1, final_pregout_mar_tte + jitter, final_pregout_mar_tte)
    )
  
  # Split the dataset by severity and pnc_wk
  split_data <- split(mar, list(mar$severity, mar$pnc_wk))
  
  # Apply the aj_analysis function to each subset
  results_list <- lapply(split_data, aj_mar)
  
  # Combine the results into a single data frame
  results <- bind_rows(results_list, .id = "group")
  
  # Extract severity and pnc_wk from the group column
  results <- results %>%
    separate(group, into = c("severity", "pnc_wk"), sep = "\\.")
  
  # Convert severity and pnc_wk to their original data types if necessary
  results <- results %>%
    mutate(severity = as.numeric(severity),
           pnc_wk = as.numeric(pnc_wk))

  
  ### MNAR Outcomes
  
}







## AJ Function for MAR

aj_mar <- function(data_subset) {
  # Jitter ties
  data_jittered <- data_subset %>% 
    group_by(final_pregout_mar_tte) %>% 
    add_tally(name = "count") %>%  # Use a different name for the tally counts
    ungroup() %>% 
    mutate(
      # Jitter event times
      jitter = runif(nrow(data_subset), min = -.01, max = .01),
      time = ifelse(count > 1, final_pregout_mar_tte + jitter, final_pregout_mar_tte)
    )
  
  # Run the AJ model
  aj <- survfit(Surv(time, factor(final_pregout_mar)) ~ trt, data = data_jittered)
  
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
    rowwise() %>% 
    mutate(rr = risk1 / risk0,
           rd = risk1 - risk0
    ) %>%
    select(risk0, risk1, rr, rd)
  
  return(summod)
}


## AJ Estimator

aj_estimator <- function(dataset, outcome_var, outcome_var_t, t_val, rename){
  
  z <- qnorm(0.975)
  
  #jitter ties
  #set.seed(1234)
  # Not sure if should set seed for jittering outcomes
  dataset2 <- dataset %>% 
    group_by(get(outcome_var_t)) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(
      # Jitter event times
      jitter = runif(nrow(dataset), min = -.01, max = .01), # N in the samp
      time = ifelse(n>1, get(outcome_var_t) + jitter, get(outcome_var_t))
    )
  
  # Run the AJ model
  aj <- survfit(Surv(time, factor(get(outcome_var))) ~ trt, data = dataset2)
  
  mod <- summary(aj)
  
  summod <- data.frame(t = mod$time,
                       r = mod$pstate[,2], 
                       se = mod$std.err[,2],
                       trt = c(rep(0, length(mod[["strata"]][mod[["strata"]] == "trt=0"])), 
                               rep(1, length(mod[["strata"]][mod[["strata"]] == "trt=1"])))
  ) %>%
    filter(t < t_val+0.5) %>%  # Deal with the jittering of outcomes
    group_by(trt) %>% 
    summarize(r = last(r),
              se = last(se),
              .groups = 'drop') %>% 
    pivot_wider(names_from = trt,
                values_from = c(r, se), 
                names_glue = "{.value}_{trt}") %>% 
    rowwise() %>% 
    mutate(rr = r_1/r_0,
           selnrr = sqrt((1/r_0)^2 * se_0^2 + (1/r_1)^2*se_1^2),
           rr_lcl = exp(log(rr) - z*selnrr),
           rr_ucl = exp(log(rr) + z*selnrr),
           rd = r_1 - r_0,
           se_rd = sqrt(se_1^2 + se_0^2),
           rd_lcl = rd - 1.96*se_rd,
           rd_ucl = rd + 1.96*se_rd,
           Estimator = "Aalen-Johanssen") %>%
    select(r_0, r_1, rr, rr_lcl, rr_ucl,
           rd, rd_lcl, rd_ucl
    )  %>% 
    rename_all((~ paste(rename, .)))
  
  return(summod)
  
  
}




