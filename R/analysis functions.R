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
  
  # Convert data to a data.table (if it's not already)
  data <- as.data.table(data)
  
  # Perform data cleaning and manipulation using data.table syntax
  data[, `:=`(
    
    # Make indicator variables for observed delivery
    obs_delivery_mar_mnar = fifelse(ltfu_mar_mnar == 'not' & pregout_t_mar_mnar >= 18, 1, 0),
    
    # Make indicator for observed pregnancy outcome
    obs_outcome_mar_mnar = fifelse(ltfu_mar_mnar == 'not', 1, 0),
    
    # Outcome indicator - multivariate, important for AJ estimator
    final_pregout_mar_mnar = fcase(
      pregout_mar_mnar == 'unknown', 0,
      preeclampsia_mar_mnar == 1, 1,
      pregout_mar_mnar == 'fetaldeath', 2,
      pregout_mar_mnar == 'livebirth', 3
    ),
    
    # Time to event - Establish the time-to-event for the outcome indicator
    final_pregout_marmnar_tte = fifelse(
      pregout_mar_mnar == 'unknown',
      t_ltfu_mar_mnar - pnc_wk,
      pregout_t_mar_mnar - pnc_wk
    ),
    
    # Sensitivity analyses
    preeclampsia_marmnar_immediate_outc_all = fifelse(pregout_mar_mnar == 'unknown', 1, preeclampsia_mar_mnar),
    
    preeclampsia_marmnar_immediate_outc_trt = fcase(
      pregout_mar_mnar == 'unknown' & trt == 1, 1,
      pregout_mar_mnar == 'unknown' & trt == 0, 0,
      pregout_mar_mnar != 'unknown', preeclampsia_mar_mnar
    ),
    
    preeclampsia_marmnar_immediate_outc_untrt = fcase(
      pregout_mar_mnar == 'unknown' & trt == 0, 1,
      pregout_mar_mnar == 'unknown' & trt == 1, 0,
      pregout_mar_mnar != 'unknown', preeclampsia_mar_mnar
    ),
    
    preeclampsia_marmnar_no_outc = fifelse(pregout_mar_mnar == 'unknown', 0, preeclampsia_mar_mnar)
  )]
  
  
  # Add small constant to follow-up time for immediate censors
  data[, final_pregout_marmnar_tte := fifelse(final_pregout_marmnar_tte == 0, 0.0001, final_pregout_marmnar_tte)]
  
  return(data)
}







#########################################
# FUNCTION: conduct_analysis()
# PURPOSE: Conduct the analyses for the
# dataset.
#########################################

conduct_analysis <- function(data){
  
  # Calculate the distribution of severity using data.table syntax
  sev_dist <- data[, .(prop = .N / nrow(data)), by = severity]
  
  # Get the risks from the potential outcomes
  pot_out_risks <- data[, .(
    risk0 = sum(preeclampsia0) / .N,
    risk1 = sum(preeclampsia1) / .N,
    rd = (sum(preeclampsia1) / .N) - (sum(preeclampsia0) / .N),
    rr = (sum(preeclampsia1) / sum(preeclampsia0))
  )] %>% 
    as_tibble()
  
  # Get the risks among observed deliveries
  observed_deliveries <- data[obs_delivery_mar_mnar == 1]
  
  observed_deliveries_mar_mnar <- calculate_risks(observed_deliveries, sev_dist) %>% 
    as_tibble()
  
  # Get the risks among observed pregnancy outcomes
  observed_outcomes <- data[obs_outcome_mar_mnar == 1]
  
  observed_outcomes_mar_mnar <- calculate_risks(observed_outcomes, sev_dist) %>% 
    as_tibble()
  
  # Time-to-event analyses
    
    tte_mar_mnar <- calculate_aj_risks(data, 
                                       sev_dist)
    
  # Sensitivity analyses
    sens_anal_risks <- data[, .(
      ## Assume that all missing outcomes have an outcome
      risk_all_outc = sum(preeclampsia_marmnar_immediate_outc_all) / .N,
      ## Assume that treated missing have an outcome
      risk_trt_outc = sum(preeclampsia_marmnar_immediate_outc_trt) / .N,
      ## Assume that untreated missing have an outcome
      risk_untrt_out = sum(preeclampsia_marmnar_immediate_outc_untrt) / .N,
      ## Assume that all missing outcomes do not have an outcome
      risk_no_outc = sum(preeclampsia_marmnar_no_outc) / .N
    ), by = .(trt, severity)]
    
    ## Merge severity distribution
    merge_sens_anal <- merge(sens_anal_risks, sev_dist, by = "severity", all.x = TRUE)
    
    ## Multiply through
    sens_analyses <- merge_sens_anal[, `:=`(
      all_outc = risk_all_outc * prop,
      trt_outc = risk_trt_outc * prop,
      untrt_outc = risk_untrt_out * prop,
      no_outc = risk_no_outc * prop
    )] %>% 
      as_tibble() %>% 
      group_by(trt) %>% 
      summarize(
        all_outc = sum(all_outc),
        trt_outc = sum(trt_outc),
        untrt_outc = sum(untrt_outc),
        no_outc = sum(no_outc)
      )
    
    
    ## Assume that all missing outcomes have an outcome
    sens_anal_mar_mnar_all_outc <- sens_analyses %>% 
      select(trt, all_outc) %>% 
      pivot_wider(names_from = trt,
                  values_from = all_outc,
                  names_prefix = "risk") %>% 
      mutate(
        rd = risk1 - risk0,
        rr = risk1 / risk0
      )
    
    ## Assume that treated missing have an outcome
    sens_anal_mar_mnar_trt_outc <- sens_analyses %>% 
      select(trt, trt_outc) %>% 
      pivot_wider(names_from = trt,
                  values_from = trt_outc,
                  names_prefix = "risk") %>% 
      mutate(
        rd = risk1 - risk0,
        rr = risk1 / risk0
      )
      
    ## Assume that untreated missing have an outcome
    sens_anal_mar_mnar_untrt_outc <- sens_analyses %>% 
      select(trt, untrt_outc) %>% 
      pivot_wider(names_from = trt,
                  values_from = untrt_outc,
                  names_prefix = "risk") %>% 
      mutate(
        rd = risk1 - risk0,
        rr = risk1 / risk0
      )
      
    ## Assume that all missing outcomes do not have an outcome
    sens_anal_mar_mnar_no_outc <- sens_analyses %>% 
      select(trt, no_outc) %>% 
      pivot_wider(names_from = trt,
                  values_from = no_outc,
                  names_prefix = "risk") %>% 
      mutate(
        rd = risk1 - risk0,
        rr = risk1 / risk0
      )

  
  return(tibble(
    sev_dist = list(sev_dist),
    potential_risks = list(pot_out_risks),
    observed_deliveries = list(observed_deliveries_mar_mnar),
    observed_outcomes = list(observed_outcomes_mar_mnar),
    tte = list(tte_mar_mnar),
    sens_anal_all_outc = list(sens_anal_mar_mnar_all_outc),
    sens_anal_trt_outc = list(sens_anal_mar_mnar_trt_outc),
    sens_anal_untrt_out = list(sens_anal_mar_mnar_untrt_outc),
    sens_anal_no_outc = list(sens_anal_mar_mnar_no_outc)
  ))

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
calculate_risks <- function(dataset, sev_dist){
  
  # Calculate the risks within strata of treatment and severity using the static variable preeclampsia_pre_miss
  strat_risks <- dataset[, .(
    risk = sum(preeclampsia_pre_miss) / .N
  ), by = .(trt, severity)]
  
  # Merge the severity distribution and calculate the standardized risks
  merge_sev_prop <- merge(strat_risks, sev_dist, by = "severity", all.x = TRUE)
  
  merge_sev_prop[, std_risk := risk * prop]
  
  # Summarize risks within each treatment group
  final_risks <- merge_sev_prop[, .(
    risk = sum(std_risk)
  ), by = trt]
  
  # Reshape data from long to wide format
  final_risks_wide <- dcast(final_risks, . ~ trt, value.var = "risk")[, .(risk0 = `0`, risk1 = `1`)]
  
  # Calculate risk difference and risk ratio
  final_risks_wide[, `:=`(
    rd = risk1 - risk0,
    rr = risk1 / risk0
  )]
  
  # Return the final results
  return(final_risks_wide)
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

aj_estimator <- function(data_subset) {
  
  # Run the AJ model
  aj <- survfit(Surv(time, outcome) ~ trt, data = data_subset)
  
  # Extract the summary of the survival fit
  mod <- summary(aj)
  
  # Create a data.table from the summary output
  summod <- data.table(
    t = mod$time,
    r = mod$pstate[, 2],
    trt = c(rep(0, length(mod$strata[mod$strata == "trt=0"])), 
            rep(1, length(mod$strata[mod$strata == "trt=1"])))
  )
  
  # Group by treatment and summarize the last risk value
  summod <- summod[, .(risk = last(r)), by = trt]
  
  # Reshape data from long to wide format
  summod_wide <- dcast(summod, . ~ trt, value.var = "risk")
  
  # Rename the columns to risk0 and risk1
  summod_wide <- summod_wide[, .(risk0 = `0`, risk1 = `1`)]
  
  return(summod_wide)
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
  
  # Calculate the count for each group
  dataset[, count := .N, by = .(trt, final_pregout_marmnar_tte)]

  # Jitter event times and calculate the adjusted time and outcome
  dataset[, jitter := runif(.N, min = -.01, max = .01)]
  dataset[, time := ifelse(count > 1, final_pregout_marmnar_tte + jitter, final_pregout_marmnar_tte)]
  dataset[, outcome := as.factor(final_pregout_mar_mnar)]
  
  results <- dataset[, aj_estimator(.SD), by = severity] 
  
  results[, severity := as.numeric(severity)]
  
  # Left merge the severity distribution onto the risks
  merge_sev <- merge(results, sev_dist, by = "severity", all.x = TRUE)
  
  merge_sev[, `:=`(
    risk0_std = risk0 * prop,
    risk1_std = risk1 * prop
  )]
  
  combined <- merge_sev[, .(
    risk0 = sum(risk0_std),
    risk1 = sum(risk1_std),
    rd = sum(risk1_std) - sum(risk0_std),
    rr = sum(risk1_std) / sum(risk0_std)
  )]
  
  return(combined)
}

