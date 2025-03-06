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
# tar_data <- tar_read(cohort_data_0.1_1.3_0.49_.0.2_7_1_0.5_Parameters_Abortion1_Preeclampsia05_EMM.xlsx_1_10000)
# data <- tar_data %>%
#   as.data.table()
# 
# data <- tar_read(data_prep_0.226_0.8_0_.0.2_7_0.5_0.5_Parameters_Abortion05_Preeclampsia05_EMM.xlsx_1_10000)



#########################################
# FUNCTION: prep_data_for_analysis()
# PURPOSE: Prepare the data for analysis
# functions.

# Checked: 08.09.2024
#########################################

prep_data_for_analysis <- function(data, pnc_wk){
  
  # Perform data cleaning and manipulation using data.table syntax
  data[, `:=`(
    
    # Make indicator variables for observed delivery
    obs_delivery_mar_mnar = fifelse(ltfu_mar_mnar == 'not' & pregout_t_mar_mnar >= 18, 
                                    1, 
                                    0),
    
    # Make indicator for observed pregnancy outcome
    obs_outcome_mar_mnar = fifelse(ltfu_mar_mnar == 'not', 
                                   1, 
                                   0),
    
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
    
    ## Among those LTFU, assume all have an outcome
    preeclampsia_marmnar_immediate_outc_all = fifelse(pregout_mar_mnar == 'unknown', 1, preeclampsia_mar_mnar),
    
    ## Among those LTFU, assume only the treated have an outcome
    preeclampsia_marmnar_immediate_outc_trt = fcase(
      pregout_mar_mnar == 'unknown' & trt == 1, 1,
      pregout_mar_mnar == 'unknown' & trt == 0, 0,
      pregout_mar_mnar != 'unknown', preeclampsia_mar_mnar
    ),
    
    ## Among those LTFU, assume only the untreated have an outcome
    preeclampsia_marmnar_immediate_outc_untrt = fcase(
      pregout_mar_mnar == 'unknown' & trt == 0, 1,
      pregout_mar_mnar == 'unknown' & trt == 1, 0,
      pregout_mar_mnar != 'unknown', preeclampsia_mar_mnar
    ),
    
    ## Among those LTFU, assume that none have an outcome
    preeclampsia_marmnar_no_outc = fifelse(pregout_mar_mnar == 'unknown', 0, preeclampsia_mar_mnar)
  )]
  
  
  # Add small constant to follow-up time for immediate censors
  data[, final_pregout_marmnar_tte := fifelse(final_pregout_marmnar_tte == 0, 0.0001,
                                              final_pregout_marmnar_tte)]
  
  return(data)
}








#########################################
# FUNCTION: sev_dist()
# PURPOSE: Get the baseline severity 
# distribution.
#########################################

sev_dist <- function(data){
  
  # Calculate the distribution of severity and rurality using data.table syntax
  severity_dist <- data[, .(prop = .N / nrow(data)), by = .(severity,rural)]
  
  return(severity_dist)
  
}











#########################################
# FUNCTION: potential_risks()
# PURPOSE: Get the risks using the 
# potential outcomes.
#########################################

potential_risks <- function(data){
  
  pot_out_risks <- data[, .(
    risk0 = sum(preeclampsia0) / .N,
    risk1 = sum(preeclampsia1) / .N,
    rd = (sum(preeclampsia1) / .N) - (sum(preeclampsia0) / .N),
    rr = (sum(preeclampsia1) / .N) / (sum(preeclampsia0) / .N),
    # Added OR -- Think that they should converge in these cases
    or = ((sum(preeclampsia1) / .N) / (1 - (sum(preeclampsia1) / .N))) / 
      (((sum(preeclampsia0) / .N) / (1 - (sum(preeclampsia0) / .N))))
  )] %>% 
    as_tibble()
  
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
calculate_risks <- function(dataset, severity_dist){
  
  # dataset <- tar_data[obs_outcome_mar_mnar == 1]
  
  # Calculate the risks within strata of treatment, severity, and rurality using the static variable preeclampsia_pre_miss
  strat_risks <- dataset[, .(
    risk = sum(preeclampsia_pre_miss == 1) / .N
  ), by = .(trt, severity, rural)]
  
  # Merge the severity distribution and calculate the standardized risks
  merge_sev_prop <- merge(strat_risks, severity_dist, by = c("severity", "rural"), all.x = TRUE)
  
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
    rr = risk1 / risk0,
    # Added OR to see if they converge as expected.
    or = (risk1 / (1-risk1)) / (risk0 / (1-risk0))
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
  # aj <- survfit(Surv(time, outcome) ~ trt, data = data_subset)
  aj <- cuminc(ftime = data_subset$time, fstatus = data_subset$outcome, group = data_subset$trt)
  
  # Extract the summary of the survival fit
  trt0_risk <- last(aj$`0 1`$est)
  trt1_risk <- last(aj$`1 1`$est)
  
  # Calculate risk difference and risk ratio directly
  risks <- data.table(
    risk0 = trt0_risk,
    risk1 = trt1_risk,
    rd = trt1_risk - trt0_risk,
    rr = trt1_risk / trt0_risk
  )
  
  return(risks[, .(risk0, risk1, rd, rr)])
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

# First, calculate the risks within each severity strata using an AJ estimator.

calculate_aj_risks_sev <- function(dataset, severity_dist) { 
  
  # dataset is a subset of the original dataset, subset by severity value
  
  # Calculate the count for each group
  dataset[, count := .N, by = .(trt, final_pregout_marmnar_tte)]
  
  # Jitter event times and calculate the adjusted time and outcome
  dataset[, jitter := runif(.N, min = -.01, max = .01)]
  dataset[, time := ifelse(count > 1, final_pregout_marmnar_tte + jitter, final_pregout_marmnar_tte)]
  dataset[, outcome := as.factor(final_pregout_mar_mnar)]
  
  results <- dataset[, aj_estimator(.SD), by = .(severity, rural)]
  
  results[, severity := as.numeric(severity)]
  results[, rural := as.numeric(rural)]
  
  # Left merge the severity distribution onto the risks
  merge_sev <- merge(results, severity_dist, by = c("severity", "rural"), all.x = TRUE)
  
  merge_sev[, `:=`(
    risk0_std = risk0 * prop,
    risk1_std = risk1 * prop
  )]
  
  return(merge_sev)
  
}


## Then, combine the standardized risks.
## Originally, these were together but separated to ensure all computation 
## could take place within one server session.

calculate_aj_risks <- function(high_sev, middle_sev, low_sev){
  
  hold <- rbind(high_sev, middle_sev, low_sev)
  
  combined <- hold[, .(
    risk0 = sum(risk0_std),
    risk1 = sum(risk1_std),
    rd = sum(risk1_std) - sum(risk0_std),
    rr = sum(risk1_std) / sum(risk0_std)
  )]
  
  return(combined)
  
}










##############################################
# FUNCTION: sensitivity_analysis_risks()
# PURPOSE: The purpose of this function is to 
# calculate risks and their contrasts for the 
# sensitivity analyses making different
# assumptions about pregnancies with missing
# outcomes.
##############################################

sensitivity_analysis_risks <- function(data, sev_dist){
  
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
  ), by = .(trt, severity, rural)]
  
  ## Merge severity distribution
  merge_sens_anal <- merge(sens_anal_risks, sev_dist, by = c("severity", "rural"), all.x = TRUE)
  
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
  
  return(sens_analyses)
  
}



















