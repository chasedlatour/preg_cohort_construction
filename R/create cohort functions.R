#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: This program contains all the R functions
# for identifying the study cohorts.
#####################################################

# Testing:
# tar_data <- tar_read(generated_data_2_0.5_Parameters_Abortion2_Preeclampsia05_EMM.xlsx_1_10000)
# dataset <- tar_data
# data <- tar_data
# marginal_p_miss_severity <- 0.02825
# beta12 <- 0.8
# marginal_p_miss_miscarriage <- 0 #0.1225
# gamma1 <- -0.2





#####################################################
# FUNCTION: generate_cohort()
# PURPOSE:
# Generate each analytic cohort for a set of
# missingness parameters.
#####################################################

generate_cohort <- function(data, marginal_p_miss_severity, beta12, 
                            marginal_p_miss_miscarriage, gamma1,
                            pnc_wk){
  
  hold = create_cohort(data, marginal_p_miss_severity, beta12,
                       marginal_p_miss_miscarriage, 
                       gamma1, pnc_wk)
  
  return(hold)
  
}











#####################################################
# FUNCTION: create_cohort()
# PURPOSE: This function takes the raw generated DGM
# data and imposes on it the missingness, assigns trt,
# and then creates analytic cohorts from the raw data.
#
# INPUTs:
# - dataset = raw generated DGM data
# - p_sev_beta = betas for the missingness logistic
#   regression based upon severity
# - gamma1 = beta for the missingness 
#   logistic regression based upon preg outcome
#####################################################




create_cohort <- function(dataset, marginal_p_miss_severity, beta12,
                          marginal_p_miss_miscarriage, 
                          gamma1, pnc_wk){
  
  # Create all the potential outcomes and revise the prenatal care encounters
  data = dataset %>% # all_outcomes %>%#  
    dplyr::mutate(
      
      ######################################
      ### Untreated Potential Outcomes
      
      untreated_po = pmap(list(preg_outcomes_untrt, 
                               preec_outcomes_untrt,
                               revised_preg, 
                               pnc_wk),
                          untreated_outcomes),
      
      ######################################
      ### Treated Potential Outcomes
      
      treated_po = pmap(list(preg_outcomes_trt, 
                             preec_outcomes_trt, 
                             revised_preg, 
                             pnc_wk),
                        treated_outcomes), 
      
      ######################################
      ### Prenatal Encounters
      ### Identify prenatal encounters after
      ### the first.
      
      pnc_enc_rev = pmap(list(pnc_enc, 
                              pnc_wk),
                         revise_pnc) 
      
    ) %>% 
    unnest_wider(c(untreated_po, treated_po)) %>% 
    filter(include == 1) %>% 
    # Assign treatment -- Wait until after filter to those who are included in the sample
    mutate(trt = rbinom(n = length(p_trt), size = 1, prob = p_trt)) %>% 
    # Create variables for the observed outcomes - ignore missingness at this stage
    mutate(
      preeclampsia_pre_miss = ifelse(trt == 0,
                                     preeclampsia0,
                                     preeclampsia1),
      pregout_pre_miss = ifelse(trt == 0,
                                final_preg0,
                                final_preg1),
      pregout_t_pre_miss = ifelse(trt == 0,
                                  final_preg0_t,
                                  final_preg1_t)
    ) 
  
  
  ## Calculate the betas for the missingness due to severity.
  ## Because this is a single simulation and not a full Monte Carlo simulation, this average 
  ## is based upon the observed data. In a multi-repetition simulation, we would calculate
  ## this in a large cohort and use that value.
  ## Because the seed is maintained and this is only calculated among the non-initiators, this 
  ## should be the same for all
  
  # Get the expected value of severity in the missing cohort
  expected_sev = as.double(subset(data, trt == 0) %>% 
                             ungroup() %>% 
                             summarize(avg = mean(severity)))

  # Calculate the weekly marginal severity for a total marginal probability of marginal_p_miss_severity
  ## Number of weeks between the end of follow-up and pnc_wk
  num_weeks <- 41 - pnc_wk
  ## Weekly probability that ltfu to maintain the marginal probability
  weekly_prob <- 1 - ((1-marginal_p_miss_severity)^(1/num_weeks))

  ## Calculate intercept using balancing intercept for model
  b0_sev = -log((1/weekly_prob)-1) - (beta12 * expected_sev)
  
  # Create p_sev_beta - for logistic regression to determine LTFU due to severity
  p_sev_beta = c(b0_sev, beta12, 2*beta12)
  
  # Get the probabilities associated with the betas for each
  # severity level - calculated via logistic regression
  p_missing_sev = c(logodds_to_p(p_sev_beta[1]),
                    logodds_to_p(p_sev_beta[1] + p_sev_beta[2]),
                    logodds_to_p(p_sev_beta[1] + p_sev_beta[3]))
  
  ## Calculate the betas for the missingness due to miscarriage.
  ## Same logic as missingness due to severity.
  expected_ga_miscarriages = as.double(subset(data, trt == 0 & pregout_t_pre_miss < 18) %>% 
                                         ungroup() %>% 
                                         summarize(avg = mean(pregout_t_pre_miss)))
  # Balancing intercept term
  gamma0 = -log((1/marginal_p_miss_miscarriage)-1) - (gamma1 * expected_ga_miscarriages)
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  ## This is the vector of beta values
  p_miss_outcome = c(gamma0, 
                     gamma1)
  
  # Finally, incorporate/create the missing outcomes
  data2b = data %>% 
    # Now determine if missing outcome based upon outcome and gw of outcome
    mutate(
      
      ######################################
      ## LTFU by Severity
      
      # Generate missingness probabilities for each person
      p_missing_by_sev = p_missing_sev[severity+1],
      
      # Create a weekly indicator variable for is someone was missing due to severity
      missing_by_sev = purrr::map(p_missing_by_sev, ~rbinom(n = 41, size = 1, prob = .x)),
      
      # # Determine if someone was missing due to severity
      # missing_by_sev = rbinom(n = length(p_missing_by_sev), size = 1, prob = p_missing_by_sev),
      
      # Determine all the weeks that someone was determined to be ltfu due to severity
      # Censoring determined the week before it is observed. This is accomodated by R by indexing from 1.
      wks_missing_by_sev = purrr::map(missing_by_sev, ~which(.x == 1)),
      
      # Identify the first ltfu due to severity indicator that occurs after PNC week
      ## If none, then this is NA
      ltfu_sev = purrr::map_dbl(wks_missing_by_sev,
                            ~ifelse(is.na(which(.x > pnc_wk)[1]),
                                    NA,
                                    .x[which(.x > pnc_wk)[1]])),
      
      # # Determine when the person was LTFU due to severity
      # ltfu_sev = runif(n = length(missing_by_sev), min = pnc_wk, max = pregout_t_pre_miss),
      # 
      # ## Make missing if not actually LTFU due to severity
      # ltfu_sev = ifelse(missing_by_sev == 1,
      #                   ltfu_sev,
      #                   NA),
      
      ######################################
      ## LTFU by Outcome
      
      # Calculate the probability that missing due to outcome
      p_out_ltfu = ifelse(pregout_pre_miss == 'fetaldeath' & 
                            pregout_t_pre_miss < 18, # Abortion definition
                          logodds_to_p(p_miss_outcome[1] + 
                                         (pregout_t_pre_miss * 
                                            p_miss_outcome[2])),
                          0),
      
      # Determine if they were LTFU due to the true pregnancy outcome
      ltfu_out = rbinom(n = length(p_out_ltfu), size = 1, prob = p_out_ltfu)
      
    ) %>% 
    # Create the observed outcomes, incorporating missing data
    mutate(
      
      ## Create an indicator variable for if the person was LTFU and how
      # -- This would occur if either:
      # -- - ltfu_sev is not missing and less than or equal to pregout_t_pre_miss, OR
      # -- - ltfu_out is equal to 1
      # Create an indicator just for MAR only and MAR+MNAR.
      ltfu_mar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
                        "sev",
                        "not"),
      ltfu_mar_mnar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss, # & ltfu_sev <= pregout_t_pre_miss
                             "sev",
                             ifelse(ltfu_out == 1,
                                    "out",
                                    "not")),
      
      # Now determine timing when LTFU
      t_ltfu_mar = pmap_dbl(list(pnc_enc_rev,
                                 ltfu_mar,
                                 41,
                                 ltfu_sev), # Just set as a constant-ignoring missing due to outcome
                            pnc_miss),
      
      t_ltfu_mar_mnar = pmap_dbl(list(pnc_enc_rev,
                                      ltfu_mar_mnar,
                                      pregout_t_pre_miss,
                                      ltfu_sev),
                                 pnc_miss)
      
    ) %>%
    # Finally, create their observed outcomes, incorporating the LTFU
    mutate(
      
      ## Only missing by severity
      preeclampsia_mar = ifelse(ltfu_mar == 'not',
                                preeclampsia_pre_miss,
                                0),
      pregout_mar = ifelse(ltfu_mar == 'not',
                           pregout_pre_miss,
                           'unknown'),
      pregout_t_mar = ifelse(ltfu_mar == 'not',
                             pregout_t_pre_miss,
                             t_ltfu_mar),
      
      ## Missing by severity and outcome
      preeclampsia_mar_mnar = ifelse(ltfu_mar_mnar == 'not',
                                     preeclampsia_pre_miss,
                                     0),
      pregout_mar_mnar = ifelse(ltfu_mar_mnar == 'not',
                                pregout_pre_miss,
                                'unknown'),
      pregout_t_mar_mnar = ifelse(ltfu_mar_mnar != 'not',
                                  t_ltfu_mar_mnar,
                                  pregout_t_pre_miss)
    )
  
  return(data2b)
  
}













#####################################################
# FUNCTION: find_first_not_contpreg()
# Find the first element of a preg outcomes list that 
# is not continuing pregnancy.
#####################################################

find_first_not_contpreg <- function(vec) {
  idx = which(vec != "contpreg_next")
  if (length(idx) == 0) { # This should never occur, a check
    return(NA)
  } else {
    return(idx[1])
  }
}





#####################################################
# FUNCTION: find_first_preeclampsia()
# Find the first element of a preg outcomes list that 
# is preeclampsia
#####################################################

find_first_preeclampsia <- function(vec) {
  idx = which(vec == 1)
  if (length(idx) == 0) {
    return(NA) # This might occur if they don't develop preeclampsia ever
  } else {
    return(idx[1])
  }
}






#####################################################
#### FUNCTION: untreated_outcomes()
# This function identified the untreated
# potential outcomes
# INPUTS:
# - preg_outcomes - the list of pregnancy outcomes
#   for a person
# - preeclampsia_outcomes - the list of preeclampsia
#   outcomes for a person
# - revised_outcomes - the list of pregnancy outcomes
#   if the person were to develop preeclampsia
# - pnc_wk - the week of the person's initial 
#   prenatal care encounter
#####################################################

untreated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes, pnc_wk){
  
  ## Find the first pregnancy outcome that is not continuing pregnancy
  # Gestational week that the outcome occurred (i.e., 1 + selected week)
  # Don't need to add 1 because Excel indexes from gest week 0 and R from 1
  first_preg_t = find_first_not_contpreg(unlist(preg_outcomes))
  
  ## Determine if the person developed preeclampsia if untreated
  first_preec_t = find_first_preeclampsia(unlist(preeclampsia_outcomes))
  
  ## A person only actually develops preeclampsia if it occurs prior
  ## to or the same week as the untreated pregnancy outcome.
  preeclampsia = ifelse(is.na(first_preec_t),
                        0,
                        as.numeric(first_preec_t <= first_preg_t))
  
  ## Now, get the revised outcome based upon their preeclampsia value
  
  # Determine the gestational week that the outcome occurs
  final_preg_t = ifelse(preeclampsia == 0, # If no preeclampsia
                        first_preg_t, # Then potential outcome
                        first_preec_t) # Else preeclampsia timing
  
  # Get the actual outcome value
  final_preg = ifelse(preeclampsia == 0,
                      unlist(preg_outcomes)[final_preg_t], 
                      unlist(revised_outcomes)[final_preg_t]) 
  
  final_preg = gsub("_next", "", final_preg)
  
  
  # Determine if the person would be indexed into the cohort at 4, 7, 
  # and 16 weeks of gestation -- A person had to NOT have had a fetal
  # death prior to that gestational week to be included
  # Make an indicator variable for if they should be included
  if (final_preg_t > pnc_wk){
    include = 1
  } else {
    include = 0
  }
  
  
  # Return a list of the desired values
  list(
    preeclampsia0 = preeclampsia,
    final_preg0_t = final_preg_t,
    final_preg0 = final_preg,
    include = include 
  )
  
}






#######################################################
### FUNCTION: find_first_not_contpreg_wk()
# Find the first element of a preg outcomes list that 
# is not continuing pregnancy.
# However, only return those after the treated week
# INPUTS:
# - vec: The vector of pregnancy outcomes
# - wk: Week of the initial prenatal encounter
#######################################################

find_first_not_contpreg_wk <- function(vec, wk) {
  
  idx = which(vec != "contpreg_next")
  
  # First index of the not continuing preg weeks that is greater than wk
  index_greater_than_wk = which(idx > wk)[1]
  
  if (length(index_greater_than_wk) == 0) {
    return(NA) # This should not occur for pregnancy outcomes, a check.
  } else {
    return(idx[index_greater_than_wk])
  }
}






#######################################################
# FUNCTION: treated_outcomes()
# This function identifies the treated
# potential outcomes
#
# INPUTS:
# - preg_outcomes: vector with the pregnancy outcomes
#   at each gestational week
# - preeclampsia_outcomes: vector with the preeclampsia
#   outcomes at each gestational week
# - revised_outcomes: vector with the pregnancy outcomes
#   post-preeclampsia at each gestational week
# - wk: week of the initial prenatal care visit
#######################################################

treated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes, wk){
  
  ## Determine if the person developed preeclampsia if treated
  first_preec_t = find_first_preeclampsia(unlist(preeclampsia_outcomes))
  
  ## Find the first pregnancy outcome that occurred after the indexing
  ## prenatal encounter -- Only see these outcomes after someone received
  ## trt at first PNC encounter
  first_preg_t = find_first_not_contpreg_wk(unlist(preg_outcomes), wk) 
  
  ## Create preeclampsia indicator: Must occur prior
  ## to or the same week as the treated pregnancy outcome.
  preeclampsia = ifelse(is.na(first_preec_t),
                        0,
                        as.numeric(first_preec_t <= first_preg_t))
  
  ## Now, get the revised outcome based upon their preeclampsia value
  # Gestational week that the outcome occurs
  final_preg_t = ifelse(preeclampsia == 0, # If no preeclampsia
                        first_preg_t, # Then potential outcome
                        first_preec_t) # Else preeclampsia timing
  
  # Get the actual outcome value
  final_preg = ifelse(preeclampsia == 0,
                      unlist(preg_outcomes)[final_preg_t], 
                      unlist(revised_outcomes)[final_preg_t])
  
  final_preg = gsub("_next", "", final_preg)
  
  
  # Return a list of the desired values
  list(
    preeclampsia1 = preeclampsia,
    final_preg1_t = final_preg_t,
    final_preg1 = final_preg
  )
  
}









###################################################
# FUNCTION: revise_pnc()
# Revise prenatal care encounters so that they 
# only occur on and after the first PNC encounter
#
# INPUTs:
# - pnc_encounters - vector of the prenatal care
#   encounter indicators
# - wk - week of the first prenatal care encounter
###################################################

revise_pnc <- function(pnc_encounters, wk){
  
  n = length(unlist(pnc_encounters))
  
  # No prenatal encounters prior to indexing prenatal encounter
  pre_week = rep(0, wk) # wk instead of wk-1 because pnc encs start at 0
  
  # Prenatal encounter on indexing prenatal encounter week
  week = 1
  
  # Prenatal encounters after the indexing prenatal encounter
  wk_plus = unlist(pnc_encounters)[(wk+2):n] 
  
  # Return the final list of prenatal encounters
  list(
    c(pre_week,
      week, 
      wk_plus)
  )
  
}









###################################################
### FUNCTION: pnc_miss()
# This function takes a vector of PNC encounters 
# and the ltfu indicator and determines the timing
# of the most recent PNC.
# Someone will be LTFU at their last PNC encounter
# PRIOR to that indicator.
###################################################


pnc_miss <- function(pnc_enc, ltfu_ind, t_preg_out, t_sev){
  
  # Unlist for easy usage
  pnc_enc = unlist(pnc_enc)
  
  # Identify the gw of PNC encounters
  # Subtract 1 because pnc enc indexed from week 0
  which1 = which(pnc_enc == 1) - 1
  
  if(ltfu_ind == "out"){
    subset = which1[which1 < t_preg_out]
    last_pnc = last(subset)
  } else if (ltfu_ind == "sev"){
    subset = which1[which1 < t_sev]
    last_pnc = last(subset)
  } else {
    last_pnc = NA
  }
  
  return(last_pnc)
  
}







#####################################################
### FUNCTION: logodds_to_p()
# This function converts the log-odds into a
# probability. Important for deriving probabilities
# from the logit function.
#####################################################

logodds_to_p <- function(logodds){
  
  p <- exp(logodds)/(1 + exp(logodds))

  return(p)
  
}
