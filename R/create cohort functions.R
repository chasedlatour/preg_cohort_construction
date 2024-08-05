#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: This program contains all the R functions
# for identifying the study cohorts.
#####################################################

# Testing:
# data <- tar_read(generated_data_0.8_0.8_Parameters_Abortion08_Preeclampsia08.xlsx_1_20000)
# Issue with createion of: cohort_data_0_0.7_1_.0.2_0.8_0.8_Parameters_Abortion08_Preeclampsia08.xlsx_1_20000

#####################################################
# FUNCTION: generate_cohort()
# PURPOSE:
# Generate each analytic cohort for a set of
# missingness parameters.
#####################################################

generate_cohort <- function(data, marginal_p_miss_severity, beta12, 
                            marginal_p_miss_miscarriage, gamma1,
                            pnc_wk){
  
  # Get the expected value of severity in the missing cohort
  ## Distribution of severity at baseline : 1/3, 1/3, 1/3
  ## P(Trt|Severity) : 0.4, 0.5, 0.6
  ## Expected distribution of severity among untreated: 4/15 = low, 5/15 = med, 6/15 = high
  ## Manually calculated based upon the values above.
  expected_sev <- (4/15*0) + (5/15*1) + (6/15*2)
  
  ## Balancing intercept for model
  b0_sev <- -log((1/marginal_p_miss_severity)-1) - (beta12 * expected_sev)
  
  # Create p_sev_beta - for logistic regression to determine LTFU due to severity
  p_sev_beta = c(b0_sev, beta12, 2*beta12)
  
  ## Split the data by sim_id to maintain balance within the datasets
  split_data <- split(data, data$sim_id)
  
  ## Select the observed cohort
  hold <- mapply(
    create_cohort,
    split_data, 
    MoreArgs = list(p_sev_beta = p_sev_beta,
                    marginal_p_miss_miscarriage = marginal_p_miss_miscarriage,
                    gamma1 = gamma1, 
                    pnc_wk = pnc_wk),
    SIMPLIFY = FALSE
  ) %>% bind_rows()
  
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

create_cohort <- function(dataset, p_sev_beta, 
                          marginal_p_miss_miscarriage, gamma1, pnc_wk){
  
  # Get the probabilities associated with the betas for each
  # severity level - calculated via logistic regression
  p_missing_sev <- c(logodds_to_p(p_sev_beta[1]),
                     logodds_to_p(p_sev_beta[1] + p_sev_beta[2]),
                     logodds_to_p(p_sev_beta[1] + p_sev_beta[3]))
  
  # Create all the potential outcomes and revise the prenatal care encounters
  data <- dataset %>% # all_outcomes %>%#  
    dplyr::mutate(
      
      ######################################
      ### Untreated Potential Outcomes
      
      untreated_po = pmap(list(preg_outcomes_untrt, preec_outcomes_untrt,
                               revised_preg, pnc_wk),
                          untreated_outcomes),
      
      ######################################
      ### Treated Potential Outcomes
      
      treated_po = pmap(list(preg_outcomes_trt, preec_outcomes_trt, 
                             revised_preg, pnc_wk),
                        treated_outcomes), 
      
      ######################################
      ### Prenatal Encounters
      ### Identify prenatal encounters after
      ### the first.
      
      pnc_enc_rev = pmap(list(pnc_enc, pnc_wk),
                         revise_pnc) #,
      
      
      ######################################
      ### Missing Indicators - COME BACK -- THIS WONT WORK ANYMORE
      # 
      # ltfu_sev = first_missing_sev(unlist(missing_by_sev),
      #                              pnc_wk)
      
    ) %>% 
    unnest_wider(c(untreated_po, treated_po)) 
  
  
  
  # Assign the treatment values among those where include ne 0
  data2 <- assign_trt(subset(data, include == 1)) %>% 
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
  
  ## Calculate the betas for the missingness due to miscarriage.
  ## Because this is a single simulation and not a full Monte Carlo simulation, this average 
  ## is based upon the observed data. In a multi-repetition simulation, we would calculate
  ## this in a large cohort and use that value.
  ## Because the seed is maintained and this is only calculated among the non-initiators, this 
  ## should be the same for all
  expected_ga_miscarriages <- as.double(subset(data2, trt == 0 & pregout_t_pre_miss < 18) %>% 
                                          ungroup() %>% 
                                          summarize(avg = mean(pregout_t_pre_miss)))
  gamma0 <- -log((1/marginal_p_miss_miscarriage)-1) - (gamma1 * expected_ga_miscarriages)
  
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  ## This is the vector of beta values
  p_miss_outcome <- c(gamma0, 
                      gamma1)
  
  # Finally, incorporate/create the missing outcomes
  data2b <- data2 %>% 
    # Now determine if missing outcome based upon outcome and gw of outcome
    mutate(
      
      ######################################
      ## LTFU by Severity
      
      # Generate missingness probabilities for each person
      p_missing_by_sev = p_missing_sev[severity+1],
      
      # Determine if someone was missing due to severity
      missing_by_sev = rbinom(n = length(p_missing_by_sev), size = 1, prob = p_missing_by_sev),
      
      # Determine when the person was LTFU (e.g., decided not to come back to the data)
      ltfu_sev = runif(n = length(missing_by_sev), min = pnc_wk, max = pregout_t_pre_miss),
      
      ## Make missing if not actually LTFU due to severity
      ltfu_sev = ifelse(missing_by_sev == 1,
                        ltfu_sev,
                        NA),
      
      ######################################
      ## LTFU by Severity
      
      # Calculate the probability that missing due to outcome
      p_out_ltfu = ifelse(pregout_pre_miss == 'fetaldeath' & 
                            pregout_t_pre_miss < 20,
                          logodds_to_p(p_miss_outcome[1] + 
                                         (pregout_t_pre_miss * 
                                            p_miss_outcome[2])),
                          0),
      
      # Determine if they were LTFU due to the true pregnancy outcome
      ltfu_out = rbinom(n = length(p_out_ltfu), size = 1, prob = p_out_ltfu)
      
    )
  
  # Create a MAR and MNAR scenario
  data3 <- data2b %>% 
    mutate(
      
      ## Create an indicator variable for if the person was LTFU and how
      # -- This would occur if either:
      # -- - ltfu_sev is not missing and less than or equal to pregout_t_pre_miss, OR
      # -- - ltfu_out is equal to 1
      # Create an indicator just for MAR only and MAR+MNAR.
      ltfu_mar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
                        "sev",
                        "not"),
      ltfu_mar_mnar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
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
    mutate(preeclampsia_mar = ifelse(ltfu_mar != 'not',
                                     0,
                                     preeclampsia_pre_miss),
           pregout_mar = ifelse(ltfu_mar != 'not',
                                'unknown',
                                pregout_pre_miss),
           pregout_t_mar = ifelse(ltfu_mar != 'not',
                                  t_ltfu_mar,
                                  pregout_t_pre_miss),
           preeclampsia_mar_mnar = ifelse(ltfu_mar_mnar != 'not',
                                          0,
                                          preeclampsia_pre_miss),
           pregout_mar_mnar = ifelse(ltfu_mar_mnar != 'not',
                                     'unknown',
                                     pregout_pre_miss),
           pregout_t_mar_mnar = ifelse(ltfu_mar_mnar != 'not',
                                       t_ltfu_mar_mnar,
                                       pregout_t_pre_miss))
 
  return(data3)
  
}






### FUNCTION: logodds_to_p()
# This function converts the log-odds into a
# probability. Important for deriving probabilities
# from the logit function.
# INPUT:
# - logodds = log-odds value that needs to be 
#   converted to a probability

logodds_to_p <- function(logodds){
  
  p <- exp(logodds)/(1 + exp(logodds))
  
  return(p)
  
}







### FUNCTION: find_first_not_contpreg()
# Find the first element of a preg outcomes list that 
# is not continuing pregnancy.
find_first_not_contpreg <- function(vec) {
  idx <- which(vec != "contpreg_next")
  if (length(idx) == 0) {
    return(NA)
  } else {
    return(idx[1])
  }
}




### FUNCTION: find_first_preeclampsia()
# Find the first element of a preg outcomes list that 
# is preeclampsia
find_first_preeclampsia <- function(vec) {
  idx <- which(vec == 1)
  if (length(idx) == 0) {
    return(NA)
  } else {
    return(idx[1])
  }
}




#### FUNCTION: untreated_outcomes()
# This function the untreated
# potential outcomes
untreated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes, pnc_wk){
  
  ## Find the first pregnancy outcome
  # Gestational week that the outcome occurred (i.e., 1 + selected week)
  # Don't need to add 1 because Excel indexes from gest week 0 and R from 1
  first_preg_t = find_first_not_contpreg(unlist(preg_outcomes)) #+ 1
  
  ## Determine if the person developed preeclampsia if untreated
  first_preec_t = find_first_preeclampsia(unlist(preeclampsia_outcomes)) # + 1
  
  ## A person only actually develops preeclampsia if it occurs prior
  ## to or the same week as the untreated pregnancy outcome.
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







### FUNCTION: find_first_not_contpreg_wk()
# Find the first element of a preg outcomes list that 
# is not continuing pregnancy.
# However, only return those after the treated week
find_first_not_contpreg_wk <- function(vec, wk) {
  idx <- which(vec != "contpreg_next")
  
  # First index of the first value greater than wk
  # This is indexed at a prenatal care visit at gestational week wk
  # which is indexed by R as wk+1.
  # We want to identify the first outcome that occurs after the first
  # prenatal care encounter at gestational week wk or R index wk+1.
  # Outcomes are determined the week prior to when they occur. Thus,
  # we want the first outcome where it is determined at or after the 
  # prenatal care visit (week+1). Thus, it needs to be greater than wk
  index_greater_than_wk <- which(idx > wk)[1] 
  
  if (length(index_greater_than_wk) == 0) {
    return(NA)
  } else {
    return(idx[index_greater_than_wk])
  }
}



#### FUNCTION: treated_outcomes()
# This function identifies the treated
# potential outcomes
treated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes, wk){
  
  ## Determine if the person developed preeclampsia if treated
  first_preec_t = find_first_preeclampsia(unlist(preeclampsia_outcomes))
  # Note: Preeclampsia cannot develop until after the prenatal timings of interest
  
  ## Find the first pregnancy outcome that occurred after the indexing
  ## prenatal encounter -- Only way that we will see these outcomes
  first_preg_t = find_first_not_contpreg_wk(unlist(preg_outcomes), wk) 
  
  ## A person only actually develops preeclampsia if it occurs prior
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



#### FUNCTION: revise_pnc()
# Revise the occurrence of prenatal care
# encounters depending upon their gestational
# week of indexing (i.e., their first prenatal
# care encounter)
# INPUTs:
# - pnc_encounters - vector of the prenatal care
#   encounter indicators
# - wk - week of the indexing prenatal care encounter

revise_pnc <- function(pnc_encounters, wk){
  
  n <- length(unlist(pnc_encounters))
  
  # No prenatal encounters prior to indexing prenatal encounter
  pre_week <- rep(0, wk) # Changed wk-1 to wk because starts at 0
  
  
  # Prenatal encounter on indexing prenatal encounter week
  #week <- 1
  
  # Prenatal encounters after the indexing prenatal encounter
  wk_plus <- unlist(pnc_encounters)[(wk+1):n] 
  
  # Return the final list of prenatal encounters
  list(
    c(pre_week,
      #week, 
      wk_plus)
    )
  
}






##### FUNCTION assign_trt()
# This function assigns a treatment value
# to each person in a cohort dependent upon their 
# severity values.
# This is done within each cohort to only assign
# treatment to those that make it into the 
# study, as would be done in a real study.
# This also maintains statistical balance.

assign_trt <- function(dataset){
  
  # Pull out the p_trt vector
  p_trt <- dataset$p_trt
  # p_trt <- data$p_trt
  
  # Generate treatments
  trt <- rbinom(n = length(p_trt),
                size = 1,
                prob = p_trt)
  
  dataset$trt <- trt
  
  return(dataset)
  
}






### FUNCTION: pnc_miss()
# This function takes a vector of PNC encounters 
# and the ltfu indicator and determines the timing
# of the most recent PNC.
# Someone will be LTFU at their last PNC encounter
# PRIOR to that indicator.
pnc_miss <- function(pnc_enc, ltfu_ind, t_preg_out, t_sev){
  
  # Unlist for easy usage
  pnc_enc <- unlist(pnc_enc)
  
  # Identify the gw of PNC encounters
  # Need to subtract 1 because LTFU due to sev and 
  # LTFU due to outcome are occurring the week after
  # they are determined.
  # However, we want the last PNC BEFORE that visit.
  # Thus, need to subtract 1.
  which1 <- which(pnc_enc == 1) - 1
  
  if(ltfu_ind == "out"){
    subset <- which1[which1 < t_preg_out]
    last_pnc <- last(subset)
  } else if (ltfu_ind == "sev"){
    subset <- which1[which1 < t_sev]
    last_pnc <- last(subset)
  } else {
    last_pnc <- NA
  }
  
  return(last_pnc)
  
}



### FUNCTION: logodds_to_p()
# This function converts the log-odds into a
# probability. Important for deriving probabilities
# from the logit function.
logodds_to_p <- function(logodds){
  
  p <- exp(logodds)/(1 + exp(logodds))
  
  return(p)
  
}


