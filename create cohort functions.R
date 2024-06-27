#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: This program contains all the R functions
# for identifying the study cohorts.
#####################################################




#####################################################
# FUNCTION: generate_cohort()
# PURPOSE:
# Generate each analytic cohort for a set of
# missingness parameters.
#####################################################

generate_cohort <- function(data, beta12, gamma0, gamma1, save_name_cohort){
  
  # Set values that are going to hold for ALL DGMs and 
  # analytic samples
  p_sev <- 1 - 0.8^(1/40) # Weekly prob (over 40 weeks) of being LTFU at lowest severity
  b0_sev <- log(p_sev / (1-p_sev)) # beta 0 in the logistic regression for LTFU due to severity is NOT varied across scenarios
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  p_miss_outcome <- c(log(gamma0/(1-gamma0)), 
                      log(gamma1/(1-gamma1))
  )
  
  # Create p_sev_beta - for logistic regression to determine LTFU due to severity
  p_sev_beta = c(b0_sev, log(beta12 / (1-beta12)), log((2*beta12 / (1 - (2*beta12)))))
  
  ## Select the observed cohort
  all_outcomes2 <- create_cohort(data, p_sev_beta, p_miss_outcome) %>% 
    mutate(diff_beta1_beta2 = beta12,
           gamma_0 = gamma0,
           gamma_1 = gamma1)
  
  # Save the RDS file
  saveRDS(all_outcomes2, save_name_cohort)
  
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
# - p_miss_outcome = betas for the missingness 
#   logistic regression based upon preg outcome
#####################################################

### CHASE -- RUN SOME SUMMARIES OF EACH OF THESE DATASETS TO MAKE SURE THEY LOOK AS EXPECTED
create_cohort <- function(dataset, p_sev_beta, p_miss_outcome){
  
  # Get the probabilities associated with the betas for each
  # severity level - calculated via logistic regression
  p_missing_sev <- c(logodds_to_p(p_sev_beta[1]),
                     logodds_to_p(p_sev_beta[1] + p_sev_beta[2]),
                     logodds_to_p(p_sev_beta[1] + p_sev_beta[3]))
  
  data <- dataset %>% #   all_outcomes %>%#
    dplyr::group_by(sim_id, id) %>% 
    dplyr::mutate(
      
      ######################################
      ## LTFU by Severity
      
      # Generate missingness probabilities for each person
      p_missing_by_sev = p_missing_sev[severity+1],
      
      # Select missingness values for each person's gestational week
      # manually input as 41
      missing_by_sev = list(rbinom(n=41, size=1, prob=p_missing_by_sev)),
      
      ######################################
      ### Untreated Potential Outcomes
      
      untreated_po = list(untreated_outcomes(unlist(preg_outcomes_untrt), 
                                             unlist(preec_outcomes_untrt), 
                                             unlist(revised_preg),
                                             pnc_wk)),
      
      ######################################
      ### Treated Potential Outcomes
      
      treated_po = list(treated_outcomes(unlist(preg_outcomes_trt),
                                         unlist(preec_outcomes_trt),
                                         unlist(revised_preg),
                                         pnc_wk)), 
      
      ######################################
      ### Prenatal Encounters
      ### Identify prenatal encounters after
      ### the first.
      
      pnc_enc_rev = list(revise_pnc(unlist(pnc_enc),
                                    pnc_wk)) ,
      
      
      ######################################
      ### Missing Indicators
      
      ltfu_sev = first_missing_sev(unlist(missing_by_sev),
                                   pnc_wk)
      
    ) %>% 
    unnest_wider(c(untreated_po, treated_po)) 
  
  
  
  # Assign the treatment values among those where include ne 0
  data2 <- assign_trt(subset(data, include != 0)) %>% 
    # Create variables for the observed outcomes - ignore missingness at this stage
    mutate(preeclampsia_pre_miss = ifelse(trt == 0,
                                          preeclampsia0,
                                          preeclampsia1),
           pregout_pre_miss = ifelse(trt == 0,
                                     final_preg0,
                                     final_preg1),
           pregout_t_pre_miss = ifelse(trt == 0,
                                       final_preg0_t,
                                       final_preg1_t)) %>% 
    # Now determine if missing outcome based upon outcome and gw of outcome
    mutate(p_out_ltfu = ifelse(pregout_pre_miss == 'fetaldeath' & 
                                 pregout_t_pre_miss < 20,
                               logodds_to_p(p_miss_outcome[1] + 
                                              (pregout_t_pre_miss * 
                                                 p_miss_outcome[2])),
                               0))
  
  # Assign an indicator for if someone is LTFU  due to their true outcome value
  # Create a MAR and MNAR scenario
  data3 <- assign_out_ltfu(data2) %>% 
    ## Create an indicator variable for if the person was LTFU and how
    # -- This would occur if either:
    # -- - ltfu_sev is not missing and less than or equal to pregout_t_pre_miss, OR
    # -- - ltfu_out is equal to 1
    # Create an indicator just for MAR only and MAR+MNAR.
    mutate(ltfu_mar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
                             "sev",
                             "not"),
           # ltfu_mar_mnar = ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
           #                        "sev",
           #                        ifelse(ltfu_out == 1,
           #                               "out",
           #                               "not"))) %>% 
           ltfu_mar_mnar = ifelse(ltfu_out == 1 & (is.na(ltfu_sev) || ltfu_sev > pregout_t_pre_miss),
                                  "out",
                                  ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
                                         "sev",
                                         "not"))) %>%
    # Now determine timing when LTFU
    mutate(t_ltfu_mar = pnc_miss(unlist(pnc_enc_rev),
                                 ltfu_mar,
                                 41, # Just set as a constant because ignoring missingness due to outcome
                                 ltfu_sev),
           t_ltfu_mar_mnar = pnc_miss(unlist(pnc_enc_rev),
                                      ltfu_mar_mnar,
                                      pregout_t_pre_miss,
                                      ltfu_sev)) %>% 
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
  # NOTE: Analysis problem to deal with -- some censored on the day that they are randomized (i.e., they have their indexing encounter and then never come back). Will likely need to add a small constant to all times. Check with Jess.
  
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
  first_preg_t = find_first_not_contpreg(preg_outcomes) #+ 1
  
  ## Determine if the person developed preeclampsia if untreated
  first_preec_t = find_first_preeclampsia(preeclampsia_outcomes) # + 1
  
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
                      preg_outcomes[final_preg_t], 
                      revised_outcomes[final_preg_t]) 
  
  final_preg = gsub("_next", "", final_preg)
  
  
  # Determine if the person would be indexed into the cohort at 4, 7, 
  # and 16 weeks of gestation -- A person had to NOT have had a fetal
  # death prior to that gestational week to be included
  # Make an indicator variable for if they should be included
  if (final_preg_t > 4  & pnc_wk == 4){
    include = 4
  } else if (final_preg_t > 7 & pnc_wk == 7){
    include = 7
  } else if (final_preg_t > 16 & pnc_wk == 16){
    include = 16
  } else{
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
  first_preec_t = find_first_preeclampsia(preeclampsia_outcomes) 
  # Note: Preeclampsia cannot develop until after the prenatal timings of interest
  
  ## Find the first pregnancy outcome that occurred after the indexing
  ## prenatal encounter -- Only way that we will see these outcomes
  first_preg_t = find_first_not_contpreg_wk(preg_outcomes, wk) 
  
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
                      preg_outcomes[final_preg_t], 
                      revised_outcomes[final_preg_t])
  
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
  
  n <- length(pnc_encounters)
  
  # No prenatal encounters prior to indexing prenatal encounter
  pre_week <- rep(0, wk-1)
  
  # Prenatal encounter on indexing prenatal encounter week
  week <- 1
  
  # Prenatal encounters after the indexing prenatal encounter
  wk_plus <- pnc_encounters[(wk+1):n]
  
  # Return the final list of prenatal encounters
  list(
    c(pre_week,week, wk_plus)
    )
  
}






### FUNCTION first_missing_sev()
# This function is intended to identify the
# first gestational week where the person
# is LTFU due to hypertension severity.
# INPUTS:
# - missing_vec: vector of the missing indicators by severity
# - gestwk: gestationanl week that prenatal care was initiated
#   have to be LTFU after that.

first_missing_sev <- function(missing_vec, gestwk){
  
  # Identify the instances where missing
  ind <- which(missing_vec == 1)
  
  # Gestational week of the first
  # pnc visit is indexed from week zero. However, the 
  # vector itself is indexed by R from 1, not zero.
  # Losses to follow-up occur the week after they were
  # determined, aligning with the R indexing
  ind2 <- which(ind > gestwk)[1]
  
  if (length(ind2) == 0) {
    return(NA)
  } else {
    # Do not subtract 1 because we assume that the loss to
    # follow-up occurs the next (prior to anything else occurring that week)
    return(ind[ind2])
  }
  
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







### FUNCTION: assign_out_ltfu()
# This function takes the vector of probabilities
# of being LTFU due to an outcome and 
# assigns a censoring indicator to each person in 
# a cohort dependent upon that probability.
# This is only done for observed outcomes.

assign_out_ltfu <- function(dataset){
  
  # Pull out the p_out_ltfu vector
  p_ltfu <- dataset$p_out_ltfu
  # p_ltfu <- data2$p_out_ltfu
  
  # Generate treatments
  ltfu_out <- rbinom(n = length(p_ltfu),
                     size = 1,
                     prob = p_ltfu)
  
  dataset$ltfu_out <- ltfu_out
  
  return(dataset)
  
}



### FUNCTION: pnc_miss()
# This function takes a vector of PNC encounters 
# and the ltfu indicator and determines the timing
# of the most recent PNC.
# Someone will be LTFU at their last PNC encounter
# PRIOR to that indicator.
pnc_miss <- function(pnc_enc, ltfu_ind, t_preg_out, t_sev){
  
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





# OLD VERSION CREATE_COHORT
# New version modified this old version so that all missingness 
# is created in the cohort identification stage.
# create_cohort <- function(dataset, p_miss_outcome){
#   
#   data <- dataset %>% # all_outcomes 
#     dplyr::group_by(sim_id, id) %>% 
#     #rowwise() %>% 
#     dplyr::mutate(
#       
#       ######################################
#       ### Untreated Potential Outcomes
#       
#       untreated_po = list(untreated_outcomes(unlist(preg_outcomes_untrt), 
#                                              unlist(preec_outcomes_untrt), 
#                                              unlist(revised_preg),
#                                              pnc_wk)),
#       
#       ######################################
#       ### Treated Potential Outcomes
#       
#       treated_po = list(treated_outcomes(unlist(preg_outcomes_trt),
#                                          unlist(preec_outcomes_trt),
#                                          unlist(revised_preg),
#                                          pnc_wk)), 
#       
#       ######################################
#       ### Prenatal Encounters
#       ### Identify prenatal encounters after
#       ### the first.
#       
#       # Need to map in two values from the tibble.
#       pnc_enc_rev = list(revise_pnc(unlist(pnc_enc),
#                                     pnc_wk)) ,
#       
#       
#       ######################################
#       ### Missing Indicators
#       
#       ltfu_sev = first_missing_sev(unlist(missing_by_sev),
#                                    pnc_wk),
#     ) %>% 
#     unnest_wider(c(untreated_po, treated_po)) 
#   
#   # Assign the treatment values among those where include ne 0
#   data2 <- assign_trt(subset(data, include != 0)) %>% 
#     # Create variables for the observed outcomes - ignore missingness at this stage
#     mutate(preeclampsia_pre_miss = ifelse(trt == 0,
#                                           preeclampsia0,
#                                           preeclampsia1),
#            pregout_pre_miss = ifelse(trt == 0,
#                                      final_preg0,
#                                      final_preg1),
#            pregout_t_pre_miss = ifelse(trt == 0,
#                                        final_preg0_t,
#                                        final_preg1_t)) %>% 
#     # Now determine if missing outcome based upon outcome and gw of outcome
#     mutate(p_out_ltfu = ifelse(pregout_pre_miss == 'fetaldeath' & 
#                                  pregout_t_pre_miss < 20,
#                                logodds_to_p(p_miss_outcome[1] + 
#                                               (pregout_t_pre_miss * 
#                                                  p_miss_outcome[2])),
#                                0))
#   
#   # Assign an indicator for if someone is LTFU  due to their true outcome value
#   data3 <- assign_out_ltfu(data2) %>% 
#     ## Create an indicator variable for if the person was LTFU and how
#     # -- This would occur if either:
#     # -- - ltfu_sev is not missing and less than or equal to pregout_t_pre_miss, OR
#     # -- - ltfu_out is equal to 1
#     mutate(ltfu = ifelse(ltfu_out == 1 & (is.na(ltfu_sev) || ltfu_sev > pregout_t_pre_miss),
#                          "out",
#                          ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
#                                 "sev",
#                                 "not"))) %>% 
#     # Now determine timing when LTFU
#     mutate(t_ltfu = pnc_miss(unlist(pnc_enc_rev),
#                              ltfu,
#                              pregout_t_pre_miss,
#                              ltfu_sev)) %>% 
#     # Finally, create their observed outcomes, incorporating the LTFU
#     mutate(preeclampsia = ifelse(ltfu != 'not',
#                                  0,
#                                  preeclampsia_pre_miss),
#            pregout = ifelse(ltfu != 'not',
#                             'unknown',
#                             pregout_pre_miss),
#            pregout_t = ifelse(ltfu != 'not',
#                               t_ltfu,
#                               pregout_t_pre_miss))
#   # NOTE: Analysis problem to deal with -- some censored on the day that they are randomized (i.e., they have their indexing encounter and then never come back). Will likely need to add a small constant to all times. Check with Jess.
#   
#   return(data3)
#   
# }







# OLD: No longer structured this way.
#### FUNCTION: create_study_samples_by_gw() -- OLD
# This function will identify study samples (prior to inducing missingness)
# for each of the gestational weeks of PNC initiation.

# create_study_samples_by_gw <- function(dataset){
#   
#   # Make a nested dataset for each of the cohorts
#   
#   gw4 <- dataset %>% 
#     filter(include_4wk) %>% 
#     select(-c(include_7wk, include_16wk, treated_po_wk7, 
#               treated_po_wk16, pnc_enc_wk7, pnc_enc_wk16,
#               missing_sev_wk7, missing_sev_wk16)) %>% 
#     unnest_wider(treated_po_wk4) %>% 
#     nest(wk4 = -sim_id)
#   
#   gw7 <- dataset %>% 
#     filter(include_7wk) %>% 
#     select(-c(include_4wk, include_16wk, treated_po_wk4, 
#               treated_po_wk16, pnc_enc_wk4, pnc_enc_wk16,
#               missing_sev_wk4, missing_sev_wk16)) %>% 
#     unnest_wider(treated_po_wk7) %>% 
#     nest(wk7 = -sim_id)
#   
#   gw16 <- dataset %>% 
#     filter(include_16wk) %>% 
#     select(-c(include_4wk, include_7wk, treated_po_wk4, 
#               treated_po_wk7, pnc_enc_wk4, pnc_enc_wk7,
#               missing_sev_wk4, missing_sev_wk7)) %>% 
#     unnest_wider(treated_po_wk16) %>% 
#     nest(wk16 = -sim_id)
#   
#   # Left join these datasets on sim_id
#   sim <- left_join(gw4, gw7, by = 'sim_id') %>% 
#     left_join(gw16, by = 'sim_id')
#   
#   return(sim)
#   
# }
