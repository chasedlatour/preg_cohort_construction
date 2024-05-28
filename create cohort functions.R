#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: This program contains all the R functions
# for identifying the study cohorts.
#####################################################


#### FUNCTION: identify_key_vars()
# This function will identify necessary values of key variables for each person
# 
# INPUT:
# dataset - data generated from the function generate()
# p_miss_outcome - vector of gamma 0 and gamma 1 for the missing preg outcomes according to outcome


### CHASE -- RUN SOME SUMMARIES OF EACH OF THESE DATASETS TO MAKE SURE THEY LOOK AS EXPECTED

create_cohort <- function(dataset, p_miss_outcome){
  
  data <- dataset %>% # all_outcomes 
    dplyr::group_by(sim_id, id) %>% 
    #rowwise() %>% 
    dplyr::mutate(
      
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
      
      # Need to map in two values from the tibble.
      pnc_enc_rev = list(revise_pnc(unlist(pnc_enc),
                                    pnc_wk)) ,

      
      ######################################
      ### Missing Indicators
      
      ltfu_sev = first_missing_sev(unlist(missing_by_sev),
                                      pnc_wk),
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
  data3 <- assign_out_ltfu(data2) %>% 
    ## Create an indicator variable for if the person was LTFU and how
    # -- This would occur if either:
    # -- - ltfu_sev is not missing and less than or equal to pregout_t_pre_miss, OR
    # -- - ltfu_out is equal to 1
    mutate(ltfu = ifelse(ltfu_out == 1 & (is.na(ltfu_sev) || ltfu_sev > pregout_t_pre_miss),
                         "out",
                         ifelse(!is.na(ltfu_sev) & ltfu_sev <= pregout_t_pre_miss,
                                "sev",
                                "not"))) %>% 
    # Now determine timing when LTFU
    mutate(t_ltfu = pnc_miss(unlist(pnc_enc_rev),
                             ltfu,
                             pregout_t_pre_miss,
                             ltfu_sev)) %>% 
    # Finally, create their observed outcomes, incorporating the LTFU
    mutate(preeclampsia = ifelse(ltfu != 'not',
                                 0,
                                 preeclampsia_pre_miss),
           pregout = ifelse(ltfu != 'not',
                            'unknown',
                            pregout_pre_miss),
           pregout_t = ifelse(ltfu != 'not',
                              t_ltfu,
                              pregout_t_pre_miss))
    # NOTE: Analysis problem to deal with -- some censored on the day that they are randomized (i.e., they have their indexing encounter and then never come back). Will likely need to add a small constant to all times. Check with Jess.
      
    return(data3)
  
}

### Some necessary helper functions


#### FUNCTION: create_study_samples_by_gw()
# This function will identify study samples (prior to inducing missingness)
# for each of the gestational weeks of PNC initiation.

create_study_samples_by_gw <- function(dataset){
  
  # Make a nested dataset for each of the cohorts
  
  gw4 <- dataset %>% 
    filter(include_4wk) %>% 
    select(-c(include_7wk, include_16wk, treated_po_wk7, 
              treated_po_wk16, pnc_enc_wk7, pnc_enc_wk16,
              missing_sev_wk7, missing_sev_wk16)) %>% 
    unnest_wider(treated_po_wk4) %>% 
    nest(wk4 = -sim_id)
  
  gw7 <- dataset %>% 
    filter(include_7wk) %>% 
    select(-c(include_4wk, include_16wk, treated_po_wk4, 
              treated_po_wk16, pnc_enc_wk4, pnc_enc_wk16,
              missing_sev_wk4, missing_sev_wk16)) %>% 
    unnest_wider(treated_po_wk7) %>% 
    nest(wk7 = -sim_id)
  
  gw16 <- dataset %>% 
    filter(include_16wk) %>% 
    select(-c(include_4wk, include_7wk, treated_po_wk4, 
              treated_po_wk7, pnc_enc_wk4, pnc_enc_wk7,
              missing_sev_wk4, missing_sev_wk7)) %>% 
    unnest_wider(treated_po_wk16) %>% 
    nest(wk16 = -sim_id)
  
  # Left join these datasets on sim_id
  sim <- left_join(gw4, gw7, by = 'sim_id') %>% 
    left_join(gw16, by = 'sim_id')
  
  return(sim)
  
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

### FUNCITON: find_first_preeclampsia()
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
  index_greater_than_wk <- which(idx > wk)[1] #Previously: wk-1
  
  if (length(index_greater_than_wk) == 0) {
    return(NA)
  } else {
    return(idx[index_greater_than_wk])
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
  ## to the untreated pregnancy outcome.
  preeclampsia = ifelse(is.na(first_preec_t),
                        0,
                        as.numeric(first_preec_t < first_preg_t))
  
  ## Now, get the revised outcome based upon their preeclampsia value
  # Gestational week that the outcome occurs
  final_preg_t = ifelse(preeclampsia == 0, # If no preeclampsia
                        first_preg_t, # Then potential outcome
                        first_preec_t) # Else preecalmpsia timing
  
  # Get the actual outcome value
  final_preg = ifelse(preeclampsia == 0,
                      preg_outcomes[final_preg_t], # Previously had - 1 in the index
                      revised_outcomes[final_preg_t]) # Previously had - 1
  
  final_preg = gsub("_next", "", final_preg)
  
  # Old indicator variables
  # include_4wk = final_preg_t > 4 
  # include_7wk = final_preg_t > 7
  # include_16wk = final_preg_t > 16
  
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
  # Don't need to change the indexes with +/- 1 here because outcomes already indexed
  # up with the R vectors.
  
  
  # Return a list of the desired values
  list(
    preeclampsia0 = preeclampsia,
    final_preg0_t = final_preg_t,
    final_preg0 = final_preg,
    include = include # Inclusion variable indicator
  )
  
}









#### FUNCTION: treated_outcomes()
# This function the treated
# potential outcomes
treated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes, wk){
  
  ## Determine if the person developed preeclampsia if treated
  first_preec_t = find_first_preeclampsia(preeclampsia_outcomes) #+ 1
  # Note: Preeclampsia cannot develop until after the prenatal timings of interest
  
  ## Find the first pregnancy outcome
  # Gestational week that the outcome occurred (i.e., 1 + selected week)
  first_preg_t = find_first_not_contpreg_wk(preg_outcomes, wk) #+ 1
  
  ## A person only actually develops preeclampsia if it occurs prior
  ## to the treated pregnancy outcome.
  preeclampsia = ifelse(is.na(first_preec_t),
                        0,
                        as.numeric(first_preec_t < first_preg_t))
  
  ## Now, get the revised outcome based upon their preeclampsia value
  # Gestational week that the outcome occurs
  final_preg_t = ifelse(preeclampsia == 0, # If no preeclampsia
                        first_preg_t, # Then potential outcome
                        first_preec_t) # Else preecalmpsia timing
  
  # Get the actual outcome value
  final_preg = ifelse(preeclampsia == 0,
                      preg_outcomes[final_preg_t], #  - 1 in the []
                      revised_outcomes[final_preg_t]) #  - 1
  
  final_preg = gsub("_next", "", final_preg)
  
  
  # Return a list of the desired values
  list(
    preeclampsia1 = preeclampsia,
    final_preg1_t = final_preg_t,
    final_preg1 = final_preg
    # This reset the names but did not save them as expected for unnesting.
    # setNames(list(preeclampsia), paste0("preeclampsia1_wk", wk)),
    # setNames(list(final_preg_t), paste0("final_preg1_t_wk", wk)),
    # setNames(list(final_preg), paste0("final_preg1_wk", wk))
  )
  
}



#### FUNCTION: revise_pnc()
# Revise the occurrence of prenatal care
# encounters depending upon their gestational
# week of indexing (i.e., their first prenatal)
# care encounter
revise_pnc <- function(pnc_encounters, wk){
  
  n <- length(pnc_encounters)
  
  # Note that Excel starts from week 0 (i.e., pnc encounter during week 0)
  # As such, to have no prior PNC visits prior to week, we would make all 
  # PNC encounter indicators including wk 0 -- R starts indexing at 1
  pre_wk <- rep(0, wk)
  
  # Adding 1 because out gestweeks start at 0  in the Excel file
  # -- So, wk=7 in the R indexing aligns with gestwk 6 in the Excel file
  wk_plus <- pnc_encounters[(wk+1):n]
  
  # Revised
  # return(c(pre_wk, wk_plus))
  list(
    c(pre_wk, wk_plus)
    )
  
}


### FUNCTION first_missing_sev()
# This function is intended to identify the
# first gestational week where the person
# is LTFU due to hypertension severity.
first_missing_sev <- function(missing_vec, gestwk){
  
  # Identify the instances where missing
  ind <- which(missing_vec == 1)
  
  # Do not add one. Gestational week of the first
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
  
  # Generate treatments
  trt <- rbinom(n = length(p_trt),
                size = 1,
                prob = p_trt)
  
  dataset$trt <- trt
  
  return(dataset)
  
}



### FUNCTION: logodds_to_p()
# This function converts the log-odds into a
# probability. Important for deriving probabilities
# from the logit function.
logodds_to_p <- function(logodds){
  
  p <- exp(logodds)/(1 + exp(logodds))
  
  return(p)
  
}



### FUNCTION: assign_out_ltfu()
# This function takes the vector of probabilities
# of being LTFU due to an outcome and 
# assigns a censoring indicator to each person in 
# a cohort dependent upon that probability.
# This is only done for observed outcomes.

assign_out_ltfu <- function(dataset){
  
  # Pull out the p_trt vector
  p_ltfu <- dataset$p_out_ltfu
  
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
