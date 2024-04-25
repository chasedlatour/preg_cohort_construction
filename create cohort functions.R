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
identify_key_vars <- function(dataset){
  
  # REPLACE W DATASET WHEN DONE EDITING
  data <- dataset %>% 
    dplyr::group_by(sim_id, id) %>% 
    #rowwise() %>% 
    dplyr::mutate(
      
      ######################################
      ### Untreated Potential Outcomes
      
      untreated_po = list(untreated_outcomes(unlist(preg_outcomes_untrt), 
                                             unlist(preec_outcomes_untrt), 
                                             unlist(revised_preg))),
      
      ######################################
      ### Treated Potential Outcomes
      
      treated_po_wk4 = list(treated_outcomes(unlist(preg_outcomes_trt),
                                             unlist(preec_outcomes_trt),
                                             unlist(revised_preg),
                                             4)),
      
      treated_po_wk7 = list(treated_outcomes(unlist(preg_outcomes_trt),
                                             unlist(preec_outcomes_trt),
                                             unlist(revised_preg),
                                             7)),
      
      treated_po_wk16 = list(treated_outcomes(unlist(preg_outcomes_trt),
                                             unlist(preec_outcomes_trt),
                                             unlist(revised_preg),
                                             16)),
      
      ######################################
      ### Prenatal Encounters
      ### Identify prenatal encounters after
      ### the first.
      
      pnc_enc_wk4 = list(revise_pnc(unlist(pnc_enc),
                                    4)),
      
      pnc_enc_wk7 = list(revise_pnc(unlist(pnc_enc),
                                    7)),
      
      pnc_enc_wk16 = list(revise_pnc(unlist(pnc_enc),
                                     16)),
      
      ######################################
      ### Missing Indicators
      
      missing_sev_wk4 = purrr::map_dbl(missing_by_sev, ~first_missing_sev(.x, 4)),
      
      missing_sev_wk7 = purrr::map_dbl(missing_by_sev, ~first_missing_sev(.x, 7)),
      
      missing_sev_wk16 = purrr::map_dbl(missing_by_sev, ~first_missing_sev(.x, 16))
      
    ) %>% 
    unnest_wider(untreated_po) %>% 
    create_study_samples_by_gw()
      
      
    return(data)
  
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
untreated_outcomes <- function(preg_outcomes, preeclampsia_outcomes, revised_outcomes){
  
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
  
  # Determine if the person would be indexed into the cohort at 4, 7, 
  # and 16 weeks of gestation -- A person had to NOT have had a fetal
  # death prior to that gestational week to be included
  include_4wk = final_preg_t > 4 
  include_7wk = final_preg_t > 7
  include_16wk = final_preg_t > 16
  # Don't need to change the indexes with +/- 1 here because outcomes already indexed
  # up with the R vectors.
  
  
  # Return a list of the desired values
  list(
    preeclampsia0 = preeclampsia,
    final_preg0_t = final_preg_t,
    final_preg0 = final_preg,
    include_4wk = include_4wk,
    include_7wk = include_7wk,
    include_16wk = include_16wk
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
  wk_plus <- pnc_encounters[wk+1:n]
  
  return(c(pre_wk, wk_plus))
  
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
  #ind2 <- ind[ind > gestwk+1][1]
  ind2 <- which(ind > gestwk)[1]
  
  if (length(ind2) == 0) {
    return(NA)
  } else {
    # Do not subtract 1 because we assume that the loss to
    # follow-up occurs the next (prior to anything else occurring that week)
    return(ind[ind2])
  }
  
}






#############################################
###OLD CODE
# 
# identify_key_vars <- function(dataset){
#   
#   data <- all_outcomes %>% 
#     #dplyr::group_by(sim_id, id) %>% 
#     rowwise() %>% 
#     dplyr::mutate(
#       
#       ######################################
#       ### Untreated Potential Outcomes
#       
#       ## Find the first pregnancy outcome if untreated
#       # Gestational week that the outcome occurred (i.e., 1 + selected week)
#       first_py0_t_v1 = sapply(preg_outcomes_untrt, find_first_not_contpreg) + 1,
#       # Get the value of the pregnancy outcome
#       # first_py0_v1 = sapply(preg_outcomes_untrt, 
#       #                    function(x) ifelse(is.na(find_first_not_contpreg(x)), 
#       #                                       NA, 
#       #                                       x[find_first_not_contpreg(x)])),
#       
#       ## Determine if the person developed preeclampsia if untreated
#       first_preec0_t_v1 = sapply(preec_outcomes_untrt, find_first_preeclampsia) + 1,
#       ## Get the outcome if preeclampsia
#       
#       
#       ### RETAIN These outcomes.
#       ## A person only actually develops preeclampsia if it occurs prior
#       ## to the untreated pregnancy outcome.
#       preeclampsia_0 = ifelse(is.na(first_preec0_t_v1),
#                               0,
#                               as.numeric(first_preec0_t_v1 < first_py0_t_v1)),
#       
#       ## Now, get the revised outcome based upon their preeclampsia value
#       # Gestational week that the outcome occurs
#       first_py0_t = ifelse(preeclampsia_0 == 0, # If no preeclampsia
#                            first_py0_t_v1, # Then potential outcome
#                            first_preec0_t_v1) , # Else preecalmpsia timing
#       
#       # Get the actual outcome value
#       first_py0 = if_else(preeclampsia_0 == 0,
#                           sapply(preg_outcomes_untrt,
#                                  function(x) x[first_py0_t_v1 - 1]), #first_py0_v1,
#                           sapply(revised_preg, function(x) x[first_py0_t - 1])),
#       
#       # Remove _next from the end of the variable
#       first_py0 = gsub("_next", "", first_py0),
#       
#       # Determine if the person would be indexed into the cohort at 4, 7, 
#       # and 16 weeks of gestation -- A person had to NOT have had a fetal
#       # death prior to that gestational week to be included
#       include_4wk = first_py0_t > 4,
#       include_7wk = first_py0_t > 7,
#       include_16wk = first_py0_t > 16,
#       
#       ##################################
#       
#       
#       
#       #################################
#       ### Treated Potential Outcomes
#       
#       ## Determine if the person developed preeclampsia if treated
#       first_preec1_t_v1 = sapply(preec_outcomes_trt, find_first_preeclampsia) + 1,
#       # Note: This is done only once because preeclampsia cannot develop until after
#       # the prenatal encounter timings of-interest
#       
#       
#       ### IF TREATED AT 4 WEEKS
#       
#       # Gestational week that the outcome occurred (i.e., 1 + selected week)
#       first_py1_t_v1_wk4 = sapply(preg_outcomes_trt, function(x)
#         find_first_not_contpreg_wk(x, 4)) + 1,
#       
#       ## A person only actually develops preeclampsia if it occurs prior
#       ## to the treated pregnancy outcome.
#       preeclampsia_1_wk4 = ifelse(is.na(first_preec1_t_v1),
#                                   0,
#                                   as.numeric(first_preec1_t_v1 < first_py1_t_v1_wk4)),
#       
#       ## Now, get the revised outcome based upon their preeclampsia value
#       # Gestational week that the outcome occurs
#       first_py1_t_wk4 = ifelse(preeclampsia_1_wk4 == 0, # If no preeclampsia
#                                first_py1_t_v1_wk4, # Then potential outcome
#                                first_preec1_t_v1), # Else preecalmpsia timing
#       # Get the actual outcome value
#       first_py1_wk4 = if_else(preeclampsia_0 == 0,
#                               sapply(preg_outcomes_trt, 
#                                      function(x) x[first_py1_t_v1_wk4 - 1]),
#                               sapply(revised_preg, 
#                                      function(x) x[first_py1_t_wk4 - 1])) #,
#       # # Remove _next from the end of the variable
#       # first_py0 = gsub("_next", "", first_py0)
#       
#       
#     )
#   
# }
