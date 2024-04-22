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
  
  data <- all_outcomes %>% 
    #dplyr::group_by(sim_id, id) %>% 
    rowwise() %>% 
    dplyr::mutate(
      
      ## Find the first pregnancy outcome if untreated
      # Gestational week that the outcome occurred (i.e., 1 + selected week)
      first_py0_t_v1 = sapply(preg_outcomes_untrt, find_first_not_contpreg) + 1,
      # Get the value of the pregnancy outcome
      first_py0_v1 = sapply(preg_outcomes_untrt, 
                         function(x) ifelse(is.na(find_first_not_contpreg(x)), 
                                            NA, 
                                            x[find_first_not_contpreg(x)])),
      
      ## Determine if the person developed preeclampsia if untreated
      first_preec0_t_v1 = sapply(preec_outcomes_untrt, find_first_preeclampsia) + 1,
      ## Get the outcome if preeclampsia
      
      
      ### RETAIN These outcomes.
      ## A person only actually develops preeclampsia if it occurs prior
      ## to the untreated pregnancy outcome.
      preeclampsia_0 = ifelse(is.na(first_preec0_t_v1),
                              0,
                              as.numeric(first_preec0_t_v1 < first_py0_t_v1)),
      
      ## Now, get the revised outcome based upon their preeclampsia value
      # Gestational week that the outcome occurs
      first_py0_t = ifelse(preeclampsia_0 == 0, # If no preeclampsia
                           first_py0_t_v1, # Then potential outcome
                           first_preec0_t_v1) , # Else preecalmpsia timing
      # Get the actual outcome value
      first_py0 = if_else(preeclampsia_0 == 0,
                         first_py0_v1,
                         sapply(revised_preg, function(x) x[first_py0_t - 1])),
      # Remove _next from the end of the variable
      first_py0 = gsub("_next", "", first_py0),

      # Determine if the person would be indexed into the cohort at 4, 7, 
      # and 16 weeks of gestation -- A person had to NOT have had a fetal
      # death prior to that gestational week to be included
      include_4wk = first_py0_t > 4,
      include_7wk = first_py0_t > 7,
      include_16wk = first_py0_t > 16,

      
      ## Determine pregnancy outcomes if they had been treated at 4 weeks.
      
      
    )
  
}




### Some necessary helper functions

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

