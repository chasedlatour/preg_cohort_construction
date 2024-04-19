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
    mutate(
      
      #COME BACK - THIS IS NOT WORKING
      
      # First untreated potential outcome (ignoring revised outcomes after preec)
      first_fetaldeath0 = match("fetaldeath_next", preg_outcomes_untrt)
      # first_fetaldeath = sapply(preg_outcomes_untrt, function(x){
      #   which()
      # } )
      
      # Indicator for if included in the cohort at 6 weeks of gestation
      #included_6 = 
    )
  
}