#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: The program contains all the R functions
# for generating the data.
#####################################################



# Testing

# n_sim <- 1
# n <- 100
# p_sev_beta <- c(-3, 0.1, 0.2)
# p_trt_sev <- c(0.35, 0.50, 0.65)
# p_indx_pnc <- c(0.25, 0.375, 0.375) # Think about what these probabilities should be -- likely just want reasonable distribution












#### FUNCTION: generate()
# This is the baseline function that will be 
# -- used to generate the potential outcomes
# -- for each cohort.

# Inputs:
# n_sim = simulation number (id)
# n = number of individuals in a cohort/repetition
# p_trt_sev = vector of treatment probabilities by disease severity
# p_indx_pnc = vector of probabilities for week of indexing prenatal care visit (6, 9, 18 from LMP)
# potential_preg_untrt = dataset with the weekly probabilities for untreated potential preg outcomes
# potential_preg_trt = same but treated potential pregnancy outcomes
# potential_preec_untrt = dataset with weekly probabilities for untreated potential preeclampsia outcomes
# potential_preec_trt = same but treated potential preeclampsia outcomes
# revised_preg = dataset with the weekly probabilities for pregnancy outcomes after preeclampsia
# pnc_prob = dataset with the weekly probabilities for a prenatal encounter

generate <- function(n_sim, n, p_trt_sev, p_indx_pnc,
                     potential_preg_untrt, potential_preg_trt, 
                     potential_preec_untrt, potential_preec_trt,
                     revised_preg, pnc_prob){
  
  # Per Morris et al. 2019. Save the random state at the beginning of 
  # the simulation in case want it later
  initial_seed <- list(.Random.seed)
  
  # Get the severity distribution from a multinomial random variable
  # We assume that there is an equal distribution across three levels. 
  # Otherwise, this needs to be modified.
  severity_dist = rmultinom(n=1, size=n, prob=c(1/3, 1/3, 1/3))
  
  
  ######## GENERATE PEOPLE AND BASELINE VALUES
  
  # Create the dataset that going to output
  data <- dplyr::tibble(
    
    # Simulation ID -- Only important if we simulate more than one simulation cohort
    sim_id = n_sim, 
    
    # Each individual's ID
    id = 1:n, 
    
    # Generate 0 through 40 gestational weeks - 1 vector
    # This indexed prenatal care encounters, not outcomes. Outcomes occur at 1+ this.
    # Note: R indexes vectors from 1, not 0, as we have done for gestational weeks from conception
    pnc_gw = list(seq(0,40, by = 1)),
    
    # Treatment needs to be assigned AFTER the cohort is selected. Otherwise,
    # there is too much imbalance and don't get asymptotic convergence of true
    # and observed treatment effect. -- Thus, done later.
    
    # Assign a person's hypertension severity. 
    severity = c(rep(0, as.numeric(severity_dist[1,1])), 
                 rep(1, severity_dist[2,1]),
                 rep(2, severity_dist[3,1])),
    # 0 = Low, 1 = Moderate, 2 = High Severity
    
    # Assign a treatment probability to the person dependent upon
    # their disease severity
    # Severity indexed as 0, 1, 2, but R indexes vectors from 1. Thus, added 1
    p_trt = p_trt_sev[severity+1],
    
    ##### GENERATE POTENTIAL PREGNANCY OUTCOMES
    
    ## Step 2: Untreated
    preg_outcomes_untrt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_untrt, .x)),
    
    ## Step 3: Treated
    preg_outcomes_trt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_trt, .x)),
    
    ##### GENERATE POTENTIAL PREECLAMPSIA OUTCOMES
    
    ## Step 4: Untreated
    preec_outcomes_untrt = purrr::map(severity, ~sample_preeclampsia(potential_preec_untrt, .x)),
    
    ## Step 5: Treated
    preec_outcomes_trt = purrr::map(severity, ~sample_preeclampsia(potential_preec_trt, .x)),
    
    ##### GENERATE REVISED OUTCOMES IF PREECLAMPSIA -- Step 6
    
    revised_preg = purrr::map(severity, ~sample_revised_preg(revised_preg, .x)),
    
    ##### GENERATE PNC ENCOUNTERS -- Step 7
    
    # Determine gestational week of indexing prenatal encounter
    pnc_wk = sample(x = c(4, 7, 16), size = n, prob = p_indx_pnc, replace = TRUE),
    
    # List of indicator variables for PNC encounters
    # We incorporate indexing prenatal encounter later
    pnc_enc = purrr::map(severity, ~sample_pnc(pnc_prob, .x))
    
  )
  
  # Final random state
  finish_seed <- list(.Random.seed)
  
  ## Store the random states in the final dataset.
  data <- data %>% 
    mutate(
      start_seed = initial_seed,
      end_seed = finish_seed
    )
  
  return(data)
  
}




      

### FUNCTION: sample_outcomes_for_id()
# This function will create the potential pregnancy outcomes based upon 
# the probabilities recorded in the Excel file.
# This will output a list of the potential pregnancy outcomes for one person.
#
# Inputs:
# - data = dataset with the weekly probabilities for potential pregnancy outcomes
# - sev = severity value for that person
sample_outcomes_for_id <- function(data, sev) {
  
  # Subset to the severity level of interest
  data2 <- subset(data, severity == sev)

  # Apply this function to each row of probability dataset
  outcomes <- sapply(1:nrow(data2), function(i) {

    # Indicate the potential options to select from
    options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")

    # Assign the probabilities to each of the potential pregnancy outcome options based upon
    # -- the corresponding rows in the probability dataset
    probabilities <- data2[i, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]

    # Sample an option based on probabilities
    sampled_option <- sample(options, size = 1, prob = probabilities)

    # Return that sampled option. This will be stored in the vector outcomes
    return(sampled_option)
  })

  return(list(outcomes))
}







### FUNCTION sample_preeclampsia()
# This function will create potential preeclampsia outcomes 
# This will output a list of the potential preeclampsia
# outcomes for one person.
#
# INPUT:
# - data = dataset with the weekly potential preeclampsia probabilities
# - sev = severity value for that person

sample_preeclampsia <- function(data, sev) {
  
  # Subset to the severity level of interest 
  # Get the probabilities for each gestational week
  probabilities <- subset(data, severity == sev)$p_preeclampsia
  
  # Use these probabilities with a bernoulli/binomial rv
  out <- rbinom(n = length(probabilities), 
                size=1,
                prob = probabilities)
  
  return(list(out))
}







#### FUNCTION: process_out()
# This function takes in the vector out from 
# sample_revised_preg() -- below -- and replaces its values
# with those that are compatible with the other vectors.
#
# INPUT:
# - out = vector from sample_revised_preg

process_out <- function(out, n_NA) {
  # Create an empty vector of characters with the same length as 'probabilities'
  results <- vector("character", length(out))
  
  # Set the first 17 elements to "NA"
  results[1:n_NA] <- NA
  
  # For the remaining elements, map 1 to "fetaldeath" and 0 to "delivery"
  results[(n_NA+1):length(out)] <- ifelse(out[(n_NA+1):length(out)] == 1,
                                          "fetaldeath_next", 
                                          "livebirth_next")
  
  return(results)
}

# Example use:
# probabilities <- c(rep(0, 10), rep(1, 10), rep(0, 5), rep(1, 5))
# output <- process_out(probabilities)
# print(output)



### FUNCTION sample_revised_preg()
# This function will provide revised pregnancy
# outcomes for the cases in which a pregnancy
# ends in preeclampsia.
sample_revised_preg <- function(data, sev){
  
  # Subset to the severity level of interest 
  # Get the probabilities for each gestational week
  probabilities <- subset(data, severity == sev)$p_fetaldeath
  
  # Determine how many NAs there should be
  n_NA <- length(probabilities[probabilities == 0])
  
  # Use these probabilities with a bernoulli/binomial rv
  out <- rbinom(n = length(probabilities), 
                size=1,
                prob = probabilities)
  # out is a numeric vector - we want a character vector with the
  # same values as the others
  
  # Revised here using above step.
  rev_out <- process_out(out, n_NA)
  
  return(list(rev_out))
  
}






### FUNCTION: sample_pnc()
# This function will provide prenatal care encounter
# indicators.
sample_pnc <- function(data, sev) {
  
  # Subset to the severity level of interest 
  # Get the probabilities for each gestational week
  probabilities <- subset(data, severity == sev)$p_pnc
  
  # Use these probabilities with a bernoulli/binomial rv
  out <- rbinom(n = length(probabilities), 
                size=1,
                prob = probabilities)
  
  return(list(out))
}






# # OLD HELPER FUNCTION
# ### FUNCTION: logodds_to_p()
# # This function converts the log-odds into a
# # probability. Important for deriving probabilities
# # from the logit function.
# logodds_to_p <- function(logodds){
#   
#   p <- exp(logodds)/(1 + exp(logodds))
#   
#   return(p)
#   
# }








### OLD VERSION OF GENERATE()
# This function was written so that LTFU due to severity was generated
# in this stage. However, it is less computationally intense to only
# generate the DGM once and then incorporate all missingness later.
# generate <- function(n_sim, n, p_sev_beta, p_trt_sev, p_indx_pnc){
#   
#   # Per Morris et al. 2019. Save the random state at the beginning of the simulation in case want it later.
#   initial_seed <- list(.Random.seed)
#   
#   # Get the severity distribution from a multinomial random variable
#   # We assume that there is an equal distribution across three levels. 
#   # Otherwise, this needs to be modified.
#   severity_dist = rmultinom(n=1, size=n, prob=c(1/3, 1/3, 1/3))
#   
#   # Get the probabilities associated with the betas for each
#   # severity level - calculated via logistic regression
#   p_missing_sev <- c(logodds_to_p(p_sev_beta[1]),
#                      logodds_to_p(p_sev_beta[1] + p_sev_beta[2]),
#                      logodds_to_p(p_sev_beta[1] + p_sev_beta[3]))
#   
#   
#   ######## GENERATE PEOPLE AND BASELINE VALUES
#   
#   # Create the dataset that going to output
#   data <- dplyr::tibble(
#     
#     # Simulation ID -- Only important if we simulate more than one simulation cohort
#     sim_id = n_sim, 
#     
#     # Each individual's ID
#     id = 1:n, 
#     
#     # Generate 0 through 40 gestational weeks - 1 vector
#     # This indexed prenatal care encounters, not outcomes. Outcomes occur at 1+ this.
#     # Note: R indexes vectors from 1, not 0, as we have done for these gestational weeks
#     pnc_gw = list(seq(0,40, by = 1)),
#     
#     # Treatment needs to be assigned AFTER the cohort is selected. Otherwise,
#     # there is too much imbalance and don't get asymptotic convergence of true
#     # and observed treatment effect. -- Thus, done later.
#     
#     # Assign a person's hypertension severity. 
#     severity = c(rep(0, as.numeric(severity_dist[1,1])), 
#                  rep(1, severity_dist[2,1]),
#                  rep(2, severity_dist[3,1])),
#     # 0 = Low, 1 = Moderate, 2 = High Severity
#     
#     # Assign a treatment probability to the person dependent upon
#     # their disease severity
#     # Severity indexed as 0, 1, 2, but R indexes vectors from 1. Thus, added 1
#     p_trt = p_trt_sev[severity+1],
#     
#     # Generate missingness probabilities for each person
#     p_missing_by_sev = p_missing_sev[severity+1],
#     
#     # Select missingness values for each person's gestational week
#     # manually input as 41
#     missing_by_sev = purrr::map(p_missing_by_sev,
#                                 ~rbinom(41, 1, .x)),
#     
#     ##### GENERATE POTENTIAL PREGNANCY OUTCOMES
#     
#     ## Step 2: Untreated
#     preg_outcomes_untrt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_untrt, .x)),
#     
#     ## Step 3: Treated
#     preg_outcomes_trt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_trt, .x)),
#     
#     ##### GENERATE POTENTIAL PREECLAMPSIA OUTCOMES
#     
#     ## Step 4: Untreated
#     preec_outcomes_untrt = purrr::map(severity, ~sample_preeclampsia(potential_preec_untrt, .x)),
#     
#     ## Step 5: Treated
#     preec_outcomes_trt = purrr::map(severity, ~sample_preeclampsia(potential_preec_trt, .x)),
#     
#     ##### GENERATE REVISED OUTCOMES IF PREECLAMPSIA -- Step 6
#     
#     revised_preg = purrr::map(severity, ~sample_revised_preg(revised_preg, .x)),
#     
#     ##### GENERATE PNC ENCOUNTERS -- Step 7
#     
#     # Determine gestational week of indexing prenatal encounter
#     pnc_wk = sample(x = c(4, 7, 16), size = n, prob = p_indx_pnc, replace = TRUE),
#     
#     # List of indicator variables for PNC encounters
#     # We incorporate indexing prenatal encounter later
#     pnc_enc = purrr::map(severity, ~sample_pnc(pnc_prob, .x))
#     
#   )
#   
#   # Final random state
#   finish_seed <- list(.Random.seed)
#   
#   ## Store the random states in the final dataset.
#   data <- data %>% 
#     mutate(
#       start_seed = initial_seed,
#       end_seed = finish_seed
#     )
#   
#   return(data)
#   
# }





