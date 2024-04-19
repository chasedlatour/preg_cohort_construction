#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: The program contains all the R functions
# for generating the data.
#####################################################



#### FUNCTION: generate()
# This is the baseline function that will be 
# -- used to generate the potential outcomes
# -- for each cohort.
generate <- function(n_sim, n){
  
  # Per Morris et al. 2019, we save the random state at the beginning of the simulation
  initial_seed <- list(.Random.seed)
  
  # Generate 0 through 40 gestational weeks - 1 vector
  gw = list(seq(0,40, by = 1))
  
  # Get the severity distribution from a multinomial random variable
  # We assume that there is an equal distribution across three levels. 
  # Otherwise, this needs to be modified.
  severity_dist = rmultinom(n=1, size=n, prob=c(1/3, 1/3, 1/3))
  
  # Create the dataset that going to output
  data <- dplyr::tibble(
    
    # Simulation ID -- Only important if we simulate more than one simulation cohort
    sim_id = n_sim, 
    
    # Each individual's ID
    id = 1:n, 
    
    # Treatment needs to be assigned AFTER the cohort is selected. Otherwise,
    # there is too much imbalance and don't get asymptotic convergence of true
    # and observed treatment effect.
    
    # Assign a person's hypertension severity. 
    severity = c(rep(0, as.numeric(severity_dist[1,1])), 
                 rep(1, severity_dist[2,1]),
                 rep(2, severity_dist[3,1])),
    # 0 = Low, 1 = Moderate, 2 = High Severity
    
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
    
    pnc_enc = purrr::map(severity, ~sample_pnc(pnc_prob, .x))
    
  )
  
  # Final random state
  finish_seed <- list(.Random.seed)
  
  ## Store the random states
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
# sample_revised_preg() and replaces its values
# with those that are compatible with the other vectors.
process_out <- function(out) {
  # Create an empty vector of characters with the same length as 'probabilities'
  results <- vector("character", length(out))
  
  # Set the first 18 elements to "NA"
  results[1:18] <- NA
  
  # For the remaining elements, map 1 to "fetaldeath" and 0 to "delivery"
  results[19:length(out)] <- ifelse(out[19:length(out)] == 1, 
                                    "fetaldeath_next", 
                                    "livebirth_next")
  
  return(results)
}

# Example use:
# probabilities <- c(rep(0, 10), rep(1, 10), rep(0, 5), rep(1, 5))
# output <- process_probabilities(probabilities)
# print(output)



### FUNCTION sample_revised_preg()
# This function will provide revised pregnancy
# outcomes for the cases in which a pregnancy
# ends in preeclampsia.
sample_revised_preg <- function(data, sev){
  
  # Subset to the severity level of interest 
  # Get the probabilities for each gestational week
  probabilities <- subset(data, severity == sev)$p_fetaldeath
  
  # Use these probabilities with a bernoulli/binomial rv
  out <- rbinom(n = length(probabilities), 
                size=1,
                prob = probabilities)
  rev_out <- process_out(out)
  
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

















