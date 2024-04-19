#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 01.22.2024
#
# Purpose: The program contains all the R functions
# for generating the data. The goal of doing this
# is to ensure that these functions only need to be
# editing one time and will apply throughout all
# simulations easily.
#####################################################


# Libraries that you need for these functions
library(tidyverse)
library(readxl)



# This is the baseline function that will be 
# -- used to generate each cohort.
# This function will be run to generate a cohort 
# -- for each simulation.

# Testing:
n_sim <- 1
n <- 150

## Upload the parameter Excel files.
potential_preg_untrt <- read_xlsx("Simulation Parameters.xlsx", 
                                  sheet = "potential_preg_untrt")
potential_preg_trt <- read_xlsx("Simulation Parameters.xlsx", 
                                sheet = "potential_preg_trt")
# phase2 <- read_xlsx("Simulation Parameters.xlsx", sheet = "Phase2")
# phase3 <- read_xlsx("Simulation Parameters.xlsx", sheet = "Phase3")
# phase4 <- read_xlsx("Simulation Parameters.xlsx", sheet = "Phase4")
# phase5 <- read_xlsx("Simulation Parameters.xlsx", sheet = "Phase5")

#each_sim <- function(n_sim, n){
  
  # Per Morris et al. 2019, we save the random seed at the beginning of the simulation
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
                 rep(2, severity_dist[3,1]))
    # 0 = Low, 1 = Moderate, 2 = High Severity
    
    
  ) %>% 
    rowwise %>% 
    mutate(
      
      #### STEP 1: POTENTIAL PREGNANCY OUTCOMES
      # Generate each patient's potential pregnancy outcomes at each 
      # gestational week.
      
      preg_outcomes_untrt = sample_outcomes_for_id(potential_preg_untrt, severity),
      
      preg_outcomes_trt = sample_outcomes_for_id(potential_preg_trt, severity)
      ) #, 
      


# This function will create the pregnancy outcomes based upon the probabilities
# -- recorded in phase1, phase2, phase3, and phase4.
sample_outcomes_for_id <- function(data, sev) {
  
  # Subset to the severity level of interest
  data2 <- subset(data, severity == sev)
  # data2 <- subset(potential_preg_untrt, severity == 1)

  # Apply this function to each row of phase1-phase4 data where needed.
  # -- Create vector of 1:nrow(data) and then apply the function below
  outcomes <- sapply(1:nrow(data2), function(i) {

    # Indicate the potential options to select from
    options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")

    # Assign the probabilities to each of the potential pregnancy outcome options based upon
    # -- the corresponding rows in the phase1-4 files.
    probabilities <- data2[i, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]

    # Sample an option based on probabilities
    sampled_option <- sample(options, size = 1, prob = probabilities)

    # Return that sampled option. This will be stored in the vector outcomes
    return(sampled_option)
  })

  #list(outcomes = outcomes)
  return(list(outcomes))
}

#Testing: outcomes <- test$phase1_outcomes[[1]]


# This function will create a vector of index events based upon the probabilities
# -- recorded in phase1.
sample_index_for_id <- function(data) {
  #id_data <- subset(your_data, id == id)  # Filter data for the specific id
  
  # Apply this function to each row of phase1-phase4 data where needed.
  # -- Create vector of 1:nrow(data) and then apply the function below
  indices <- sapply(1:nrow(data), function(i) { 
    
    # Determine the probability of an index at that gestational week
    prob = phase1$p_index[i]
    
    # Sample the index value based upon that probability
    
    # Indicate the potential options to select from
    options <- c("contpreg_next")
    
    # Assign the probabilities to each of the potential pregnancy outcome options based upon 
    # -- the corresponding rows in the phase1-4 files.
    probabilities <- data[i, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]
    
    # Sample an option based on probabilities
    sampled_option <- sample(options, size = 1, prob = probabilities)
    
    # Return that sampled option. This will be stored in the vector outcomes
    return(sampled_option)
  })
  
  #list(outcomes = outcomes)
  return(list(outcomes))
}






pre18 <- function(outcomes){
  # Create a subset of gestational weeks 0 through 18, though note that these
  # -- are indexed as 1:20 because of R's coding style (vector indices start at 1). 
  # -- Notably for moving forward, outcomes are observed at the following gestational week
  outcomes_sub <- outcomes[1:19]  # Represents 0:18
  
  # Identify the first gestational week where we will see a fetal death on the 
  # -- following gestational week.
  first_ga <- which(outcomes_sub == 'fetaldeath_next')[1]
  
  # Identify the last gestational week that someone could be indexed into the cohort.
  # -- This would be the last week before someone experiences a fetal death.
  last_ga <- ifelse(is.na(first_ga),
                    18,
                    first_ga - 1) 
  # Subtracted 1 because first_ga is indexed starting at 1, while our gestational weeks are indexed
  # -- starting at 0.
  
  return(last_ga)
}


# Re-sample pregnancy outcomes for those pregnancies that were treated.
# -- Want to retain their phase1 outcomes until their index event, and then 
# -- resample all pregnancy outcomes starting at their index.
# -- This assumes an immediate treatment effect.
resample_outcomes <- function(first_index, phase1_outcomes, phase2){
  
  #browser()
  
  # Add 1 to accomodate R vector indexing starting at 1
  index <- first_index + 1 
  
  # Grab all the phase1 outcomes that don't change.
  # -- Want all those values BEFORE the index
  phase1 <- phase1_outcomes[1:(index-1)]
  
  # Re-sample outcomes for phase2
  phase2 <- unlist(sample_outcomes_for_id(phase2))
  
  # Mark suggestion: Rename phase2 here.
  
  # Now merge them together
  outcomes <- c(phase1, phase2[index:length(phase2)])
  
  # Coerce outcomes to a vector
  # -- This step is necessary to make sure that the values come out as the same type as phase1_outcomes.
  outcomes <- unlist(outcomes)
  
  list(outcomes = outcomes)
  
}


## This function generates the final pregnancy outcome for each of the observed pregnancies

preg_outcome <- function(phase2_outcomes, preeclampsia_list, phase3){
  #browser()
  
  # Determine which occurred first: pregnancy outcome or preeclampsia.
  # -- If they occurred within the same gestational week, we will assume that the preeclampsia occurred first.
  
  ## First outcome from phase2_outcomes
  first_outcome_index <- which(phase2_outcomes %in% c('fetaldeath_next','livebirth_next'))[1]
  
  ## First preeclampsia from preeclampsia_list
  first_preeclampsia <- which(preeclampsia_list == 1)[1]
  
  # If no preeclampsia evidence or first outcome from phase2_outcomes, 
  # -- then return a list of pregnancy outcomes from phase2_outcomes
  # -- Do not worry about +/- 1 here because they're both indexed in a vector from 1.
  if(first_outcome_index < first_preeclampsia | is.na(first_preeclampsia)){
    return(list(
      
      # Determine what the final outcome was from phase2_outcomes
      preg_outcome_final = sub("_next$", "", phase2_outcomes[first_outcome_index]),
      
      # Determine the gestational timing of the first outcome
      # -- Originally, we added 1 here because we are assuming that the outcome is observed at the 
      # -- week after they're determined.
      # -- However, R indexes a vector from 1, not 0, and we are indexing from 0, so we would
      # -- subtract 1 to account for that. -- They equal out.
      preg_outcome_final_gw = first_outcome_index,
      
      # Make preeclampsia indicator - Preeclampsia didn't occur before the outcome and so not observed
      preeclampsia = 0
    ))
  }
  
  
  ## If first or equally first outcome is preeclampsia, then regenerate the pregnancy outcome and timing
  if(first_outcome_index >= first_preeclampsia){
    
    # Regenerate the pregnancy outcome
    p_fetaldeath <- subset(phase3, gestweek_conception == (first_preeclampsia - 1))$p_fetaldeath
    p_livebirth <- subset(phase3, gestweek_conception == (first_preeclampsia - 1))$p_livebirth
    options <- c("fetaldeath", "livebirth")
    # Assign the corresponding probabilities to each of the potential pregnancy outcomes
    probabilities <- c(p_fetaldeath, p_livebirth)
    outcome <- sample(options, size=1, prob=probabilities)
    
    return(list(
      
      # Determine the final pregnancy outcome from the probabilities in phase 3
      preg_outcome_final = outcome,
      
      # Determine the gestational timing of the outcome - same logic on not adding/subtracting 1
      preg_outcome_final_gw = first_preeclampsia,
      
      # Make preeclampsia indicator - Preeclampsia observed before the outcome so observed
      preeclampsia = 1
      
    ))
  }
  
}

sga_func <- function(final_pregnancy_outcomes, trt, phase4){
  
  # Determine probability of SGA based upon trt and preeclampsia
  prob_sga <- subset(phase4, trt_value == trt & preeclampsia_flag ==
                       final_pregnancy_outcomes$preeclampsia)$p_sga
  
  # Determine if they actually had SGA based upon this probability
  sga <- ifelse(final_pregnancy_outcomes$preg_outcome_final == "fetaldeath",
                0,
                rbinom(1,1, prob = prob_sga))
  #sga <- rbinom(1,1, prob = prob_sga)
  
  # SGA does not apply if the pregnancy did not end in a live birth
  # sga2 <- ifelse(final_pregnancy_outcomes$preg_outcome_final == "fetaldeath",
  #                0,
  #                sga)
  
  return(sga)
}


# This function will determine censoring events based upon the probabilities 
# -- recorded in phase5.
censoring <- function(data, first_index){
  
  # Create a vector of the censoring probabilities
  # -- Each element is a gestational week
  prob_censoring <- data$p_censoring
  
  # Determine the length of prob_censoring for easy reference later
  n <- length(prob_censoring)
  
  # Create a vector with the realized censoring events based
  # -- upon the prob_censoring vecotr
  censor <- rbinom(n, size = 1, prob=prob_censoring)
  
  # Determine the first_index - we are interested in the first censoring
  # -- event after this index.
  first_index2 <- first_index + 1
  # We add one because first_index gw starts at 0, but R vectors start at 1.
  
  # Determine the gw of the first censoring event
  censoring <- which(censor == 1)
  #first_censor <- which(censor[first_index2:n]==1)[1] - 1
  # Subtract one to get back to gw indexed at 0 or conception
  
  first_censor_index <- which(censoring > first_index2)
  
  first_censor <- ifelse(all(censor == 0),
                         NA,
                         censoring[first_censor_index] - 1)
  
  # Return the gestational week of the first censoring event
  # -- Note that this is independent of outcome generation and exposure
  # -- Thus, uninformative.
  return(first_censor)
  
}



