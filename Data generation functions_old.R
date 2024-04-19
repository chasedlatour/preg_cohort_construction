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



# This is the baseline function that will be 
# -- used to generate each cohort.
# This function will be run to generate a cohort 
# -- for each simulation.
each_sim <- function(n_sim, n){
  
  initial_seed <- list(.Random.seed)
  
  # Set the seed for the simulation as equal to n_sim.
  #set.seed(n_sim)
  # Decided that this should only be done once at the beginning of the simulation, per Morris et al. 2019
  
  # Generate 0 through 40 gestational weeks - 1 vector
  gw = list(seq(0,40, by = 1))
  
  # Create the dataset that going to output
  data <- dplyr::tibble(
    
    sim_id = n_sim, # Simulation ID
    
    id = 1:n, # Each individual's ID
    
    # Assign treatment here - doesn't matter when but efficient here
    trt = rbinom(n=n, size=1, prob=0.5)
    
  ) %>% 
    rowwise %>% 
    mutate(
      
      #### PHASE 1
      
      # Create the pregnancy outcomes for Phase 1 of the data generation
      #phase1_outcomes = list(sample_outcomes_for_id(phase1)), 
      phase1_outcomes = sample_outcomes_for_id(phase1), 
      
      # Identify the last gestational week that a person experiences of pregnancy up to 20 weeks of gestation
      # -- from LMP or 18 weeks from conception.
      last_GA_pre18 = pre18(phase1_outcomes),
      
      # Create the index outcomes for the trial
      index_gw = list(rbinom(nrow(phase1), size=1, prob=phase1$p_index)),
      
      # Identify the first index event -- At least one should be 1 because
      # -- Probability of index at week 20 is 1.
      first_index = which(index_gw == 1)[1] - 1,
      # Subtract 1 because we index at conception or week 0, and R indexes at 1
      # The same was done for index_gw last_GA_pre20 in the function.
      
      # Create an indicator variable as to whether the patient actually
      # -- enters the trial. If last_GA_pre20 < first_index, then they can't enter the trial.
      trial_participant = (last_GA_pre18 >= 2) && (first_index <= last_GA_pre18),
      # Use 2 here instead of 4 because from conception and not LMP.
      
      #### PHASE 2
      
      # Not going to subset to trial participants for now
      # Generate pregnancy outcomes for phase2, focusing on those pregnancies that are treated
      # phase2_outcomes = ifelse(trt == 0, # If untreated
      #                          #I(list(phase1_outcomes)),
      #                          list(phase1_outcomes),
      #                          resample_outcomes(first_index, phase1_outcomes, phase2)),
      
      # Create all of the outcomes if they had been treated at the first_index
      phase2_outcomes = resample_outcomes(first_index, phase1_outcomes, phase2),
      
      #### PHASE 3 -- determine whether a person developed preeclampsia & revised outcomes
      
      # Calculate their odds of preeclampsia from logistic regression
      odds_preeclampsia = list(exp(phase3$ln_odds_preeclampsia + (phase3$treat_effect_ln_OR * trt))),
      # Calculate the probability of preeclampsia
      prob_preeclampsia = list(c(rep(0,25), odds_preeclampsia/(1+odds_preeclampsia))),
      
      # Use binomial r.v. to determine if someone has preeclampsia
      preeclampsia_list = list(rbinom(length(prob_preeclampsia), size=1, prob=prob_preeclampsia)),
      
      # Determine the final pregnancy outcomes with preeclampsia timing 
      ## These are outcomes if all were treated (i.e., potential outcomes)
      final_pregnancy_outcomes_trt = list(preg_outcome(phase2_outcomes, preeclampsia_list, phase3)),
      ## These are outcomes if all were untreated (i.e., potential outcomes)
      final_pregnancy_outcomes_untrt = list(preg_outcome(phase1_outcomes, preeclampsia_list, phase3)),
      
      
      ##### PHASE 4 -- Determine if people have SGA
      ## These are outcomes if all were treated (i.e., potential outcomes)
      sga_trt = sga_func(final_pregnancy_outcomes_trt, 1, phase4),
      ## These are outcomes if all were untreated (i.e., potential outcomes)
      sga_untrt = sga_func(final_pregnancy_outcomes_untrt, 0, phase4),
      #sga = sga_func(final_pregnancy_outcomes, trt, phase4)
      
      
      ##### PHASE 5 -- Censoring -- STILL NEED TO INCORPORATE THIS.
      first_censoring_gw = censoring(phase5, first_index)
      
      
      
    ) %>% 
    ungroup() %>% 
    mutate(final_pregnancy_outcomes_trt = map(final_pregnancy_outcomes_trt, 
                                              ~set_names(.x, str_c(names(.x), "_trt"))),
           final_pregnancy_outcomes_untrt = map(final_pregnancy_outcomes_untrt, 
                                                ~set_names(.x, str_c(names(.x), "_untrt")))) %>% 
    unnest_wider(c(final_pregnancy_outcomes_trt, final_pregnancy_outcomes_untrt)) 
  
  finish_seed <- list(.Random.seed)
  
  ## Store the random states
  data <- data %>% 
    mutate(
      start_seed = initial_seed,
      end_seed = finish_seed
    )
  
  # #Testing
  # testtrt <- subset(data, sga_trt == 1)
  # table(testtrt$preg_outcome_final_trt)
  # testutrt <- subset(data, sga_untrt == 1)
  # table(testutrt$preg_outcome_final_untrt)
  
  return(data)
  
}



### Below, are all of the functions that will be used within the `each_sim()` function.



# This function will create the pregnancy outcomes based upon the probabilities
# -- recorded in phase1, phase2, phase3, and phase4.
sample_outcomes_for_id <- function(data) {
  #id_data <- subset(your_data, id == id)  # Filter data for the specific id
  
  # Apply this function to each row of phase1-phase4 data where needed.
  # -- Create vector of 1:nrow(data) and then apply the function below
  outcomes <- sapply(1:nrow(data), function(i) { 
    
    # Indicate the potential options to select from
    options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")
    
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



