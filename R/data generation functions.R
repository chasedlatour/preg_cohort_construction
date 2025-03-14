#####################################################
# Program: Data generation functions. R
# Programmer: Chase
# Date last modified: 04.18.2024
#
# Purpose: The program contains all the R functions
# for generating the data.
#####################################################






#####################################################
# Function can be used to generate each scenario with 
# the specified inputs.
#####################################################

generate_dgm <- function(n_sim, n, param_file,
                         rr_abortion, rr_preec){
  
  
  ## Upload the parameter Excel files.
  potential_preg_untrt <- read_xlsx(param_file, sheet = "potential_preg_untrt")
  potential_preg_trt <- read_xlsx(param_file, sheet = "potential_preg_trt")
  potential_preec_untrt <- read_xlsx(param_file, sheet = "potential_preec_untrt")
  potential_preec_trt <- read_xlsx(param_file, sheet = "potential_preec_trt")
  revised_preg <- read_xlsx(param_file, sheet = "postpreec_preg")
  pnc_prob <- read_xlsx(param_file, sheet = "pnc_prob")
  
  ## Create the parameter list that are the same across simulations
  params_list_gen <- list(
    n_sim = n_sim,
    n = n, #5000,
    p_sev_dist = c(1/3, 1/3, 1/3), # Even distribution
    # p_trt_sev = c(0.4, 0.50, 0.6),
    p_trt_sev = c(0.25, 0.5, 0.75),
    p_indx_pnc = c(0.23, 0.33, 0.44),
    rrmage_abor = 1.5,
    ormage_preec = 2,
    potential_preg_trt = potential_preg_trt,
    potential_preg_untrt = potential_preg_untrt,
    potential_preec_untrt = potential_preec_untrt,
    potential_preec_trt = potential_preec_trt,
    revised_preg = revised_preg,
    pnc_prob = pnc_prob
  )
  
  ## Generate all the potential outcomes
  
  all_outcomes <- do.call(generate, params_list_gen) %>% 
    # Save the data generation values
    mutate(rr_trt_abortion = rr_abortion,
           rr_trt_preec = rr_preec)
  
  return(all_outcomes)
  
}















###################################################################################
#### FUNCTION: generate()
# This is the baseline function that will be used to generate the potential outcomes
# -- for each cohort.

# Inputs: (most saved in params_list_gen in the step before.)
# n_sim = simulation number (id)
# n = number of individuals in a cohort/repetition
# p_sev_dist = vector of patient distribution by disease severity
# p_trt_sev = vector of treatment probabilities by disease severity
# p_indx_pnc = vector of probabilities for week of indexing prenatal care visit (6, 9, 18 from LMP)
# potential_preg_untrt = dataset with the weekly probabilities for untreated potential preg outcomes
# potential_preg_trt = same but treated potential pregnancy outcomes
# potential_preec_untrt = dataset with weekly probabilities for untreated potential preeclampsia outcomes
# potential_preec_trt = same but treated potential preeclampsia outcomes
# revised_preg = dataset with the weekly probabilities for pregnancy outcomes after preeclampsia
# pnc_prob = dataset with the weekly probabilities for a prenatal encounter
###################################################################################

generate <- function(n_sim, n, p_sev_dist, p_trt_sev, p_indx_pnc,
                     rrmage_abor, ormage_preec,
                     potential_preg_untrt, potential_preg_trt, 
                     potential_preec_untrt, potential_preec_trt,
                     revised_preg, pnc_prob){
  
  # Get the severity distribution from a multinomial random variable
  # We assume that there is an equal distribution across three levels. 
  # Otherwise, this needs to be modified.
  severity_dist = rmultinom(n=1, size=n, prob=p_sev_dist)
  
  ######## GENERATE PEOPLE AND BASELINE VALUES
  
  # Create the dataset that going to output
  data <- dplyr::tibble(
    
    # Simulation ID -- Only important if we simulate more than one simulation cohort. Retaining in case we do
    sim_id = n_sim, 
    
    # Each individual's ID
    id = 1:n, 
    
    # Maternal age less than 35 or ge 35 - ADDED 2.28.25
    mage35 = rbinom(n=n, size=1, prob = 0.3),
    
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
    # preg_outcomes_untrt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_untrt, .x)),
    
    preg_outcomes_untrt = purrr::map2(
      severity, mage35,
      ~ {
        # Subset to the severity level of interest
        # prob_data <- potential_preg_untrt[potential_preg_untrt$severity == sev,] %>% 
        #   select(p_fetaldeath_next, p_livebirth_next, p_contpreg_next)
        prob_data <- potential_preg_untrt[potential_preg_untrt$severity == .x,]
        
        # Adjust probabilities if mage35 == 1
        prob_data$p_fetaldeath_next <- ifelse(.y == 1 & prob_data$gestweek_conception <= 16,
                                              prob_data$p_fetaldeath_next * rrmage_abor,
                                              prob_data$p_fetaldeath_next)
        # Now adjust associated probabilities of continuing pregnancy
        prob_data$p_contpreg_next <- ifelse(.y == 1 & prob_data$gestweek_conception <= 16,
                                            1 - prob_data$p_fetaldeath_next - prob_data$p_livebirth_next,
                                            prob_data$p_contpreg_next)
        
        prob_data2 <- prob_data[, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]
        
        # Define the potential options to select from
        options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")
        
        # Apply the sampling function to each row of probabilities
        outcomes <- apply(
          prob_data2, 1, function(prob) {
            sample(options, size = 1, prob = prob)
          }
        )
        
        # Return the outcomes as a list
        list(outcomes)
      }
    ),
    
    ## Step 3: Treated
    # preg_outcomes_trt = purrr::map(severity, ~sample_outcomes_for_id(potential_preg_trt, .x)),
    
    preg_outcomes_trt = purrr::map2(
      severity, mage35,
      ~ {
        # Subset to the severity level of interest
        prob_data <- potential_preg_trt[potential_preg_trt$severity == .x,]
        
        # Adjust probabilities if mage35 == 1
        prob_data$p_fetaldeath_next <- ifelse(.y == 1 & prob_data$gestweek_conception <= 16,
                                              prob_data$p_fetaldeath_next * rrmage_abor,
                                              prob_data$p_fetaldeath_next)
        # Now adjust associated probabilities of continuing pregnancy
        prob_data$p_contpreg_next <- ifelse(.y == 1 & prob_data$gestweek_conception <= 16,
                                            1 - prob_data$p_fetaldeath_next - prob_data$p_livebirth_next,
                                            prob_data$p_contpreg_next)
        
        prob_data2 <- prob_data[, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]
        
        # Define the potential options to select from
        options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")
        
        # Apply the sampling function to each row of probabilities
        outcomes <- apply(
          prob_data2, 1, function(prob) {
            sample(options, size = 1, prob = prob)
          }
        )
        
        # Return the outcomes as a list
        list(outcomes)
      }
    ),
    
    ##### GENERATE POTENTIAL PREECLAMPSIA OUTCOMES
    
    ## Step 4: Untreated
    # preec_outcomes_untrt = purrr::map(severity, ~sample_preeclampsia(potential_preec_untrt, .x)),
    
    preec_outcomes_untrt = purrr::map2(
      severity, mage35,
      ~{
        probabilities <- potential_preec_untrt[potential_preec_untrt$severity == .x,]$p_preeclampsia
        
        # Add the OR for mage: convert to log odds, add log-transformed OR then exponentiate
        if(.y == 1){
          probabiltiies <- exp(log(probabilities/(1-probabilities)) + log(ormage_preec))
        }
        
        list(rbinom(n = length(probabilities), size = 1, prob = probabilities))
      }
    ),
    
    ## Step 5: Treated
    # preec_outcomes_trt = purrr::map(severity, ~sample_preeclampsia(potential_preec_trt, .x)),
    
    preec_outcomes_trt = purrr::map2(
      severity, mage35,
      ~{
        probabilities <- potential_preec_trt[potential_preec_trt$severity == .x,]$p_preeclampsia
        
        # Add the OR for mage: convert to log odds, add log-transformed OR then exponentiate
        if(.y == 1){
          probabiltiies <- exp(log(probabilities/(1-probabilities)) + log(ormage_preec))
        }
        
        list(rbinom(n = length(probabilities), size = 1, prob = probabilities))
      }
    ),
    
    ##### GENERATE REVISED OUTCOMES IF PREECLAMPSIA -- Step 6
    
    # revised_preg = purrr::map(severity, ~sample_revised_preg(revised_preg, .x)),
    
    revised_preg = purrr::map(
      severity,
      ~ {
        probabilities <- revised_preg[revised_preg$severity == .x,]$p_fetaldeath
        n_NA <- sum(probabilities == 0)
        out <- c(rep(0, n_NA), rbinom(n = length(probabilities) - n_NA,
                                      size = 1,
                                      prob = probabilities[probabilities > 0]))
        results <- vector("character", length(probabilities))
        results[1:n_NA] <- NA
        results[(n_NA + 1):length(probabilities)] <- 
          ifelse(out[(n_NA + 1):length(out)] == 1, "fetaldeath_next", "livebirth_next")
        
        list(results)
      }
    ),
    
    ##### GENERATE PNC ENCOUNTERS -- Step 7
    
    # List of indicator variables for PNC encounters
    # We incorporate indexing prenatal encounter later
    # pnc_enc = purrr::map(severity, ~sample_pnc(pnc_prob, .x))
    
    pnc_enc = purrr::map(
      severity,
      ~{
        probabilities <- pnc_prob[pnc_prob$severity == .x,]$p_pnc
        list(rbinom(n = length(probabilities), size = 1, prob = probabilities))
      }
    )
    
  )
  
  return(data)
  
}





###################################################################################
### FUNCTION: sample_outcomes_for_id()
# This function will create the potential pregnancy outcomes based upon 
# the probabilities recorded in the Excel file.
# This will output a list of the potential pregnancy outcomes for one person.
#
# Inputs:
# - data = dataset with the weekly probabilities for potential pregnancy outcomes
# - sev = severity value for that person
###################################################################################


# sample_outcomes_for_id <- function(data, sev) {
# 
# 
#   data2 <- subset(data, severity == sev)
#   # data2 <- data %>% filter(severity == sev)
# 
#   # Define the potential options to select from
#   options <- c("fetaldeath_next", "livebirth_next", "contpreg_next")
# 
#   # Extract the matrix of probabilities
#   probabilities <- data2[, c("p_fetaldeath_next", "p_livebirth_next", "p_contpreg_next")]
# 
#   # Apply the sampling function to each row of probabilities using apply
#   outcomes <- apply(probabilities, 1, function(prob) {
# 
#     sample(options, size = 1, prob = prob)
# 
#   })
# 
#   return(list(outcomes))
# }










###################################################################################
### FUNCTION sample_preeclampsia()
# This function will create potential preeclampsia outcomes 
# This will output a list of the potential preeclampsia
# outcomes for one person.
#
# INPUT:
# - data = dataset with the weekly potential preeclampsia probabilities
# - sev = severity value for that person
###################################################################################

# sample_preeclampsia <- function(data, sev) {
#   
#   # Subset to the severity level of interest 
#   # Get the probabilities for each gestational week
#   probabilities <- subset(data, severity == sev)$p_preeclampsia
#   # probabilities <- (data %>% filter(severity == sev))$p_preeclampsia
#   
#   # Use these probabilities with a bernoulli/binomial rv
#   out <- rbinom(n = length(probabilities), 
#                 size=1,
#                 prob = probabilities)
#   
#   return(list(out))
# }






# ###################################################################################
# #### FUNCTION: process_out()
# # This function takes in the vector out from 
# # sample_revised_preg() -- below -- and replaces its values
# # with those that are compatible with the other vectors.
# # Further, those where the probability of a revised pregnancy outcome was 0
# # are replaced with NA.
# #
# # INPUT:
# # - out = vector from sample_revised_preg
# ###################################################################################
# 
# process_out <- function(out, n_NA) {
#   
#   # Create an empty vector of characters with the same length as 'probabilities'
#   results <- vector("character", length(out))
#   
#   # Set the first 17 elements to "NA"
#   results[1:n_NA] <- NA
#   
#   # For the remaining elements, map 1 to "fetaldeath" and 0 to "delivery"
#   results[(n_NA+1):length(out)] <- ifelse(out[(n_NA+1):length(out)] == 1,
#                                           "fetaldeath_next", 
#                                           "livebirth_next")
#   
#   return(results)
# }
# 
# # Example use:
# # probabilities <- c(rep(0, 10), rep(1, 10), rep(0, 5), rep(1, 5))
# # output <- process_out(probabilities)
# # print(output)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###################################################################################
# ### FUNCTION sample_revised_preg()
# # This function will provide revised pregnancy
# # outcomes for the cases in which a pregnancy
# # ends in preeclampsia.
# ###################################################################################
# 
# sample_revised_preg <- function(data, sev){
#   
#   # Subset to the severity level of interest 
#   # Get the probabilities for each gestational week
#   probabilities <- subset(data, severity == sev)$p_fetaldeath
#   # probabilities <- (data %>% filter(severity == sev))$p_fetaldeath
#   
#   # Determine how many NAs there should be
#   n_NA <- length(probabilities[probabilities == 0])
#   
#   # Use these probabilities with a bernoulli/binomial rv
#   out <- c(rep(0, n_NA), rbinom(n = length(probabilities) - n_NA, 
#                 size=1,
#                 prob = probabilities))
#   # out is a numeric vector - we want a character vector with NA
#   # for those values where prob was 0
#   
#   # Revised here using above step.
#   rev_out <- process_out(out, n_NA)
#   
#   return(list(rev_out))
#   
# }





# 
# 
# ### FUNCTION: sample_pnc()
# # This function will provide prenatal care encounter
# # indicators.
# sample_pnc <- function(data, sev) {
#   
#   # Subset to the severity level of interest 
#   # Get the probabilities for each gestational week
#   probabilities <- subset(data, severity == sev)$p_pnc
#   # probabilities <- (data %>% filter(severity == sev))$p_pnc
#   
#   # Use these probabilities with a bernoulli/binomial rv
#   out <- rbinom(n = length(probabilities), 
#                 size=1,
#                 prob = probabilities)
#   
#   return(list(out))
# }



