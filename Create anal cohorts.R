#####################################################
# Program: Generate data.R
# Programmer: Chase
#
# Purpose: The program generates the final analytic 
# cohort defined by different values to create 
# missing outcomes (i.e., LTFU).
#####################################################



#####################################################
# Pull in the necessary libraries and functions
#####################################################

library(tidyverse)

# Pull in the cohort creation functions - These functions actually create the analytic cohorts from the raw data
source("create cohort functions.R")



generate_cohort <- function(b0_sev, beta12, gamma0, gamma1, save_name_cohort){
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  p_miss_outcome <- c(log(gamma0 / (1-gamma0)), log(gamma1/(1-gamma1)))
  
  # Create p_sev_beta - for logistic regression to determine LTFU due to severity
  p_sev_beta = c(b0_sev, log(beta12 / (1-beta12)), log((2*beta12 / (1 - (2*beta12)))))
  
}









#####################################################
# Function can be used to generate each scenario with 
# the specified inputs.
#####################################################

generate_scenario <- function(param_file, b0_sev, beta12, save_name_gen, 
                              gamma0, gamma1, save_name_cohort,
                              rr_abortion, rr_preec){
  
  ## Create the parameter list
  params_list_gen <- list(
    n_sim = 1,
    n = 1200,
    p_trt_sev = c(0.35, 0.50, 0.65),
    p_indx_pnc = c(0.25, 0.375, 0.375)
  )
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  p_miss_outcome <- c(log(gamma0 / (1-gamma0)), log(gamma1/(1-gamma1)))
  
  p_sev_beta = c(b0_sev, log(beta12 / (1-beta12)), log((2*beta12 / (1 - (2*beta12)))))
  
  ## Upload the parameter Excel files.
  potential_preg_untrt <- read_xlsx(param_file, sheet = "potential_preg_untrt")
  potential_preg_trt <- read_xlsx(param_file, sheet = "potential_preg_trt")
  potential_preec_untrt <- read_xlsx(param_file, sheet = "potential_preec_untrt")
  potential_preec_trt <- read_xlsx(param_file, sheet = "potential_preec_trt")
  revised_preg <- read_xlsx(param_file, sheet = "postpreec_preg")
  pnc_prob <- read_xlsx(param_file, sheet = "pnc_prob")
  
  ## Generate all the potential outcomes
  
  all_outcomes <- do.call(generate, params_list_gen) %>% 
    # Save the data generation values
    mutate(rr_trt_abortion = rr_abortion,
           rr_trt_preec = rr_preec)
  
  # Save the RDS file -- Potentially important for bootstrapping SEs
  saveRDS(all_outcomes, save_name_gen)
  
  ## Select the observed cohort
  all_outcomes2 <- create_cohort(all_outcomes, p_sev_beta, p_miss_outcome) %>% 
    mutate(diff_beta1_beta2 = beta12,
           gamma_0 = gamma0,
           gamma_1 = gamma1)
  
  # Save the RDS file
  saveRDS(all_outcomes2, save_name_cohort)
  
}







##################################################
# RR, Trt-Abortion = 0.8
# RR, Trt-Preeclampsia = 0.8
#################################################

# The Excel file with the simulation parameters will be the same
# for all analytic cohorts created through this mechanism

param_file <- "Parameters_Abortion0_8_Preeclampsia0_8.xlsx"
rr_abortion <- 0.8
rr_preec <- 0.8


###### 

