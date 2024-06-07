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





#####################################################
# Set values that are going to hold for ALL DGMs and 
# analytic samples
#####################################################
p_sev <- 1 - 0.9^(1/40) # Weekly prob (over 40 weeks) of being LTFU at lowest severity
b0_sev <- log(p_sev / (1-p_sev)) # beta 0 in the logistic regression for LTFU due to severity is NOT varied across scenarios







#####################################################
# Function can be used to generate each analytic 
# cohort.
#####################################################

generate_cohort <- function(data, b0_sev, beta12, gamma0, gamma1, save_name_cohort){
  
  ## Create p_miss_outcome - for logistic regression to determine LTFU due to missing outcome
  p_miss_outcome <- c(log(gamma0/(1-gamma0)), 
                      log(gamma1/(1-gamma1))
                      )
  
  # Create p_sev_beta - for logistic regression to determine LTFU due to severity
  p_sev_beta = c(b0_sev, log(beta12 / (1-beta12)), log((2*beta12 / (1 - (2*beta12)))))
  
  ## Select the observed cohort
  all_outcomes2 <- create_cohort(data, p_sev_beta, p_miss_outcome) %>% 
    mutate(diff_beta1_beta2 = beta12,
           gamma_0 = gamma0,
           gamma_1 = gamma1)
  
  # Save the RDS file
  saveRDS(all_outcomes2, save_name_cohort)
  
}












##################################################
# RR, Trt-Abortion = 0.8
# RR, Trt-Preeclampsia = 0.8

# Same DGM data for all of these
data <- readRDS('DGM - Abortion08Preeclampsia08.rds')
#################################################



#### Missing: Beta1 = 0.01, Gamma0 = 0.1, Gamma1 = 0.001

beta12 <- 0.01
gamma0 <- 0.1
gamma1 <- 0.001
save_name_cohort <- "ab08preec08_beta001_gamma01_001.rds"

set.seed(1234) # Same seed throughout
set.seed(5678)
generate_cohort(data, b0_sev, beta12, gamma0, gamma1, save_name_cohort)





#### Missing: Beta1 = 0.05, Gamma0 = 0.1, Gamma1 = 0.02

beta12 <- 0.05
gamma0 <- 0.1
gamma1 <- 0.02
save_name_cohort <- "ab08preec08_beta005_gamma01_002.rds"

set.seed(1234) # Same seed throughout
set.seed(5678)
generate_cohort(data, b0_sev, beta12, gamma0, gamma1, save_name_cohort)

