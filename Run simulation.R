#####################################################
# Program: Run simulation.R
# Programmer: Chase
# Date created: 06.27.2024
#
# Purpose: This program contains the code to run the 
# simulation. This will eventually be deprecated
# with the targets package. But, using for checking.
#####################################################


#####################################################
# Pull in the necessary libraries and functions
#####################################################
library(tidyverse)
library(readxl)
library(survival)
library(rlang)

# Pull in the data generation functions - These functions create the raw data
source("R/Data generation functions.R")

# Pull in the cohort creation functions - These functions actually create the analytic cohorts from the raw data
source("R/create cohort functions.R")





#####################################################
# Generate the data
#####################################################
n_sim <- 1
n <- 50000
rr_abortion <- 0.8
rr_preec <- 0.8
param_file <- "Parameters_Abortion0_8_Preeclampsia0_8.xlsx"


all_outcomes <- generate_dgm(n_sim, n, param_file, rr_abortion, rr_preec)

# How to do for multiple
# all_outcomes <- purrr::map_dfr(1:n_sim, ~generate_dgm(.x, n, param_file, rr_abortion, rr_preec))

saveRDS(all_outcomes, "data/DGM - Abortion08Preeclampsia08.rds")


#####################################################
# Create the cohort data
# all_outcomes <- readRDS("data/DGM - Abortion08Preeclampsia08.rds")
#####################################################

#### Missing: Beta1 = 0.05, Gamma0 = 0.1, Gamma1 = 0.02
# 
marginal_p_miss_severity <- 0.1
beta12 <- 0.7
  #c(0.7519006, 1.006749) # Severity less predictive
# beta12 <- c(1.222581, 1.742095) # Severity more predictive
marginal_p_miss_miscarriage <- 0.5
gamma1 <- 0.01
save_name_cohort <- "data/ab08preec08_beta001_gamma01_001.rds"

set.seed(2094857309) # Same seed throughout
test <- generate_cohort(all_outcomes, marginal_p_miss_severity, beta12, marginal_p_miss_miscarriage, gamma1)

saveRDS(test, save_name_cohort)





##################################################
# Run ANALYSES
#################################################

#test2 <- run_analysis(test)








