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


#### Missing: Beta1 = 0.05, Gamma0 = 0.1, Gamma1 = 0.02
# 
# beta12 <- 0.05
# gamma0 <- 0.1
# gamma1 <- 0.02
# save_name_cohort <- "ab08preec08_beta005_gamma01_002.rds"
# 
# set.seed(2094857309) # Same seed throughout
# generate_cohort(data, beta12, gamma0, gamma1, save_name_cohort)
# 




##################################################
# Run ANALYSES
#################################################

#test2 <- run_analysis(test)








