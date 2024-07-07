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
source("Data generation functions.R")

# Pull in the cohort creation functions - These functions actually create the analytic cohorts from the raw data
source("create cohort functions.R")





#####################################################
# Generate the data
#####################################################

n_sim <- 1
n <- 50000

##################################################
# RR, Trt-Abortion = 0.8
# RR, Trt-Preeclampsia = 0.8
#################################################

# The Excel file with the simulation parameters will be the same
# for all analytic cohorts created through this mechanism

param_file <- "Parameters_Abortion0_8_Preeclampsia0_8.xlsx"
rr_abortion <- 0.8
rr_preec <- 0.8
save_name_gen <- "DGM - Abortion08Preeclampsia08.rds"

# All RDS files are saved
set.seed(2094857309)
all_outcomes <- generate_dgm(n_sim, n, param_file, rr_abortion, rr_preec)

# For multiple, optimize this somehow
#all_outcomes <- map_dfr(1:n_sim, ~generate_dgm(.x, n, param_file, rr_abortion, rr_preec))

# Save the RDS file -- Potentially important for bootstrapping SEs
saveRDS(all_outcomes, save_name_gen)


###### 




#####################################################
# Create the analytic cohorts.
#####################################################




##################################################
# RR, Trt-Abortion = 0.8
# RR, Trt-Preeclampsia = 0.8
#################################################



#### Missing: Beta1 = 0.01, Gamma0 = 0.1, Gamma1 = 0.001

beta12 <- 0.01
gamma0 <- 0.1
gamma1 <- 0.001
save_name_cohort <- "ab08preec08_beta001_gamma01_001.rds"

set.seed(2094857309) # Same seed throughout
test <- generate_cohort(all_outcomes, beta12, gamma0, gamma1)#, save_name_cohort)

saveRDS(test, "ab08preec08_beta001_gamma01_001.rds")




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








