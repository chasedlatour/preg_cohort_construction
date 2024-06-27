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

# Pull in the data generation functions - These functions create the raw data
source("Data generation functions.R")

# Pull in the cohort creation functions - These functions actually create the analytic cohorts from the raw data
source("create cohort functions.R")





#####################################################
# Generate the data
#####################################################


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
set.seed(1234)
set.seed(5678)
generate_dgm(param_file, save_name_gen, rr_abortion, rr_preec)


###### 




#####################################################
# Create the analytic cohorts.
#####################################################




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
generate_cohort(data, beta12, gamma0, gamma1, save_name_cohort)





#### Missing: Beta1 = 0.05, Gamma0 = 0.1, Gamma1 = 0.02

beta12 <- 0.05
gamma0 <- 0.1
gamma1 <- 0.02
save_name_cohort <- "ab08preec08_beta005_gamma01_002.rds"

set.seed(1234) # Same seed throughout
set.seed(5678)
generate_cohort(data, beta12, gamma0, gamma1, save_name_cohort)






