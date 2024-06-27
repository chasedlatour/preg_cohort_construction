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
generate_cohort(data, b0_sev, beta12, gamma0, gamma1, save_name_cohort)





#### Missing: Beta1 = 0.05, Gamma0 = 0.1, Gamma1 = 0.02

beta12 <- 0.05
gamma0 <- 0.1
gamma1 <- 0.02
save_name_cohort <- "ab08preec08_beta005_gamma01_002.rds"

set.seed(1234) # Same seed throughout
set.seed(5678)
generate_cohort(data, b0_sev, beta12, gamma0, gamma1, save_name_cohort)

