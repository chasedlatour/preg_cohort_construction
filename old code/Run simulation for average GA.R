#####################################################
# Program: Run simulation for average GA.R
# Programmer: Chase
# Date created: 07.26.2024
#
# Purpose: This program generates a cohort of
# 5 million people who become pregnant so that 
# we can have a good idea of the expected gestational
# age among those individuals who have a pregnancy
# survive to their initial prenatal encounter.
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
n <- 5000000
rr_abortion <- 0.8
rr_preec <- 0.8
param_file <- "Parameters_Abortion0_8_Preeclampsia0_8.xlsx"

set.seed(873654)
all_outcomes <- generate_dgm(n_sim, n, param_file, rr_abortion, rr_preec)

# How to do for multiple
# all_outcomes <- purrr::map_dfr(1:n_sim, ~generate_dgm(.x, n, param_file, rr_abortion, rr_preec))

saveRDS(all_outcomes, "DGM - 5 million for expected GA.rds")



###################################################
# Now run the analysis function so that we can 
# get the distribution by GA
###################################################

set.seed(873654)

# Betas don't matter for this exercise - just need something 
marginal_p_miss_severity <- 0.1
beta12 <- 0.7
#c(0.7519006, 1.006749) # Severity less predictive
# beta12 <- c(1.222581, 1.742095) # Severity more predictive
marginal_p_miss_miscarriage <- 0.5
gamma1 <- 0.01
save_name_cohort <- "DGM - COHORT - 5 million for expected GA.rds"

set.seed(2094857309) # Same seed throughout
test <- generate_cohort(all_outcomes, marginal_p_miss_severity, beta12, marginal_p_miss_miscarriage, gamma1)

saveRDS(test, save_name_cohort)







###################################################
# Get the GA distribution of miscarriages
# among the untreated
###################################################

data <- readRDS("DGM - COHORT - 5 million for expected GA.rds")

data <- bind_rows(data)

# Subset to non-initiators
non_initiators <- subset(data, trt == 0)

miscarriages <- subset(non_initiators, pregout_t_pre_miss < 18) ## Approx 260k

#Average GA
miscarriages %>% 
  ungroup() %>% 
  summarize(avg_ga = mean(pregout_t_pre_miss))


expected <- as.double(subset(data, trt == 0 & pregout_t_pre_miss < 18) %>% 
                        ungroup() %>% 
                        