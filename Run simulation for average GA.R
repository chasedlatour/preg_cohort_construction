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

saveRDS(all_outcomes, "data/DGM - 5 million for expected GA.rds")

