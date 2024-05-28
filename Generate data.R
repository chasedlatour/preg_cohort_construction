#####################################################
# Program: Generate data.R
# Programmer: Chase
#
# Purpose: The program generates the data for each
# DGM defined by different RR values for the effect
# of treatment on abortion and preeclampsia.
#####################################################


#####################################################
# Pull in the necessary libraries and functions
#####################################################
library(tidyverse)
library(readxl)

# Pull in the data generation functions - These functions create the raw data
source("Data generation functions.R")







#####################################################
# Function can be used to generate each scenario with 
# the specified inputs.
#####################################################

generate_dgm <- function(param_file, save_name_gen,
                         rr_abortion, rr_preec){
  
  
  ## Upload the parameter Excel files.
  potential_preg_untrt <- read_xlsx(param_file, sheet = "potential_preg_untrt")
  potential_preg_trt <- read_xlsx(param_file, sheet = "potential_preg_trt")
  potential_preec_untrt <- read_xlsx(param_file, sheet = "potential_preec_untrt")
  potential_preec_trt <- read_xlsx(param_file, sheet = "potential_preec_trt")
  revised_preg <- read_xlsx(param_file, sheet = "postpreec_preg")
  pnc_prob <- read_xlsx(param_file, sheet = "pnc_prob")
  
  ## Create the parameter list
  params_list_gen <- list(
    n_sim = 1,
    n = 5000, #1200,
    p_trt_sev = c(0.35, 0.50, 0.65),
    p_indx_pnc = c(0.25, 0.375, 0.375),
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
  
  # Save the RDS file -- Potentially important for bootstrapping SEs
  saveRDS(all_outcomes, save_name_gen)
  
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
save_name_gen <- "DGM - Abortion08Preeclampsia08.rds"

# All RDS files are saved
set.seed(1234)
generate_dgm(param_file, save_name_gen, rr_abortion, rr_preec)


###### 

