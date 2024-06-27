#####################################################
# Program: Generate data.R
# Programmer: Chase
#
# Purpose: The program generates the data for each
# DGM defined by different RR values for the effect
# of treatment on abortion and preeclampsia.
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


 