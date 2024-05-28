########################################
# Program: DGM Scenario 1.R
# Programmer: Chase Latour
# Date last Modified: 01/22/2024
#
# Purpose: Generate the cohort for 
# scenario 1.
########################################

# Pull in the necessary libraries
library(tidyverse)
library(readxl)

# Pull in the data generation functions.
source("Data generation functions.R")
source("create cohort functions.R")



################################
# TESTING
################################

## Upload the parameter Excel files.
potential_preg_untrt <- read_xlsx("Simulation Parameters.xlsx", sheet = "potential_preg_untrt")
potential_preg_trt <- read_xlsx("Simulation Parameters.xlsx", sheet = "potential_preg_trt")
potential_preec_untrt <- read_xlsx("Simulation Parameters.xlsx", sheet = "potential_preec_untrt")
potential_preec_trt <- read_xlsx("Simulation Parameters.xlsx", sheet = "potential_preec_trt")
revised_preg <- read_xlsx("Simulation Parameters.xlsx", sheet = "postpreec_preg")
pnc_prob <- read_xlsx("Simulation Parameters.xlsx", sheet = "pnc_prob")

# Generate the Data
all_outcomes <- generate(n_sim = 1, n = 1200, p_sev_beta = c(-3, 0.1, 0.2))

# Identify key variables and then output a dataset with 3 columns - one for each trial.
all_outcomes2 <- identify_key_vars(all_outcomes)


### Look at how things look

# Week 4
wk4 <- all_outcomes2[[2]][[1]]

wk4_sum <- wk4 %>% 
  ungroup() %>% 
  summarize(n = n(),
            n_include = sum(include_4wk, na.rm = FALSE),
            n_preeclampsia0 = sum(preeclampsia0, na.rm = FALSE),
            final_preg0_t_min = min(final_preg0_t, na.rm = FALSE),
            final_preg0_t_med = median(final_preg0_t, na.rm = FALSE),
            final_preg0_t_max = max(final_preg0_t, na.rm = FALSE),
            n_preeclampsia1 = sum(preeclampsia1, na.rm = FALSE),
            final_preg1_t_min = min(final_preg0_t, na.rm = FALSE),
            final_preg1_t_med = median(final_preg0_t, na.rm = FALSE),
            final_preg1_t_max = max(final_preg0_t, na.rm = FALSE),
            missing_min = min(missing_sev_wk4, na.rm = TRUE),
            missing_med = median(missing_sev_wk4, na.rm = TRUE),
            missing_max = max(missing_sev_wk4, na.rm = TRUE),
            missing_n = sum(is.na(missing_sev_wk4))
            )

View(wk4_sum)

# Look at the distribution of pregnancy outcome timing
hist(wk4$final_preg0_t)

hist(wk4$final_preg1_t)




















################################
# SCENARIO 1
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 1.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 1.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 1.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 1.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 1.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario1.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)




################################
# SCENARIO 2
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 2.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 2.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 2.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 2.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 2.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario2.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)





################################
# SCENARIO 3
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 3.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 3.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 3.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 3.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 3.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario3.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)





################################
# SCENARIO 4
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 4.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 4.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 4.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 4.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 4.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario4.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)




################################
# SCENARIO 5
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 5.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 5.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 5.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 5.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 5.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario5.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)




################################
# SCENARIO 6
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 6.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 6.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 6.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 6.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 6.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario6.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)





################################
# SCENARIO 7
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 7.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 7.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 7.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 7.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 7.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario7.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)





################################
# SCENARIO 8
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 8.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 8.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 8.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 8.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 8.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario8.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)





################################
# SCENARIO 9
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 9.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 9.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 9.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 9.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 9.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario9.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)






################################
# SCENARIO 10
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 10.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 10.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 10.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 10.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 10.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario10.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)








################################
# SCENARIO 11
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 11.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 11.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 11.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 11.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 11.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario11.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)








################################
# SCENARIO 12
################################

## Upload the parameter Excel files.
phase1 <- read_xlsx("Scenario 12.xlsx", sheet = "Phase1")
phase2 <- read_xlsx("Scenario 12.xlsx", sheet = "Phase2")
phase3 <- read_xlsx("Scenario 12.xlsx", sheet = "Phase3")
phase4 <- read_xlsx("Scenario 12.xlsx", sheet = "Phase4")
phase5 <- read_xlsx("Scenario 12.xlsx", sheet = "Phase5")

# Generate the Data

# Specify the simulation settings
for_sim <- c(list(.x = 1:n_sim, .f=each_sim), settings)
# Setting seed once at the beginning of the simulation
set.seed(1234)
# Simulate the data
all_sims <- do.call(purrr::map, args = for_sim)

### Save the file as an RDS file

saveRDS(all_sims, file = "scenario12.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("phase1","phase2","phase3","phase4","phase5",
                 "for_sim", "all_sims")

rm(list = data_delete)
