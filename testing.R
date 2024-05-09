
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
all_outcomes <- generate(n_sim = 1, 
                         n = 1200, 
                         p_sev_beta = c(-3, 0.1, 0.2),
                         p_trt_sev = c(0.35, 0.50, 0.65),
                         p_indx_pnc = c(0.25, 0.375, 0.375)
                         )

p_miss_outcome <- c(0.2, -0.2)

# Identify key variables and then output a dataset with 3 columns - one for each trial.
all_outcomes2 <- create_cohort(all_outcomes, p_miss_outcome)



#####################################
### Look at dataset summaries


# Week 4
wk4 <- all_outcomes2[[2]][[1]]
# Week 7
wk7 <- all_outcomes2[[3]][[1]]
# Week 16
wk16 <- all_outcomes2[[4]][[1]]

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

# Look at the distribution of treatment
table(wk4$trt)
table(wk4$trt, wk4$severity)


