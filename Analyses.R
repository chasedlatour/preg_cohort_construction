##############################################
# Program: Analyses.R
# Programmer: Chase Latour
# Date Last Modified:
#
# Purpose: Run the analyses for the different
# scenarios.
##############################################

## Load in the Libraries that need
library(tidyverse)
library(survival)
library(Epi)
library(kableExtra)

## Load in the functions that need.
source('analysis functions.R')





#############################################
## SCENARIO 1
#############################################

## Load in the data
data <- readRDS('scenario1.rds')
# For testing: data <- readRDS('test.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario1.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario1.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario1.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





#############################################
## SCENARIO 2
#############################################

## Load in the data
data <- readRDS('scenario2.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario2.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario2.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario2.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)






#############################################
## SCENARIO 3
#############################################

## Load in the data
data <- readRDS('scenario3.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario3.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario3.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario3.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





#############################################
## SCENARIO 4
#############################################

## Load in the data
data <- readRDS('scenario4.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario4.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario4.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario4.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)






#############################################
## SCENARIO 5
#############################################

## Load in the data
data <- readRDS('scenario5.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario5.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario5.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario5.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)







#############################################
## SCENARIO 6
#############################################

## Load in the data
data <- readRDS('scenario6.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario6.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario6.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario6.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)




#############################################
## SCENARIO 7
#############################################

## Load in the data
data <- readRDS('scenario7.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario7.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario7.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario7.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)






#############################################
## SCENARIO 8
#############################################

## Load in the data
data <- readRDS('scenario8.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario8.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario8.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario8.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





#############################################
## SCENARIO 9
#############################################

## Load in the data
data <- readRDS('scenario9.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario9.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario9.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario9.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





#############################################
## SCENARIO 10
#############################################

## Load in the data
data <- readRDS('scenario10.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario10.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario10.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario10.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





#############################################
## SCENARIO 11
#############################################

## Load in the data
data <- readRDS('scenario11.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario11.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario11.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario11.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)






#############################################
## SCENARIO 12
#############################################

## Load in the data
data <- readRDS('scenario12.rds')

# Clean the data for analyses
set.seed(1234)
trial <- trial_cohort(data)

## Describe the cleaned data
trial_descriptive <- describe_trial(trial)
saveRDS(trial_descriptive, "describe trial_scenario12.rds")

## Describe the outcomes in the cleaned data
trial_outcomes <- describe_outcomes(trial)
saveRDS(trial_outcomes, "describe outcomes_scenario12.rds")

## Analyze the cleaned data
analyses <- clean_analyze(trial)
saveRDS(analyses, "analyses_scenario12.rds")

## Remove datasets from the environment to ensure no issues

data_delete <- c("data", "trial", "trial_descriptive", 
                 "trial_outcomes","analyses")

rm(list = data_delete)





