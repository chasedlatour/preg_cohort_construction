# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.


## For Installing packages on N2 server:
# install.packages("igraph", repos = c("https://igraph.r-universe.dev", "https://cloud.r-project.org"))
# install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz", repos=NULL, type="source")
# install.packages("https://cran.r-project.org/src/contrib/etm_1.1.1.tar.gz", repos = NULL, type="source")


# Set target options:
tar_option_set(
  # Packages that your targets need for their tasks.
  packages = c("tibble", "tidyverse", "readxl", "survival", "rlang"),
  seed = 13049857,
  # format = "qs", # Optionally set the default storage format. qs is fast.
  controller = crew::crew_controller_local(workers = 2, seconds_idle = 3),
  # memory = "transient",
  storage = "worker",
  retrieval = "worker",
  garbage_collection = TRUE
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.



# Create the data generation parameters
treatment_effects <- expand.grid(
  rr_abortion = c(0.7, 1, 1.5),# c(0.8), #
  rr_preec = c(0.7, 1) #c(0.8) #
) 
treatment_effects$param_file = paste0("Parameters_Abortion",
                                      gsub("\\.", "", as.character(treatment_effects$rr_abortion)), 
                                      "_Preeclampsia", 
                                      gsub("\\.", "", as.character(treatment_effects$rr_preec)), 
                                      ".xlsx")
treatment_effects$n_sim = 1
treatment_effects$n = 10000


# Create a dataset with the missing data parameters
missing_params <- data.frame(
  marginal_p_miss_severity = 1.13 * c(0, 0.025, 0.05,
                                      0, 0.100, 0.20),
  beta12 = 0.7,
  marginal_p_miss_miscarriage = c(4.9 * c(0.05, 0.025, 0,
                                          0.20, 0.100, 0)),
  # Cut it down from 5 so not too large
  gamma1 = -0.2,
  pnc_wk = rep(7)
)




# Replace the target list below with your own:
list(
  
  # Generate the data and create the cohorts
  mapped1 <- tarchetypes::tar_map(
    unlist = FALSE,
    values = treatment_effects,
    list(
      # Generate datasets
      targets::tar_target(
        generated_data,
        generate_dgm(
          n_sim = n_sim,
          n = n,
          param_file = param_file,
          rr_abortion = rr_abortion,
          rr_preec = rr_preec
        )
      ),
      # Create the cohorts
      mapped2 <- tarchetypes::tar_map(
        unlist = FALSE,
        values = missing_params,
        targets::tar_target(
          cohort_data,
          generate_cohort(
            generated_data,
            marginal_p_miss_severity = marginal_p_miss_severity,
            beta12 = beta12,
            marginal_p_miss_miscarriage = marginal_p_miss_miscarriage,
            gamma1 = gamma1,
            pnc_wk = pnc_wk
          )
          ),
        # Describe the data
        targets::tar_target(
          described_data,
          describe_cohort(
            cohort_data,
            rr_abortion = rr_abortion,
            rr_preec = rr_preec,
            marginal_p_miss_severity = marginal_p_miss_severity,
            beta12 = beta12,
            marginal_p_miss_miscarriage = marginal_p_miss_miscarriage,
            gamma1 = gamma1,
            pnc_wk = pnc_wk
          )
        ),
        # Analyze the data
        targets::tar_target(
          analyzed_data,
          run_analysis(
            cohort_data,
            rr_abortion = rr_abortion,
            rr_preec = rr_preec,
            marginal_p_miss_severity = marginal_p_miss_severity,
            beta12 = beta12,
            marginal_p_miss_miscarriage = marginal_p_miss_miscarriage,
            gamma1 = gamma1,
            pnc_wk = pnc_wk
          )
        )
      ),
      # Combine all of the descriptive statistics together
      tarchetypes::tar_combine(
        scenario_descriptives,
        mapped2[["described_data"]],
        command = dplyr::bind_rows(!!!.x)
      ),
      # Combine all of the analyses together
      tarchetypes::tar_combine(
        scenario_analyses,
        mapped2[["analyzed_data"]],
        command = dplyr::bind_rows(!!!.x)
      )
    )
  ),
  
  # Combine all of the descriptive results into one dataset
  tar_combine(
    all_descriptives,
    mapped1[["scenario_descriptives"]],
    command = dplyr::bind_rows(!!!.x)
  ),
  
  # Combine all of the analysis results into one dataset 
  tar_combine(
    all_analyses,
    mapped1[["scenario_analyses"]],
    command = dplyr::bind_rows(!!!.x)
  ),
  
  # Output a RMD file with all of the output.
  tar_render(report, "Describe Cohorts.Rmd")
  
  
)
