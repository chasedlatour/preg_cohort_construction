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
  seed = 13049857
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.



# Create the data generation parameters
treatment_effects <- expand.grid(
  rr_abortion = c(0.8), #c(0.5, 0.8, 1, 1.25, 2),
  rr_preec = c(0.8) #c(0.5, 0.8, 1)
) 
treatment_effects$param_file = paste0("Parameters_Abortion",
                                      gsub("\\.", "", as.character(treatment_effects$rr_abortion)), 
                                      "_Preeclampsia", 
                                      gsub("\\.", "", as.character(treatment_effects$rr_preec)), 
                                      ".xlsx")
treatment_effects$n_sim = 1
treatment_effects$n = 20000


# Create a dataset with the missing data parameters
missing_params <- data.frame(
  marginal_p_miss_severity = 1.13 * c(0.0125, 0.025, 0.0375, 0.05, 0.1, 0.15),
  beta12 = 0.7,
  marginal_p_miss_miscarriage = 5 * c(0.0375, 0.025, 0.0125, 0.15, 0.1, 0.05),
  gamma1 = -0.2
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
            gamma1 = gamma1
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
            gamma1 = gamma1
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
            gamma1 = gamma1
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
