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
  packages = c("tibble", "tidyverse", "readxl", "survival", "rlang", "data.table",
               "cmprsk"),
  seed = 13049857 #,
  # format = "qs", # Optionally set the default storage format. qs is fast.
  # controller = crew::crew_controller_local(workers = 2, seconds_idle = 3),
  # memory = "transient",
  # storage = "worker",
  # retrieval = "worker",
  # garbage_collection = TRUE
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()



# Create the data generation parameters
treatment_effects <- expand.grid(
  rr_abortion = c(0.5, 1, 2),
  rr_preec = c(0.5, 1) 
) 
treatment_effects$param_file = paste0("Parameters_Abortion",
                                      gsub("\\.", "", as.character(treatment_effects$rr_abortion)), 
                                      "_Preeclampsia", 
                                      gsub("\\.", "", as.character(treatment_effects$rr_preec)), 
                                      "_EMM.xlsx")
treatment_effects$n_sim = 1
treatment_effects$n = 10000


# Create a dataset with the missing data parameters
missing_params <- data.frame(
  # marginal_p_miss_severity = 1.13 * c(0, 0.025, 0.05,
  #                                     0, 0.100, 0.20),
  marginal_p_miss_severity = c(0, 0.025, 0.05,
                               0, 0.100, 0.20) * 0.3,
  beta12 = 2.55, 
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
        # Prep the data for analysis - put in data.table structure for speed.
        targets::tar_target(
          data_prep,
          as.data.table(cohort_data) %>% 
            prep_data_for_analysis(pnc_wk = pnc_wk)
        ),
        # Derive the severity distribution among all pregnancies
        targets::tar_target(
          severity_dist_total,
          sev_dist(data_prep)
        ),
        # Derive the severity distribution among all pregnancies with observed outcomes
        targets::tar_target(
          severity_dist_outcomes,
          sev_dist(data_prep[obs_outcome_mar_mnar == 1])
        ),
        # Derive the severity distribution among all pregnancies with observed deliveries
        targets::tar_target(
          severity_dist_deliveries,
          sev_dist(data_prep[obs_delivery_mar_mnar == 1])
        ),
        # Get the risks using the potential outcomes
        targets::tar_target(
          potential_outcomes, 
          potential_risks(data_prep)
        ),
        # Get the risks among the observed deliveries
        targets::tar_target(
          observed_deliveries,
          calculate_risks(data_prep[obs_delivery_mar_mnar == 1], severity_dist_deliveries) %>% 
            as_tibble()
        ),
        # Get the risks among the pregnancies with observed outcomes
        targets::tar_target(
          observed_outcomes,
           calculate_risks(data_prep[obs_outcome_mar_mnar == 1], severity_dist_outcomes) %>% 
            as_tibble()
        ),
        # Get the risks among high-severity pregnancies using Aalen-Johansen estimator
        targets::tar_target(
          high_sev_aj,
          calculate_aj_risks_sev(data_prep[severity == 2], severity_dist_total)
        ),
        # Get the risks among moderate severity pregnancies using Aalen-Johansen estimator
        targets::tar_target(
          middle_sev_aj,
          calculate_aj_risks_sev(data_prep[severity == 1], severity_dist_total)
        ),
        # Get the risks among low severity pregnancies using Aalen-Johansen estimator
        targets::tar_target(
          low_sev_aj,
          calculate_aj_risks_sev(data_prep[severity == 0], severity_dist_total)
        ),
        # Get the risks among all pregnancies using the AJ estimator - standardized by severity distribution
        targets::tar_target(
          tte_mar_mnar,
          calculate_aj_risks(high_sev_aj, middle_sev_aj, low_sev_aj) %>%
            as_tibble()
        ),
        # Get the risks for the sensitivity analyses
        targets::tar_target(
          sensitivity_analyses,
          sensitivity_analysis_risks(data_prep, severity_dist_total)
        ),
        # Get the tibble for all outcomes - Assume that all missing outcomes have an outcome
        targets::tar_target(
          sens_anal_mar_mnar_all_outc,
          sensitivity_analyses %>% 
            select(trt, all_outc) %>% 
            pivot_wider(names_from = trt,
                        values_from = all_outc,
                        names_prefix = "risk") %>% 
            mutate(
              rd = risk1 - risk0,
              rr = risk1 / risk0
            )
        ),
        # Get the tibble for all trt outcomes - Assume that treated missing have an outcome
        targets::tar_target(
          sens_anal_mar_mnar_trt_outc,
          sensitivity_analyses %>% 
            select(trt, trt_outc) %>% 
            pivot_wider(names_from = trt,
                        values_from = trt_outc,
                        names_prefix = "risk") %>% 
            mutate(
              rd = risk1 - risk0,
              rr = risk1 / risk0
            )
        ),
        # Get the tibble for all untrt outcomes - Assume that untreated missing have an outcome
        targets::tar_target(
          sens_anal_mar_mnar_untrt_outc,
          sensitivity_analyses %>% 
            select(trt, untrt_outc) %>% 
            pivot_wider(names_from = trt,
                        values_from = untrt_outc,
                        names_prefix = "risk") %>% 
            mutate(
              rd = risk1 - risk0,
              rr = risk1 / risk0
            )
        ),
        # Get the tibble for all no outcomes - Assume that all missing outcomes do not have an outcome
        targets::tar_target(
          sens_anal_mar_mnar_no_outc,
          sensitivity_analyses %>% 
            select(trt, no_outc) %>% 
            pivot_wider(names_from = trt,
                        values_from = no_outc,
                        names_prefix = "risk") %>% 
            mutate(
              rd = risk1 - risk0,
              rr = risk1 / risk0
            )
        ),
        # Combine the analyzed data
        targets::tar_target(
          analyzed_data,
          tibble(
            sev_dist_all_pregnancies = list(severity_dist_total),
            sev_dist_observed_outcomes = list(severity_dist_outcomes),
            sev_dist_observed_deliveries = list(severity_dist_deliveries),
            potential_risks = list(potential_outcomes),
            observed_deliveries = list(observed_deliveries),
            observed_outcomes = list(observed_outcomes),
            tte = list(tte_mar_mnar),
            sens_anal_all_outc = list(sens_anal_mar_mnar_all_outc),
            sens_anal_trt_outc = list(sens_anal_mar_mnar_trt_outc),
            sens_anal_untrt_out = list(sens_anal_mar_mnar_untrt_outc),
            sens_anal_no_outc = list(sens_anal_mar_mnar_no_outc)
          ) %>% 
            mutate(
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

