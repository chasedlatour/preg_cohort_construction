# Calculate expected values among the untreated

## Expected values for the severity distribution among the untreated

# Severity distribution in overall cohort
p_sev_dist = c(1/3, 1/3, 1/3)

# Probability of treatment dependent upon severity value
p_trt_sev = c(0.4, 0.50, 0.6)

# Untreated severity distribution over the whole cohort
p_sev_whole_cohort <- p_sev_dist * p_trt_sev

# Severity distribution among the untreated

p_sev_untreated <- p_sev_whole_cohort / sum(p_sev_whole_cohort)

# Expected PNC distribution among those who are untreated
## Note: Timing of PNC initiation is unrelated to severity

p_indx_pnc <- c(0.23, 0.33, 0.44)

# Probability of low severity indices
p_low_indx <- p_sev_untreated[1]*p_indx_pnc
p_med_indx <- p_sev_untreated[2]*p_indx_pnc
p_high_indx <- p_sev_untreated[3]*p_indx_pnc

# Confirm that sum to 1: sum(p_low_indx, p_med_indx, p_high_indx) -- Good
param_file <- "Parameters_Abortion0_8_Preeclampsia0_8.xlsx"
potential_preg_untrt <- read_xlsx(param_file, sheet = "potential_preg_untrt")

