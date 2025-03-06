###############################################
# Program: describe cohorts functions.R
# Programmer: Chase Latour
#
# Purpose: The purpose of this program is to 
# store functions that will be use to describe
# the individuals included in the analyses.
###############################################


describe_cohort <- function(dataset, rr_abortion, rr_preec,
                            marginal_p_miss_severity, beta12,
                            marginal_p_miss_miscarriage, gamma1,
                            pnc_wk){
  
  hold <- dataset %>%
    group_by(trt) %>% 
    summarize(
      
      # Calculate the number per category
      n_low = sum(severity == 0),
      n_low_perc = round(100*n_low/n(), 0),
      n_med = sum(severity == 1),
      n_med_perc = round(100*n_med/n(), 0),
      n_high = sum(severity == 2),
      n_high_perc = round(100*n_high/n(), 0),
      
      # Calcualte the number of rural per category
      n_rural = sum(rural == 1),
      n_rural_perc = round(100*n_rural/n(), 0),
      
      # Preterm birth
      n_preterm_all = sum(pregout_t_pre_miss >= 18 & pregout_t_pre_miss < 35),
      preterm_p = round(100*n_preterm_all / n(), 0),
      
      # Distribution of outcomes
      n_miscarriage_all = sum(pregout_pre_miss == "fetaldeath" & pregout_t_pre_miss < 18),
      miscarriageall_p = round(100*n_miscarriage_all / n(), 0),
      n_stillbirth_all = sum(pregout_pre_miss == "fetaldeath" & pregout_t_pre_miss >= 18),
      stillbirthall_p = round(100*n_stillbirth_all/n(), 0),
      n_livebirth_all = sum(pregout_pre_miss == "livebirth"),
      livebirthall_p = round(100*n_livebirth_all/n(), 0),
      
      # Timing of those outcomes
      miscarriageall_p25 = quantile(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                         pregout_t_pre_miss < 18], 
                                    probs = 0.25)+2,
      miscarriageall_p75 = quantile(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                         pregout_t_pre_miss < 18],
                                    probs = 0.75)+2,
      miscarriageall_med = median(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                       pregout_t_pre_miss < 18])+2,
      
      stillbirthall_p25 = quantile(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                        pregout_t_pre_miss >= 18],
                                   probs = 0.25)+2,
      stillbirthall_p75 = quantile(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                        pregout_t_pre_miss >= 18],
                                   0.75)+2,
      stillbirthall_med = median(pregout_t_pre_miss[pregout_pre_miss == "fetaldeath" &
                                                      pregout_t_pre_miss >= 18])+2,
      
      livebirth_p25 = quantile(pregout_t_pre_miss[pregout_pre_miss == "livebirth"],
                               probs = 0.25)+2,
      livebirth_p75 = quantile(pregout_t_pre_miss[pregout_pre_miss == "livebirth"],
                               probs = 0.75)+2,
      livebirth_med = median(pregout_t_pre_miss[pregout_pre_miss == "livebirth"])+2,
      
      # Distribution of Preeclampsia
      n_preec_l32 = sum(preeclampsia_pre_miss == 1 & pregout_t_pre_miss < 30),
      n_preec_l32_p = round(100*n_preec_l32/n(), 0),
      n_preec = sum(preeclampsia_pre_miss == 1),
      n_preec_perc = round(100*n_preec/n(),0),
      preec_25 = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.25)+2,
      preec_75 = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.75)+2,
      preec_med = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.50)+2,
      
      # Distribution of uncensored outcomes
      missing_mar = sum(ltfu_mar != "not"),
      missing_mar_perc = round(100*missing_mar/n(), 0),
      missing_mar_p25 = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.25)+2,
      missing_mar_p75 = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.75)+2,
      missing_mar_med = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.5)+2,
      
      missing_mnar = sum(ltfu_mar_mnar != "not") - sum(ltfu_mar != "not"),
      missing_mnar_perc = round(100*missing_mnar/n(), 0),
      missing_mnar_p25 = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not" & ltfu_mar == "not"], probs=0.25)+2,
      missing_mnar_p75 = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not" & ltfu_mar == "not"], probs = 0.75)+2,
      missing_mnar_med = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not" & ltfu_mar == "not"], probs = 0.5)+2,
      
      .groups = 'keep'
    ) %>% 
    mutate(
      # Number of pregnancies per severity
      n_low = paste0(n_low," (", n_low_perc, "%)"),
      n_med = paste0(n_med," (", n_med_perc, "%)"),
      n_high = paste0(n_high," (", n_high_perc, "%)"),
      # Rurality
      n_rural = paste0(n_rural, " (", n_rural_perc, "%)"),
      # Preterm birth
      n_preterm_all = paste0(n_preterm_all, " (", preterm_p, "%)"),
      # Pregnancy outcomes
      n_miscarriage_all = paste0(n_miscarriage_all, " (", miscarriageall_p,"%)"),
      n_stillbirth_all = paste0(n_stillbirth_all, " (", stillbirthall_p, "%)"),
      n_livebirth_all = paste0(n_livebirth_all, " (", livebirthall_p, "%)"),
      # Timing of Preg Outcomes
      miscarriage_t = paste0(miscarriageall_med, " (", miscarriageall_p25, ", ", miscarriageall_p75, ")"),
      stillbirth_t = paste0(stillbirthall_med, " (", stillbirthall_p25, ", ", stillbirthall_p75, ")"),
      livebirth_t = paste0(livebirth_med, " (", livebirth_p25, ", ", livebirth_p75, ")"),
      # Preeclampsia
      n_preeclampsia = paste0(n_preec, " (", n_preec_perc, "%)"),
      preeclampsia_t = paste0(preec_med, " (", preec_25, ", ", preec_75, ")"),
      # Preeclampsia less than 32 weeks
      n_preeclampsia_l32 = paste0(n_preec_l32, " (", n_preec_l32_p, "%)"),
      # Missing by MAR
      missing_mar = paste0(missing_mar, " (", missing_mar_perc, "%)"),
      missing_mar_t = paste0(missing_mar_med, " (", missing_mar_p25, ", ", missing_mar_p75, ")"),
      # Missing by MNAR
      missing_mnar = paste0(missing_mnar, " (", missing_mnar_perc, "%)"),
      missing_mnar_t = paste0(missing_mnar_med, " (", missing_mnar_p25, ", ", missing_mnar_p75, ")") 
    ) %>% 
    select(trt, n_low, n_med, n_high, n_rural,
           n_preterm_all, n_miscarriage_all,
           n_stillbirth_all, n_livebirth_all, miscarriage_t, stillbirth_t, livebirth_t,
           n_preeclampsia, preeclampsia_t, n_preeclampsia_l32,
           missing_mar, missing_mar_t, missing_mnar, missing_mnar_t
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
  
  return(hold)
  
}


# cohort_data <- tar_read(cohort_data_0.1_1.3_0.49_.0.2_7_0.5_1_Parameters_Abortion05_Preeclampsia1_EMM.xlsx_1_10000)

describe_analytic <- function(cohort_data, rr_abortion, rr_preec,
                              marginal_p_miss_severity, beta12,
                              marginal_p_miss_miscarriage, gamma1,
                              pnc_wk){
  
  # Number of pregnancies by treatment arm among observed deliveries
  observed_deliveries <- cohort_data %>% 
    filter(ltfu_mar_mnar == 'not' & pregout_t_mar_mnar >= 18) %>% 
    group_by(trt) %>% 
    summarize(
      n_pregnancies_deliveries = n(),
      n_preeclampsia_deliveries = sum(preeclampsia_mar_mnar)
    )
  
  # Number of pregnancies by treatment arm among observed outcomes
  observed_outcomes <- cohort_data %>% 
    filter(ltfu_mar_mnar == 'not') %>% 
    group_by(trt) %>% 
    summarize(
      n_pregnancies_outcomes = n(),
      n_preeclampsia_outcomes = sum(preeclampsia_pre_miss)
    )
  
  # Number of pregnancies by treatment arm among all pregnancies
  observed_pregnancies <- cohort_data %>% 
    group_by(trt) %>% 
    summarize(
      n_pregnancies_all = n(),
      n_preeclampsia_all = sum(preeclampsia_pre_miss)
    )
  
  observed_n <- left_join(observed_deliveries, observed_outcomes, by = "trt") %>% 
    left_join(observed_pregnancies, by = "trt") %>% 
    # Add in identifying variables
    mutate(
      rr_abortion = rr_abortion, 
      rr_preec = rr_preec,
      marginal_p_miss_severity = marginal_p_miss_severity, 
      beta12 = beta12,
      marginal_p_miss_miscarriage = marginal_p_miss_miscarriage, 
      gamma1 = gamma1,
      pnc_wk = pnc_wk
    )
  
  return(observed_n)
  
}


# cohort_data %>% 
#   group_by(ltfu_mar_mnar, pregout_mar_mnar) %>% 
#   summarize(
#     n = n(),
#     min = min(pregout_t_mar_mnar),
#     max = max(pregout_t_mar_mnar),
#     preeclampsia = sum(preeclampsia_pre_miss),
#     preeclampsia_miss = sum(preeclampsia_mar_mnar)
#   )
