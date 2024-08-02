###############################################
# Program: describe cohorts functions.R
# Programmer: Chase Latour
#
# Purpose: The purpose of this program is to 
# store functions that will be use to describe
# the individuals included in the analyses.
###############################################


describe_cohort <- function(dataset){
  
  hold <- dataset %>%
    summarize(
      
      # Calculate the number per category
      n_low = sum(severity == 0),
      n_low_perc = round(100*n_low/n(),2),
      n_med = sum(severity == 1),
      n_med_perc = round(100*n_med/n(), 2),
      n_high = sum(severity == 2),
      n_high_perc = round(100*n_high/n(), 2),
      
      # Distribution by first pnc encounter
      n_pnc_4 = sum(pnc_wk == 4),
      pnc4_perc = round(100*n_pnc_4/n(),2),
      n_pnc_7 = sum(pnc_wk == 7),
      pnc7_perc = round(100*n_pnc_7/n(),2),
      n_pnc_16 = sum(pnc_wk == 16),
      pnc16_perc = round(100*n_pnc_16/n(),2),
      
      # Preterm birth
      n_preterm_all = sum(pregout_t_pre_miss >= 18 & pregout_t_pre_miss < 35),
      preterm_p = round(100*n_preterm_all / n(), 2),
      
      # Distribution of outcomes
      n_miscarriage_all = sum(pregout_pre_miss == "fetaldeath" & pregout_t_pre_miss < 18),
      miscarriageall_p = round(100*n_miscarriage_all / n(),2),
      n_stillbirth_all = sum(pregout_pre_miss == "fetaldeath" & pregout_t_pre_miss >= 18),
      stillbirthall_p = round(100*n_stillbirth_all/n(), 2),
      n_livebirth_all = sum(pregout_pre_miss == "livebirth"),
      livebirthall_p = round(100*n_livebirth_all/n()),
      
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
      n_preec_l32_p = round(100*n_preec_l32/n(), 2),
      n_preec = sum(preeclampsia_pre_miss == 1),
      n_preec_perc = round(100*n_preec/n(),2),
      preec_25 = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.25)+2,
      preec_75 = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.75)+2,
      preec_med = quantile(pregout_t_pre_miss[preeclampsia_pre_miss == 1], probs=0.50)+2,
      
      # Distribution of uncensored outcomes
      missing_mar = sum(ltfu_mar != "not"),
      missing_mar_perc = round(100*missing_mar/n(), 2),
      missing_mar_p25 = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.25)+2,
      missing_mar_p75 = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.75)+2,
      missing_mar_med = quantile(pregout_t_mar[ltfu_mar != "not"], probs = 0.5)+2,
      
      missing_mnar = sum(ltfu_mar_mnar != "not"),
      missing_mnar_perc = round(100*missing_mnar/n(), 2),
      missing_mnar_p25 = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not"], probs=0.25)+2,
      missing_mnar_p75 = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not"], probs = 0.75)+2,
      missing_mnar_med = quantile(pregout_t_mar_mnar[ltfu_mar_mnar != "not"], probs = 0.5)+2,
      
      .groups = 'keep'
    ) %>% 
    mutate(
      # Number of pregnancies per severity
      n_low = paste0(n_low," (", n_low_perc, "%)"),
      n_med = paste0(n_med," (", n_med_perc, "%)"),
      n_high = paste0(n_high," (", n_high_perc, "%)"),
      # Number of pregnancies per pnc
      n_pnc_4 = paste0(n_pnc_4, " (", pnc4_perc, "%)"),
      n_pnc_7 = paste0(n_pnc_7, " (", pnc7_perc, "%)"),
      n_pnc_16 = paste0(n_pnc_16, " (", pnc4_perc, "%)"),
      # Preterm birth
      n_preterm_all = paste0(n_preterm_all, " (", preterm_p, ")%"),
      # Pregnancy outcomes
      n_miscarriage_all = paste0(n_miscarriage_all, " (", miscarriageall_p,")%"),
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
    select(n_low, n_med, n_high, n_pnc_4, n_pnc_7, n_pnc_16,
           n_preterm_all, n_miscarriage_all,
           n_stillbirth_all, n_livebirth_all, miscarriage_t, stillbirth_t, livebirth_t,
           n_preeclampsia, preeclampsia_t, n_preeclampsia_l32,
           missing_mar, missing_mar_t, missing_mnar, missing_mnar_t
    ) %>% 
    pivot_longer(
      cols = everything(),
      names_to = "SummaryStatistic",
      values_to = "value"
    ) %>% 
    pivot_wider(
      names_from = SummaryStatistic,
      values_from = value
    )
  
  return(hold)
  
}