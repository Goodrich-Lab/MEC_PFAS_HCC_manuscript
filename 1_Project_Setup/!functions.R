# Rename PFAS -----------------------
rename_pfas <- function(pfas_names){
  x <- tibble(pfas = pfas_names)
  pfas2 <-  x %>%
    mutate(pfas = str_remove_all(pfas, "pfas_t_")) %>%
    mutate(pfas = case_when(pfas == "PFHPA" ~ "PFHpA",
                             pfas == "PFHXA" ~ "PFHxA",
                             pfas == "PFHXS" ~ "PFHxS",
                             pfas == "PFUNDA" ~ "PFUnDA",
                             pfas == "PFHPS" ~ "PFHpS",
                             pfas == "NMFOSAA" ~ "NMeFOSAA",
                             pfas == "NETFOSAA" ~"NEtFOSAA",
                             TRUE ~ pfas))
  return(pfas2$pfas)
}

# Model with interaction of prs_wt_std and pfas
model_interaction <- function(pfas){
  # Formula
  formula_with_int <- as.formula(
    paste("status ~", paste0(pfas, "*prs_wt_std"),
          "+ v1 + v2 + v3 + v4 + v5 + v6 +", 
          paste(covars[1:3], collapse = " + "),
          "+ strata(setnum)"))
  # Model
  model_with_int <- clogit(formula_with_int, 
                           data = data_hcc, 
                           method = "efron", 
                           robust = TRUE)
  # Interaction effect at high/low PRS
  res <- emtrends(model_with_int, ~prs_wt_std, var = pfas, 
                  at = list(prs_wt_std = c(cuts[1], cuts[9]), 
                            smokestatus_imputed = "Never",
                            diabetes = "No",
                            q1_edih = 0,
                            v1 = 0,
                            v2 = 0,
                            v3 = 0,
                            v4 = 0,
                            v5 = 0,
                            v6 = 0),
                  transform = NULL) %>%
    data.frame() %>%
    select(prs_wt_std, ends_with(".trend"), asymp.LCL, asymp.UCL) %>%
    rename(estimate = ends_with(".trend"),
           conf_low = asymp.LCL,
           conf_high = asymp.UCL) %>%
    mutate(type = ifelse(prs_wt_std == cuts[1], "Low Genetic Risk", "High Genetic Risk"))
  
  # Main effect without interaction
  res_main <- epiomics::owas_clogit(data_hcc , 
                                    cc_set = "setnum", 
                                    cc_status = "status", 
                                    covars = covars,
                                    omics = pfas,
                                    conf_int = TRUE) %>%
    mutate(type = "Overall (no interaction with PRS)")
  # Combine results
  res_combined <- bind_rows(res_main, res) %>%
    mutate(odds_ratio = exp(estimate),
           exp_ci_low = exp(conf_low),
           exp_ci_high = exp(conf_high),
           group = if_else(str_detect(type, "Overall"), "overall", "with genetic"),
           pfas = pfas)
  
  res_combined
}
