# Test covariates new analysis: 4/8/2025

# Set up data: 
source(fs::path(here::here("!libraries.R")))
source(fs::path(here::here("!directories.R")))
source(fs::path(here::here("!functions.R")))
source(fs::path(here::here("!load_data.R")))


# Test covariates
covariates <- c(NULL, "q1_dp_amds_totscore", "alcohol_intake", "smokestatus_imputed",
                "q1_fdgp7", #"met_lg2_Creatinine", 
                "alcohol_intake_cat", "q1_dgm165", 
                "diabetes", "blood_age", #"year_between", 
                "q1_dgm156", "q1_fdgp6", "q1_fdgp6", "q1_d_fdgp6", "q1_fdgp7", 
                "q1_d_fdgp7", "q1_dgm167", "q1_dgm168", "q1_dp_ahei2010_totscore",
                "q1_dp_amds_totscore", "q1_dp_amds_e_totscore", 
                "q1_dp_dash_totscore", "q1_dii_density", 
                "q1_dp_hei2015_totscore", "q1_edip", "q1_edih", "q1_edir", 
                "q1_elih_imputed", "q1_elir_imputed", 
                colnames(data_hcc |> tidylog::select(prs_uw:prs_uw_tertile)))



covariates <- unique(covariates)


# Function to perform the analysis with one covariate at a time
run_single_covariate_analysis <- function(covariate) {
  res <- data_hcc |> 
    # mutate(across(all_of(pfas_name), ~exp(.x) |> scale() |> as.numeric())) |>
    epiomics::owas_clogit(cc_set = "setnum",
                          cc_status = "status",
                          covars = covariate,
                          omics = pfas_name,
                          conf_int = TRUE)
  res$covariate <- covariate  # add covariate name to result for clarity
  return(res)
}

# Apply the function across all covariates
results_list <- map(covariates, run_single_covariate_analysis)

# Combine results into a single dataframe
covar_mod_results <- bind_rows(results_list)

# Calculate percent change from base model ----------
base_res <- data_hcc |> 
  # mutate(across(all_of(pfas_name), ~exp(.x) |> scale() |> as.numeric())) |>
  epiomics::owas_clogit(cc_set = "setnum",
                        cc_status = "status",
                        covars = NULL,
                        omics = pfas_name,
                        conf_int = TRUE)
# res$covariate <- covariate  # add covariate name to result for clarity

base_est <- base_res |> 
  tidylog::select(feature_name, estimate) |>
  rename(estimate_basemod = estimate)

# Join with covariate models
final_results <- tidylog::full_join(base_est, covar_mod_results)

# calculate the percent change in the estimate from baseline model
final_results <- final_results |> 
  mutate(abs_percent_change = abs(100*(estimate-estimate_basemod)/estimate_basemod))

hist(final_results$abs_percent_change)
table(final_results$abs_percent_change<10)
# covariates that change the estimates greater than 10%
potential_covars <- final_results |> tidylog::filter(abs_percent_change>10)

# There are 12 PFAS, so lets select covariates that change the estimate by >10% in at least half
covar_summary <- potential_covars |> 
  group_by(covariate) |> 
  summarise(n = n(), 
            pct_pfas = n/12,
            median_pct_chg = median(abs_percent_change),
            pfas = str_c(feature_name, collapse = "; "))|> 
  arrange(desc(n))

# Categories of covariates 
covar_summary <- covar_summary |> 
  mutate(category = case_when(
    str_detect(covariate, "prs") ~ "genetic", 
    str_detect(covariate, "q1_") ~ "diet", 
    TRUE ~ "Other")) |> 
  arrange(category, -median_pct_chg) 

# Filter to covariates that change the estimate by >10% in at least half
View(covar_summary |> tidylog::filter(pct_pfas>0.5))

# Plot results ----------
final_results |> 
  mutate(feature_name = fct_reorder(feature_name, estimate)) |>
  ggplot(aes(y = feature_name, 
             x = estimate, 
             color = p_value < 0.05)) +
  geom_point() +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2) +
  facet_wrap(~covariate) + 
  cowplot::theme_cowplot() +
  geom_vline(xintercept = 0) + 
  labs(title = "Covariate Analysis Results",
       x = "Covariate",
       y = "Estimate") +
  scale_color_manual(values = c( "black", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



# Plot results differently
ggplot(final_results, aes(x = covariate, y = estimate, color = p_value < 0.05)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2) +
  facet_wrap(~feature_name) + 
  cowplot::theme_cowplot() +
  geom_hline(yintercept = 0) + 
  labs(title = "Covariate Analysis Results",
       x = "Covariate",
       y = "Estimate") +
  scale_color_manual(values = c( "black", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 




# Test individual analyses:
# PFHpA quartiles
resout <- clogit(status ~ as.numeric(pfas_t_PFHPA_quartile) +
                   diabetes + 
                   smokestatus_imputed + 
                   q1_edih +
                   prs_wt + v1 + v2 + v3 + v4 + v5 + v6 +
                   strata(setnum), 
                 data = data_hcc)

broom::tidy(resout)
