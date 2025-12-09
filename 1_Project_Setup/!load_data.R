# Set up data for analysis
# library(tidylog)


# Read Data ---------------------------------------------------------------
data <- readRDS(fs::path(
  dir_cleaned_data,
  "hhear_complete_data_v5.RDS")) %>%
  mutate(pfas_t_PFBS_detected = ifelse(pfas_t_PFBS_detected == "Detected", 1, 0))


snp_info <- readxl::read_xlsx(fs::path(dir_data %>% dirname() %>% dirname() %>% dirname(),
                                       "0_data_processing",
                                       "0_raw_data",
                                       "diet_PRS_data",
                                       "snp_info.xlsx")) %>%
  janitor::clean_names() %>%
  mutate(snp_name = paste0("x", chr, "_", pos38))


# untargeted_met_data <- readRDS(fs::path(dir_cleaned_data, "Untargeted_feature_metadata.RDS"))

# variables--------------
untargeted_pfas_name <- colnames(data %>% 
                        dplyr::select(contains("_ut_"), -contains("_detected")))

pfas_name <- colnames(data %>% 
                        dplyr::select(contains("pfas_t_")))[-13]

prs_name <- c("prs_uw", "prs_uw_tertile")


# Add PFAS quartiles, based on controls only-------
pfas_quartile <- data %>%
  filter(casetype == "Control") 
cuts_df <- data.frame()
for(pfas in pfas_name) {
  data_temp_control <- pfas_quartile %>% dplyr::select(pfas)
  cuts <- quantile(data_temp_control[1], 
                   c(0, 1/4,1/2,3/4, 1), na.rm = TRUE)
  cuts_df <- cuts_df %>% 
    bind_rows(data.frame(pfas_name = pfas,
                         Q0 = cuts[1],
                         Q1 = cuts[2],
                         Q2 = cuts[3],
                         Q3 = cuts[4],
                         Q4 = cuts[5]))
  data_temp <- data %>% dplyr::select(pfas)
  data_temp[2] <- case_when(data_temp[1] <= cuts[2] ~ 1, 
                            data_temp[1] <= cuts[3] ~ 2,
                            data_temp[1] <= cuts[4] ~ 3,
                            !is.na(data_temp[1]) ~ 4) |> as.factor() 
  colnames(data_temp)[2] <- paste0(pfas, "_quartile")
  data <- data %>% bind_cols(data_temp %>% dplyr::select(contains("quartile")))
}
cuts_df_l <- cuts_df %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 2)))%>%
  mutate(Quantile1 = paste0("≤ ", Q1 ),
         Quantile2 = paste0("(", Q1, ", ", Q2, "]"),
         Quantile3 = paste0("(", Q2, ", ", Q3, "]"),
         Quantile4 = paste0("> ", Q3))%>%
  select(-Q0, -Q1, -Q2, -Q3, -Q4) %>%
  pivot_longer(
               cols = c(Quantile1, Quantile2, Quantile3, Quantile4),
               names_to = "quartile",
               values_to = "value")

writexl::write_xlsx(cuts_df_l,fs::path(dir_result, "quartile_cuts_df.xlsx"))
# filter out controls/cases with no pfas values
# data <- data %>% filter(!setnum %in% c("1938", "1350", "1006")) # This part is for untargeted PFAS

# dataset split by outcome
sample_id1 <- data %>% filter(status == "HCC") %>% 
  dplyr::select(setnum)
sample_id2 <- data %>% filter(casetype == "NAFLD-related HCC") %>% 
  dplyr::select(setnum)
sample_id3 <- data %>% filter(casetype == "Other HCC") %>% 
  dplyr::select(setnum)


# Comparison of targeted PFAS and Untargeted PFAS

pfas_df_l <- data %>% select(sid, pfas_name) %>%
  pivot_longer(cols = pfas_name,
               values_to = "targeted_value")%>%
  drop_na() %>%
  mutate(name = str_remove_all(name,"pfas_t_"))%>%
  mutate(name = case_when(name == "NMFOSAA" ~ "NMeFOSAA",
                          name == "PFHPA" ~ "PFHpA",
                          name == "PFHPS" ~ "PFHpS",
                          name == "PFHXA" ~ "PFHxA",
                          name == "PFHXS" ~ "PFHxS",
                          name == "PFUNDA" ~ "PFUnDA",
                          TRUE ~ name))%>%
  filter(name != "NETFOSAA")%>% # since NETFOSAA was not in untargeted PFAS
  tidylog::left_join(data %>% select(sid,untargeted_pfas_name) %>%
                       pivot_longer(cols = untargeted_pfas_name,
                                    values_to = "untargeted_value")%>%
                       mutate(name = str_remove(name, "pfas_ut_") %>%
                                str_remove("Me_")),
                     by = c("sid", "name"))



ggplot(pfas_df_l, aes(x = untargeted_value, y = targeted_value)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~name, scales = "free")+
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")
# 
# ggsave(filename = fs::path(dir_figures,
#                            "Fig. correlation of targeted and untargeted PFAS.jpeg"),
#        width = 10,
#        height = 6,
#        dpi = 300,
#        bg = "white")

# Subset HCC,nafld_hcc, other hcc datasets
data_hcc <- data %>% filter(setnum %in% sample_id1$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~log2(.) %>% scale(.)%>%as.vector()) # Log transformed and Normalized PFAS
data_nafld_hcc <- data %>% filter(setnum %in% sample_id2$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~log2(.) %>% scale(.)%>% as.vector())
data_other_hcc <- data %>% filter(setnum %in% sample_id3$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~ log2(.) %>% scale(.)%>% as.vector())


data_hcc <- data_hcc %>%
  mutate(eth = as.character(eth)) %>%
  mutate(eth_new = ifelse(eth %in% c("White", "African American", "Native Hawaiian"), "Others", as.character(eth)))



# categorical pfas name
pfas_name_cat <- colnames(data %>% dplyr::select(contains("quartile")))


# covariates
# covars = c("smokestatus_imputed", "q1_d_chol", "q1_elih_imputed") 
covars = c("smokestatus_imputed", "diabetes", "q1_edih", "prs_wt_std", "v1", "v2", "v3", "v4", "v5","v6")
covars_matched <- c("sex", "eth", "studyarea","q1_age_cohent")


diet <- haven::read_sas(fs::path(dir_data %>% dirname() %>% dirname() %>% dirname(),
                                 "0_data_processing",
                                 "0_raw_data",
                                 "diet_PRS_data",
                                 "pfas_food_diet.sas7bdat")) %>%
  janitor::clean_names()


diet <- diet %>% mutate_at(.vars = c(colnames(diet)[-1]),
                           .funs = list(scld = ~scale(.)))

rm(data_temp, data_temp_control, pfas_df_l, 
   sample_id1, sample_id2, sample_id3, 
   pfas_quartile, cuts_df, cuts_df_l)


# standardizing the PRS by the control values in race/ethnicity

data_hcc <- data_hcc%>%
  group_by(eth) %>%
  mutate(
    prs_wt_std = (prs_wt - mean(prs_wt[casetype == "Control"], na.rm = TRUE)) /
      sd(prs_wt[casetype == "Control"], na.rm = TRUE)
  ) %>%
  ungroup()
