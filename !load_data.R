# Set up data for analysis
# library(tidylog)


# Read Data ---------------------------------------------------------------
data <- readRDS(fs::path(
  dir_cleaned_data,
  "hhear_complete_data_v1.RDS")) 

table(is.na(data$bmi))

data |> 
  select(-diagnosis_age, -year_between, 
         -contains("chem_"), -contains("met_"), -contains("pfas_"), 
         -c(x1_220796686:prs_uw_tertile)) |>
  naniar::gg_miss_upset(nintersects = 20, nsets = 15)

# calorie <- haven::read_sas(fs::path(dir_data %>% dirname() %>% dirname() %>% dirname(),
#                                     "0_data_processing",
#                                     "0_raw_data",
#                                     "diet_PRS_data",
#                                     "pfas_calories_072624.sas7bdat")) %>%
#   janitor::clean_names()
# # 
# data <- data %>% tidylog::left_join(calorie, by = "barcode")%>%
#   mutate(lg_q1_calories = log(q1_calories))
# 
# write_rds(data, fs::path(
#   dir_cleaned_data,"hhear_complete_data_v1.RDS" ))
# Adding 


snp_info <- readxl::read_xlsx(fs::path(dir_data %>% dirname() %>% dirname() %>% dirname(),
                                       "0_data_processing",
                                       "0_raw_data",
                                       "diet_PRS_data",
                                       "snp_info.xlsx")) %>%
  janitor::clean_names() %>%
  mutate(snp_name = paste0("x", chr, "_", pos38))


untargeted_met_data <- readRDS(fs::path(dir_cleaned_data, "Untargeted_feature_metadata.RDS"))

# variables--------------
pfas_name <- colnames(data %>% 
                        dplyr::select(contains("pfas")) %>% 
                        dplyr::select(contains("lg2")) %>% 
                        dplyr::select(-contains("quartile")))

prs_name <- c("prs_uw", "prs_uw_tertile")

data <- data %>% 
  mutate_at(.vars = c(pfas_name),
            .funs = list(raw = ~exp(.)%>% as.vector(.)))%>%
  rename_at(.vars = str_c(pfas_name,"_","raw"),
            .funs = ~str_remove(., "_lg2"))

pfas_name_raw <- colnames(data%>% dplyr::select(contains("raw")))

# Impute missing smoking status ----------------------------------------

# Select only the variables involved
impute_data <- data[, c("smokestatus", "casetype", "sex",  "q1_byr", "q1_elih", "bmi", "eth")] 

impute_data <- as.data.frame(impute_data) 

for (i in 1:ncol(impute_data)) {
  haven::zap_labels(impute_data[,i])
}

# Convert categorical variables to factors
impute_data$smokestatus <- as.factor(as.character(impute_data$smokestatus))
impute_data$casetype    <- as.factor(as.character(impute_data$casetype))
impute_data$sex         <- as.factor(as.character(impute_data$sex))
impute_data$q1_byr      <- as.numeric(impute_data$q1_byr)
impute_data$bmi         <- as.numeric(impute_data$bmi)

# Run imputation
set.seed(1234)
imputed_data <- missForest::missForest(impute_data, maxiter = 10, ntree = 100)


# Get the completed data
completed_data <- imputed_data$ximp

# Replace the original smokestatus with the imputed values
data$smokestatus_imputed <- completed_data$smokestatus
data$q1_elih_imputed <- completed_data$q1_elih
data$bmi_imputed <- completed_data$bmi

data <- data |> select(-c(smokestatus, q1_elih, bmi, bmi_cat))


# Add PFAS quartiles, based on controls only-------
pfas_quartile <- data %>%
  filter(casetype == "Control") 

for(pfas in pfas_name) {
  data_temp_control <- pfas_quartile %>% dplyr::select(pfas)
  cuts <- quantile(data_temp_control[1], 
                   c(1/4,1/2,3/4), na.rm = TRUE)
  data_temp <- data %>% dplyr::select(pfas)
  data_temp[2] <- case_when(data_temp[1] <= cuts[1] ~ 1, 
                            data_temp[1] <= cuts[2] ~ 2,
                            data_temp[1] <= cuts[3] ~ 3,
                            !is.na(data_temp[1]) ~ 4) |> as.factor() 
  colnames(data_temp)[2] <- paste0(pfas, "_quartile")
  data <- data %>% bind_cols(data_temp %>% dplyr::select(contains("quartile")))
}

# filter out controls/cases with no pfas values
data <- data %>% filter(!setnum %in% c("1938", "1350", "1006"))

# dataset split by outcome
sample_id1 <- data %>% filter(status == "HCC") %>% 
  dplyr::select(setnum)
sample_id2 <- data %>% filter(casetype == "NAFLD-related HCC") %>% 
  dplyr::select(setnum)
sample_id3 <- data %>% filter(casetype == "Other HCC") %>% 
  dplyr::select(setnum)

data_hcc <- data %>% filter(setnum %in% sample_id1$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~scale(.)%>%as.vector())
data_nafld_hcc <- data %>% filter(setnum %in% sample_id2$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~scale(.)%>% as.vector())
data_other_hcc <- data %>% filter(setnum %in% sample_id3$setnum) %>%
  mutate(status = ifelse(status == "HCC", 1, 0)) %>%
  mutate_at(.vars = pfas_name,
            .funs = ~scale(.)%>% as.vector())


data_hcc <- data_hcc %>%
  mutate(eth = as.character(eth)) %>%
  mutate(eth_new = ifelse(eth %in% c("White", "African American", "Native Hawaiian"), "Others", as.character(eth)))



# categorical pfas name
pfas_name_cat <- colnames(data %>% dplyr::select(contains("quartile")))


# covariates
covars = c("smokestatus_imputed", "alcohol_intake", "q1_d_chol", "q1_elih_imputed") 
covars_matched <- c("sex", "eth", "studyarea","q1_age_cohent")


diet <- haven::read_sas(fs::path(dir_data %>% dirname() %>% dirname() %>% dirname(),
                                 "0_data_processing",
                                 "0_raw_data",
                                 "diet_PRS_data",
                                 "pfas_food_diet.sas7bdat")) %>%
  janitor::clean_names()


diet <- diet %>% mutate_at(.vars = c(colnames(diet)[-1]),
                           .funs = list(scld = ~scale(.)))

rm(data_temp, data_temp_control)