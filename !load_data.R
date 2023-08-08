# Set up data for analysis
# library(tidylog)

# Read Data ---------------------------------------------------------------

data <- readRDS(fs::path(dir_cleaned_data, "p27_hhear_pfas_covar_imputed_V2.RDS"))

targeted_met_data <- readRDS(fs::path(dir_cleaned_data, "Targeted_feature_tables_filtered_lg2.RDS"))

untargeted_met_data <- readRDS(fs::path(dir_cleaned_data, "Untargeted_feature_tables_filtered_lg2.RDS"))
