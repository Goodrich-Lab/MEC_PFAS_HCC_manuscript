# directory fs::paths for file architecture
library(tidyverse)
# home directory for project
dir_home <- here::here() %>% dirname() %>% dirname()

# raw data folder
dir_data <- fs::path(dir_home, "0_data")

# programs folder
dir_cleaned_data <- fs::path(dir_home %>% 
  dirname() %>% dirname() , "1_analysis_ready_data")

# result folder
dir_result <- fs::path(dir_home, "2_results")

# figure folder
dir_figures <- fs::path(dir_home, "3_figures")
