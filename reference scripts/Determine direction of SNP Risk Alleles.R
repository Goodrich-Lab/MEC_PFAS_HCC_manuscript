# Determining SNP Risk Alleles
source(fs::path(here::here("!libraries.R")))
source(fs::path(here::here("!directories.R")))
# source(fs::path(here::here("functions.R")))
# source(fs::path(here::here("!load_data.R")))

# Read in data ------

# Read in original SNP data
og_snp_data <- read_rds(fs::path(dir_home  %>% dirname() %>% dirname(), 
                                 '0_data_processing', '0_raw_data', 
                                 'diet_PRS_data', 'pfas_prs.RDS'))

write_csv(og_snp_data, "og_snp_data.csv")


# Read in PRS scoring data 
prs_scoring_data <- readxl::read_xlsx(fs::path(dir_home  %>% dirname() %>% dirname(), 
                                               '0_data_processing', '0_raw_data', 
                                               'diet_PRS_data', 'MEC PRS Scoring File.xlsx'), 
                                      skip = 1) |> 
  janitor::clean_names() |>
  dplyr::select(-c(1:7))

# Calculate log odds ratio
prs_scoring_data$log_or = log(prs_scoring_data$overall_or)

# 2. Clean Data -------------
# 2.a. Clean and merge data frames --------------

# Merge with snp_info DF
# Should be 21 matches
prs_scoring <- snp_info |> 
  tidylog::left_join(prs_scoring_data) |> 
  tidylog::select(snp_name, overall_or, everything())

# GWAS Data
mec_snps <- og_snp_data |> 
  column_to_rownames('barcode') |> 
  dplyr::select(x1_220796686:x22_43928850) |> 
  as.matrix()

# Check that the order of cols is the same
table(prs_scoring$snp_name == colnames(mec_snps))


# 2.b. Calculate PRS based on original values ------
og_snp_data$prs_calculated <- as.numeric(mec_snps %*% prs_scoring$overall_or)

# Clearly, there are issues (should be a 1:1 line): 
ggplot(og_snp_data, aes(x = prs_calculated, 
                        y = prs_uw)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) 


# 3. Figure out which are risk alleles ---------
## 3.a. Start with comparing risk allele freq from PRS vs. in our data ----
# Calculate summary of snps
snp_summaries_cont <- data %>%
  tidylog::filter(status == "Control") %>%
  dplyr::select(eth, x1_220796686:x22_43928850) %>%
  pivot_longer(cols = x1_220796686:x22_43928850) %>% 
  group_by(name, eth) %>%  # , eth
  summarise(#pct_0           = mean(value == 0, na.rm = TRUE),
    # pct_base_allele = sum(value < 0.5, na.rm = TRUE)/length(value),
    #variance        = var(value, na.rm = TRUE),
    #mean            = mean(value, na.rm = TRUE), 
    #median          = median(value, na.rm = TRUE)
  )

snp_summaries_cont_w <- snp_summaries_cont |> 
  pivot_wider(names_from = eth, values_from = pct_base_allele)

# Risk alleles per ethnicity
risk_allele_freq_all_of_mec <- prs_scoring_data |> 
  filter(included_in_mec_prs == 1) |> 
  tidylog::select(gene, contains("risk_allele_frequency")) 

# Combine summary with 
snp_summaries_cont2 <-  tidylog::left_join(snp_summaries_cont_w, 
                                           snp_info, 
                                           by = c("name" = "snp_name")) |>
  tidylog::left_join(risk_allele_freq_all_of_mec, 
                     by = c("gene" = "gene")) |> 
  # filter(pct_0 < 0.75) %>% #dplyr::select(name)
  dplyr::select(name, gene, risk_allele_frequency_mec, everything()) |>
  janitor::clean_names()

# Calculate summary by ethnicity
snp_summaries_cont2 <- snp_summaries_cont2 |> 
  mutate(
    white_delta  = white             - white_risk_allele_frequency,                      
    aa_delta     = african_american  - aa_risk_allele_frequency,                  
    ja_delta     = japanese_american - ja_risk_allele_frequency,                       
    latino_delta = latino            - latino_risk_allele_frequency,                 
    nh_delta     = native_hawaiian   - nh_risk_allele_frequency) |> 
  rowwise() |> 
  mutate(
    average_delta = sum(white_delta, aa_delta,ja_delta,latino_delta, nh_delta)/5)

deltas <- snp_summaries_cont2 |> 
  select(gene, average_delta, 
         c(white_delta, white, white_risk_allele_frequency,                      
           aa_delta, african_american, aa_risk_allele_frequency,                  
           ja_delta, japanese_american, ja_risk_allele_frequency,                       
           latino_delta, latino, latino_risk_allele_frequency,                 
           nh_delta, native_hawaiian, nh_risk_allele_frequency))

# Plot
ggplot(snp_summaries_cont, aes(x = risk_allele_frequency_mec, 
                               y = pct_base_allele, 
                               label = gene)) + 
  geom_point() + 
  geom_label() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylim(c(0,1)) + xlim(c(0,1)) + 
  facet_wrap(~eth) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5)



## 3.b. Calculate rescaled variables for these SNPS: ----
rev_snps <- prs_scoring |> 
  tidylog::filter(rev_order_in_data == 1)

og_snp_data_rescaled <- data |> 
  mutate_at(.vars = c(rev_snps$snp_name), 
            .fun = function(x){(-1*x) + 2})



# Check that the order of cols is the same
table(prs_scoring$snp_name == colnames(mec_snps))


# 2.b. Calculate PRS based on new values ------
# GWAS Data
mec_snps_rescaled <- og_snp_data_rescaled |> 
  column_to_rownames('barcode') |> 
  dplyr::select(x1_220796686:x22_43928850) |> 
  as.matrix()

# Calculate PRS:
og_snp_data$prs_calculated_new <- as.numeric(mec_snps_rescaled %*% prs_scoring$overall_or)


# Clearly, there are issues (should be a 1:1 line): 
base::plot(x = og_snp_data$prs_calculated_new, y = og_snp_data$prs_uw)




# 2.c. Check if differences are from eth by calculating with eth specific values ------
og_snp_data_rescaled_ja <- og_snp_data_rescaled |> 
  tidylog::filter(eth == 'Japanese American')

mec_snps_rescaled <- og_snp_data_rescaled_ja |>
  column_to_rownames('barcode') |> 
  dplyr::select(x1_220796686:x22_43928850) |> 
  as.matrix()

# Calculate PRS:
og_snp_data_rescaled_ja$prs_calculated_new <- as.numeric(mec_snps_rescaled %*% prs_scoring$ja_or)


# Clearly, there are issues (should be a 1:1 line): 
ggplot(og_snp_data_rescaled_ja, aes(x = prs_calculated_new, 
                                    y = prs_uw)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) 
# ylim(c(0,1)) + xlim(c(0,1)) + 
