# Rename PFAS -----------------------
rename_pfas_new <- function(pfas_names){
  x <- tibble(pfas = pfas_names)
  pfas2 <-  x %>%
    mutate(pfas2 = str_split_fixed(pfas, "pfas_", 2)[,2]) %>%
    mutate(pfas2 = str_split_fixed(pfas2, "_raw", 2)[,1]) %>%
    mutate(pfas2 = case_when(pfas2 == "4H_PFBA" ~ "4H-PFBA",
                                        pfas2 == "10_2_FToH" ~ "10:2 FTOH",
                                        pfas2 == "2_aminohexafluoropropan_2_ol" ~ "2-Aminohexafluoropropan-2-ol",
                                        pfas2 == "Me_PFHxA" ~ "MePFHxA",
                                        TRUE ~ pfas2) %>% 
             as.factor() %>% 
             fct_relevel(., 
                    "PFHpA", 
                    "4H-PFBA", 
                    "PFBS", 
                    "NMeFOSAA", 
                    "PFDA",
                    "10:2 FTOH", 
                    "PFPrA",
                    "PFDoS",
                    "PFNA",
                    "PFHxS",
                    "PFUnDA",
                    "PFOS",
                    "POSF",
                    "PFOA", 
                    "PFHpS",
                    "2-Aminohexafluoropropan-2-ol",
                    "MePFHxA"))
  
  return(pfas2$pfas2)
}
