rename_pfas <- function(df){
  df<- df %>% 
    mutate(exposure = case_when(grepl("pfda", exposure)~"PFDA",
                               grepl("pfos", exposure)~"PFOS",
                               grepl("pfhpa", exposure)~"PFHpA",
                               grepl("pfhxs", exposure)~"PFHxS",
                               grepl("pfhps", exposure)~"PFHpS",
                               grepl("pfna", exposure)~"PFNA",
                               grepl("pfoa", exposure)~"PFOA",
                               grepl("pfuda", exposure)~ "PFUdA",
                               grepl("pfds", exposure)~ "PFDS",
                               grepl("pfpes", exposure)~ "PFPeS",
                               grepl("pftrda", exposure)~ "PFTrDA"
                            ))
  return(df)
}

