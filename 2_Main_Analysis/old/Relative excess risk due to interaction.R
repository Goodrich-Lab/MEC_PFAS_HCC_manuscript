# PRS dichotomous interaction -------------
library(interactionR)

# Median split for PFHpA in cases only:
pfhpa_median <- quantile(data_hcc$pfas_t_PFHPA[data_hcc$status == 0], na.rm = TRUE, 
                         probs = 0.5)

pfhpa_median

data_hcc$pfas_t_PFHPA_median <- factor(
  ifelse(data_hcc$pfas_t_PFHPA <= pfhpa_median, 1, 2),
  levels = c(1, 2),
  labels = c("Low", "High")
)

table(data_hcc$pfas_t_PFHPA_median)[1]/table(data_hcc$pfas_t_PFHPA_median)[2]


mod1 <- clogit(status ~ pfas_t_PFHPA_median * prs_wt_median + 
                 smokestatus_imputed + 
                 diabetes + 
                 q1_edih + 
                 v1 + v2 + v3 + v4 + v5 + v6 + 
                 strata(setnum), 
               data = data_hcc,  
               method = "efron", 
               robust = TRUE)

summary(mod1)

epiR::epi.interaction(mod1, coef = c(1,2,13), param = "product")

out <- interactionR::interactionR(mod1, exposure_names = c("pfas_t_PFHPA_median", "prs_wt_median"), 
                           ci.type = "mover",
                           ci.level = 0.95,
                           em = F, # FALSE for interaction
                           recode = F)

out


# Plot data ----
plot_dat <- out$dframe %>%
  filter(Measures %in% c("OR00", "OR01", "OR10", "OR11")) %>%
  mutate(
    group = recode(
      Measures,
      "OR00" = "Low PFHPA / Low PRS",
      "OR01" = "Low PFHPA / High PRS",
      "OR10" = "High PFHPA / Low PRS",
      "OR11" = "High PFHPA / High PRS"
    ),
    group = factor(
      group,
      levels = c(
        "Low PFHPA / Low PRS",
        "Low PFHPA / High PRS",
        "High PFHPA / Low PRS",
        "High PFHPA / High PRS"
      )
    ),
    ci_label = if_else(
      Measures == "OR00",
      "Reference",
      sprintf("%.2f (%.2f–%.2f)", Estimates, CI.ll, CI.ul)
    ),
    CI.ll_plot = if_else(Measures == "OR00", Estimates, CI.ll),
    CI.ul_plot = if_else(Measures == "OR00", Estimates, CI.ul)
  )

ggplot(plot_dat, aes(x = group, y = Estimates)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbar(
    aes(ymin = CI.ll_plot, ymax = CI.ul_plot),
    width = 0,
    linewidth = 0.8
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(label = ci_label),
    hjust = -0.05,
    vjust = -0.9,
    size = 3.6,
    color = "black"
  ) +
  # scale_y_log10(
  #   breaks = c(0.5, 1, 2, 5, 10, 20),
  #   labels = label_number(accuracy = 0.1),
  #   limits = c(0.5, 25)
  # ) +
  coord_flip(clip = "off") +
  labs(
    x = NULL,
    y = "Odds ratio for HCC",
    # title = "Joint associations of PFHPA exposure and polygenic risk score",
    # subtitle = "Odds ratios relative to the low-PFHPA / low-PRS reference group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 95, 10, 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

