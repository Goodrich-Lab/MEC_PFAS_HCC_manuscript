# Regenerate an overall supplemental table for sex-stratified analysis by
# combining:
#   Table S4. PFAS HCC association by Sex_081126.xlsx  (continuous PFAS-HCC, by sex)
#   Table S4. Dose response_by_sex_081026.xlsx          (categorical/quartile dose-response, by sex)
#   Table S5. Heterogeneity test_081026.xlsx, type == "Sex"
#     (Cochran's Q test for heterogeneity between Female and Male)
#


source(fs::path(here::here("1_Project_Setup/!libraries.R")))
source(fs::path(here::here("1_Project_Setup/!directories.R")))
source(fs::path(here::here("1_Project_Setup/!functions.R")))

dir_supp     <- fs::path(dir_result, "Supplemental Tables")
dir_modified <- fs::path(dir_supp, "Tables_modified")

s4  <- readxl::read_xlsx(fs::path(dir_supp, "Table S4. PFAS HCC association by Sex_081126.xlsx"))
s2  <- readxl::read_xlsx(fs::path(dir_supp, "Table S4. Dose response_by_sex_081026.xlsx"))
s5  <- readxl::read_xlsx(fs::path(dir_supp, "Table S5. Heterogeneity test_081026.xlsx"))

# The continuous table (s4) labels sex "female"/"male"; the dose-response
# table (s2) labels it "Female"/"Male". The sample sizes below are a fixed
# label (as with the "(n = ...)" labels in reorganize.R, not derivable from
# the files used here) provided directly by the study team.
groups <- list(
  list(cont_key = "female", dose_suffix = "Female", label = "Female\n(n = 138)"),
  list(cont_key = "male",   dose_suffix = "Male",   label = "Male\n(n = 308)")
)

pfas_dose_order <- unique(s2$pfas_name)
pfas_no_dose    <- setdiff(s4$pfas_name, pfas_dose_order)

het <- s5 %>%
  dplyr::filter(type == "Sex") %>%
  dplyr::select(pfas_name, Het_Q = `Q test Statistic`, Het_P = `P for heterogeneity`,
                Het_Qval = `Q for heterogeneity`)

term_levels <- c("__name__", "Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4",
                  "p_trend", "Continuous")

skeleton_dose <- tidyr::expand_grid(pfas_name = pfas_dose_order,
                                     term = term_levels) %>%
  dplyr::mutate(pfas_name = factor(pfas_name, levels = pfas_dose_order),
                term = factor(term, levels = term_levels)) %>%
  dplyr::arrange(pfas_name, term)

skeleton_no_dose <- tibble::tibble(pfas_name = factor(pfas_no_dose, levels = pfas_no_dose),
                                    term = factor("Continuous", levels = term_levels))

skeleton <- dplyr::bind_rows(skeleton_dose, skeleton_no_dose)

# For a given group, build the OR/P/Q columns in skeleton row order
build_group_columns <- function(cont_key, dose_suffix) {

  or_col    <- paste0("Odds Ratio[95%CI]_", dose_suffix)
  p_col     <- paste0("P-Value_", dose_suffix)
  trend_col <- paste0("P trend_", dose_suffix)
  qval_col  <- paste0("Q value_", dose_suffix)

  dose <- s2 %>%
    dplyr::transmute(
      pfas_name = factor(pfas_name, levels = pfas_dose_order),
      term      = factor(term, levels = term_levels),
      OR        = ifelse(term == "Quartile 1", "1", .data[[or_col]]),
      P         = .data[[p_col]],
      p_trend   = .data[[trend_col]],
      Q_trend   = .data[[qval_col]]
    )

  name_rows <- tibble::tibble(pfas_name = factor(pfas_dose_order, levels = pfas_dose_order),
                               term = factor("__name__", levels = term_levels),
                               OR = NA_character_, P = NA_real_, Q = NA_real_)

  quartile_rows <- dose %>%
    dplyr::transmute(pfas_name, term, OR, P, Q = NA_real_)

  trend_rows <- dose %>%
    dplyr::distinct(pfas_name, p_trend, Q_trend) %>%
    dplyr::transmute(pfas_name, term = factor("p_trend", levels = term_levels),
                      OR = NA_character_, P = p_trend, Q = Q_trend)

  cont_rows <- s4 %>%
    dplyr::filter(pfas_name %in% pfas_dose_order, group == cont_key) %>%
    dplyr::transmute(pfas_name = factor(pfas_name, levels = pfas_dose_order),
                      term = factor("Continuous", levels = term_levels),
                      OR = `Odds Ratio[95%CI]`, P = `P-Value`, Q = `Q-Value`)

  dose_block <- dplyr::bind_rows(name_rows, quartile_rows, trend_rows, cont_rows) %>%
    dplyr::arrange(pfas_name, term) %>%
    dplyr::select(OR, P, Q)

  no_dose_block <- s4 %>%
    dplyr::filter(pfas_name %in% pfas_no_dose, group == cont_key) %>%
    dplyr::mutate(pfas_name = factor(pfas_name, levels = pfas_no_dose)) %>%
    dplyr::arrange(pfas_name) %>%
    dplyr::transmute(OR = `Odds Ratio[95%CI]`, P = `P-Value`, Q = `Q-Value`)

  dplyr::bind_rows(dose_block, no_dose_block)
}

group_columns <- lapply(groups, function(g) build_group_columns(g$cont_key, g$dose_suffix))

# Heterogeneity test is a PFAS-level result, shown once per PFAS on the
# Continuous row (PFBS_detected's single row uses "Continuous" as its
# skeleton term, so this covers it too)
het_columns <- skeleton %>%
  dplyr::mutate(pfas_name = as.character(pfas_name)) %>%
  dplyr::left_join(het, by = "pfas_name") %>%
  dplyr::transmute(
    Het_Q    = ifelse(term == "Continuous", Het_Q, NA_real_),
    Het_P    = ifelse(term == "Continuous", Het_P, NA_real_),
    Het_Qval = ifelse(term == "Continuous", Het_Qval, NA_real_)
  )

pfas_col <- dplyr::case_when(
  as.character(skeleton$pfas_name) %in% pfas_no_dose ~ as.character(skeleton$pfas_name),
  skeleton$term == "__name__" ~ as.character(skeleton$pfas_name),
  TRUE ~ as.character(skeleton$term)
)

out <- dplyr::bind_cols(
  tibble::tibble(PFAS = pfas_col),
  group_columns[[1]], tibble::tibble(sp1 = NA),
  group_columns[[2]], tibble::tibble(sp2 = NA),
  het_columns,
  .name_repair = "unique_quiet"
)

# ---- Write out with a two-row header: group/heterogeneity labels merged
# across their 3 columns, then the (group-agnostic) measure labels ----

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Sheet1")

group_header <- c("", groups[[1]]$label, "", "", "",
                   groups[[2]]$label, "", "", "",
                   "Test for Heterogeneity", "", "")

col_header <- c("PFAS",
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Cochran's \U1D444 test statistic", "P-value", "Q-value")

openxlsx::writeData(wb, "Sheet1", t(group_header), startRow = 1, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", t(col_header), startRow = 2, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", out, startRow = 3, colNames = FALSE)

openxlsx::mergeCells(wb, "Sheet1", cols = 2:4,  rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 6:8,  rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 10:12, rows = 1)

fs::dir_create(dir_modified)
openxlsx::saveWorkbook(wb, fs::path(dir_modified, "Table S4. Dose response_by_sex_M.xlsx"), overwrite = TRUE)
