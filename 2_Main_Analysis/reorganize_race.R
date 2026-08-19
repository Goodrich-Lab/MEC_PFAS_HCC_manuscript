# Regenerate "Table S2. Dose response_0811_M.xlsx" by combining
#   Table S1. PFAS and HCC_081026.xlsx        (continuous PFAS-HCC associations)
#   Table S2. Dose response_081026.xlsx       (categorical/quartile dose-response)
#   Table S5. Heterogeneity test_081026.xlsx  (Cochran's Q heterogeneity test, by race)

source(fs::path(here::here("1_Project_Setup/!libraries.R")))
source(fs::path(here::here("1_Project_Setup/!directories.R")))

dir_supp     <- fs::path(dir_result, "Supplemental Tables")
dir_modified <- fs::path(dir_supp, "Tables_modified")

s1 <- readxl::read_xlsx(fs::path(dir_supp, "Table S1. PFAS and HCC_081026.xlsx"))
s2 <- readxl::read_xlsx(fs::path(dir_supp, "Table S2. Dose response_081026.xlsx"))
s5 <- readxl::read_xlsx(fs::path(dir_supp, "Table S5. Heterogeneity test_081026.xlsx"))

# Population group suffixes as they appear in the source column names (Table S1
# uses "\n" line breaks, Table S2 uses "\r\n"), paired with the "(n = ...)"
# labels used for these groups across the project's Supplemental Tables outputs
# (these labels are a fixed project convention and are not derivable from the
# case/control counts in the source column names).
pop_groups <- list(
  list(suffix_s1 = "Overall \n(n cases: 223\nn controls: 223)",
       suffix_s2 = "Overall \r\n(n cases: 223\r\nn controls: 223)",
       label     = "Overall \n(n = 464)"),
  list(suffix_s1 = "Japanese American \n(n cases: 88\nn controls: 88)",
       suffix_s2 = "Japanese American \r\n(n cases: 88\r\nn controls: 88)",
       label     = "Japanese American \n(n = 176)"),
  list(suffix_s1 = "Latino \n(n cases: 73\nn controls: 73)",
       suffix_s2 = "Latino \r\n(n cases: 73\r\nn controls: 73)",
       label     = "Latino \n(n = 158)"),
  list(suffix_s1 = "Others \n(n cases: 62\nn controls: 62)",
       suffix_s2 = "Others \r\n(n cases: 62\r\nn controls: 62)",
       label     = "Others\n(n = 124)")
)

pfas_dose_order <- unique(s2$pfas_name)
pfas_no_dose    <- setdiff(s1$pfas_name, pfas_dose_order)

# Heterogeneity test is stratified by race, matching the Overall/Japanese
# American/Latino/Others population breakdown used throughout this table
het <- s5 %>%
  dplyr::filter(type == "Race") %>%
  dplyr::select(pfas_name, Het_Q = `Q test Statistic`, Het_P = `P for heterogeneity`,
                Het_Qval = `Q for heterogeneity`)

term_levels <- c("__name__", "Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4",
                  "p_trend", "Continuous")

# Row skeleton (PFAS name + row type) shared by every population's OR/P/Q columns
skeleton_dose <- tidyr::expand_grid(pfas_name = pfas_dose_order,
                                     term = term_levels) %>%
  dplyr::mutate(pfas_name = factor(pfas_name, levels = pfas_dose_order),
                term = factor(term, levels = term_levels)) %>%
  dplyr::arrange(pfas_name, term)

skeleton_no_dose <- tibble::tibble(pfas_name = factor(pfas_no_dose, levels = pfas_no_dose),
                                    term = factor("Continuous", levels = term_levels))

skeleton <- dplyr::bind_rows(skeleton_dose, skeleton_no_dose)

# For a given population group, build the OR/P/Q columns in skeleton row order
build_pop_columns <- function(suffix_s1, suffix_s2) {

  or_col_s1 <- paste0("Odds Ratio[95%CI]_", suffix_s1)
  p_col_s1  <- paste0("P-Value_", suffix_s1)
  q_col_s1  <- paste0("Q-Value_", suffix_s1)

  or_col_s2 <- paste0("Odds Ratio[95%CI]_", suffix_s2)
  p_col_s2  <- paste0("P-Value_", suffix_s2)
  trend_col <- paste0("P trend_", suffix_s2)
  qval_col  <- paste0("Q value_", suffix_s2)

  dose <- s2 %>%
    dplyr::transmute(
      pfas_name = factor(pfas_name, levels = pfas_dose_order),
      term      = factor(term, levels = term_levels),
      OR        = ifelse(term == "Quartile 1", "1", .data[[or_col_s2]]),
      P         = .data[[p_col_s2]],
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

  cont_rows <- s1 %>%
    dplyr::filter(pfas_name %in% pfas_dose_order) %>%
    dplyr::transmute(pfas_name = factor(pfas_name, levels = pfas_dose_order),
                      term = factor("Continuous", levels = term_levels),
                      OR = .data[[or_col_s1]], P = .data[[p_col_s1]], Q = .data[[q_col_s1]])

  dose_block <- dplyr::bind_rows(name_rows, quartile_rows, trend_rows, cont_rows) %>%
    dplyr::arrange(pfas_name, term) %>%
    dplyr::select(OR, P, Q)

  no_dose_block <- s1 %>%
    dplyr::filter(pfas_name %in% pfas_no_dose) %>%
    dplyr::mutate(pfas_name = factor(pfas_name, levels = pfas_no_dose)) %>%
    dplyr::arrange(pfas_name) %>%
    dplyr::transmute(OR = .data[[or_col_s1]], P = .data[[p_col_s1]], Q = .data[[q_col_s1]])

  dplyr::bind_rows(dose_block, no_dose_block)
}

pop_columns <- lapply(pop_groups, function(g) build_pop_columns(g$suffix_s1, g$suffix_s2))

# The heterogeneity test is a PFAS-level result (not population-specific), so
# it only appears once per PFAS: on the Continuous row (PFBS_detected's single
# row uses "Continuous" as its skeleton term, so this covers it too)
het_columns <- skeleton %>%
  dplyr::mutate(pfas_name = as.character(pfas_name)) %>%
  dplyr::left_join(het, by = "pfas_name") %>%
  dplyr::transmute(
    Het_Q    = ifelse(term == "Continuous", Het_Q, NA_real_),
    Het_P    = ifelse(term == "Continuous", Het_P, NA_real_),
    Het_Qval = ifelse(term == "Continuous", Het_Qval, NA_real_)
  )

# The merged "PFAS" column shows the PFAS name on its own name row (and on
# PFBS_detected's single row, which has no separate name row), and the row
# type (Quartile 1-4 / p_trend / Continuous) everywhere else
pfas_col <- dplyr::case_when(
  as.character(skeleton$pfas_name) %in% pfas_no_dose ~ as.character(skeleton$pfas_name),
  skeleton$term == "__name__" ~ as.character(skeleton$pfas_name),
  TRUE ~ as.character(skeleton$term)
)

out <- dplyr::bind_cols(
  tibble::tibble(PFAS = pfas_col),
  pop_columns[[1]], tibble::tibble(sp1 = NA),
  pop_columns[[2]], tibble::tibble(sp2 = NA),
  pop_columns[[3]], tibble::tibble(sp3 = NA),
  pop_columns[[4]], tibble::tibble(sp4 = NA),
  het_columns,
  .name_repair = "unique_quiet"
)

# ---- Write out with a two-row header: population/heterogeneity group labels
# merged across their 3 columns, then the (group-agnostic) measure labels ----

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Sheet1")

group_header <- c("", pop_groups[[1]]$label, "", "", "",
                   pop_groups[[2]]$label, "", "", "",
                   pop_groups[[3]]$label, "", "", "",
                   pop_groups[[4]]$label, "", "", "",
                   "Test for Heterogeneity", "", "")

col_header <- c("PFAS",
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Odds Ratio[95%CI]", "P-Value", "Q-value", NA,
                 "Cochran's \U1D444 test statistic", "P-value", "Q-value")

openxlsx::writeData(wb, "Sheet1", t(group_header), startRow = 1, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", t(col_header), startRow = 2, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", out, startRow = 3, colNames = FALSE)

openxlsx::mergeCells(wb, "Sheet1", cols = 2:4,   rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 6:8,   rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 10:12, rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 14:16, rows = 1)
openxlsx::mergeCells(wb, "Sheet1", cols = 18:20, rows = 1)

fs::dir_create(dir_modified)
openxlsx::saveWorkbook(wb, fs::path(dir_modified, "Table S2. Dose response_M.xlsx"), overwrite = TRUE)
