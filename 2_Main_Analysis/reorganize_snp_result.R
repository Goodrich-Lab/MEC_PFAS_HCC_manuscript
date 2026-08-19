# Regenerate "Table S8. PFAS and SNP_081026_M.xlsx" from
#   Table S8. PFAS and SNP_081026.xlsx
#
# Each population group has three effect types (PFAS Main Effect, SNP Main
# Effect, SNP * PFAS Interaction), each with its own OR/P/Q. This reorganizes
# that into a 3-row header: population label (row 1, merged across its
# columns), effect-type label (row 2, merged across its 3 measure columns),
# and plain Odds Ratio[95%CI]/P-Value/Q-value labels (row 3).

source(fs::path(here::here("1_Project_Setup/!libraries.R")))
source(fs::path(here::here("1_Project_Setup/!directories.R")))

dir_supp     <- fs::path(dir_result, "Supplemental Tables")
dir_modified <- fs::path(dir_supp, "Tables_modified")

s8 <- readxl::read_xlsx(fs::path(dir_supp, "Table S8. PFAS and SNP_081026.xlsx")) %>%
  dplyr::filter(!is.na(pfas_name))

populations <- c("Overall \n(n = 446)", "Japanese American \n(n = 176)", "Latino \n(n = 146)", "Others \n(n = 124)")
effects     <- c("PFAS Main Effect", "SNP Main Effect", "SNP * PFAS Interaction")

id_cols <- s8 %>% dplyr::select(pfas_name, snp, genes, rsid)

# Build the data columns, and record the column position of each block as we
# go so the header merges can be placed correctly afterwards
data_blocks   <- list()
pop_header    <- character()
effect_header <- character()
measure_header <- character()
pop_merge_ranges    <- list()
effect_merge_ranges <- list()

col <- ncol(id_cols)

for (pop in populations) {
  pop_start <- col + 1

  for (eff in effects) {
    or_col <- paste0("Odds Ratio[95%CI]_", eff, "_", pop)
    p_col  <- paste0("P-Value_", eff, "_", pop)
    q_col  <- paste0("Q-Value_", eff, "_", pop)

    block <- s8 %>% dplyr::transmute(OR = .data[[or_col]], P = .data[[p_col]], Q = .data[[q_col]])
    data_blocks[[length(data_blocks) + 1]] <- block

    eff_start <- col + 1
    col <- col + 3
    effect_merge_ranges[[length(effect_merge_ranges) + 1]] <- eff_start:col
    effect_header <- c(effect_header, eff)
    measure_header <- c(measure_header, "Odds Ratio[95%CI]", "P-Value", "Q-value")

    if (eff != effects[length(effects)]) {
      # spacer between effect-type blocks within this population. Only
      # measure_header (written sequentially as row 3) needs a filler entry
      # here -- effect_header/pop_header are indexed by merge range, not by
      # column position, so they must stay parallel to
      # effect_merge_ranges/pop_merge_ranges with no filler entries
      data_blocks[[length(data_blocks) + 1]] <- tibble::tibble(sp = NA)
      col <- col + 1
      measure_header <- c(measure_header, "")
    }
  }

  pop_merge_ranges[[length(pop_merge_ranges) + 1]] <- pop_start:col
  pop_header <- c(pop_header, pop)

  if (pop != populations[length(populations)]) {
    # spacer between population groups
    data_blocks[[length(data_blocks) + 1]] <- tibble::tibble(sp = NA)
    col <- col + 1
    measure_header <- c(measure_header, "")
  }
}

out <- dplyr::bind_cols(id_cols, dplyr::bind_cols(data_blocks, .name_repair = "unique_quiet"))

n_col <- ncol(out)
row1 <- rep("", n_col)
row2 <- rep("", n_col)
row3 <- c("pfas_name", "snp", "genes", "rsid", measure_header)
stopifnot(
  "row3 length must match the number of data columns" = length(row3) == n_col,
  "pop_header must be parallel to pop_merge_ranges" = length(pop_header) == length(pop_merge_ranges),
  "effect_header must be parallel to effect_merge_ranges" = length(effect_header) == length(effect_merge_ranges)
)

for (i in seq_along(pop_merge_ranges)) row1[pop_merge_ranges[[i]]] <- ""
for (i in seq_along(pop_merge_ranges)) row1[pop_merge_ranges[[i]][1]] <- pop_header[i]
for (i in seq_along(effect_merge_ranges)) row2[effect_merge_ranges[[i]][1]] <- effect_header[i]

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Sheet1")

openxlsx::writeData(wb, "Sheet1", t(row1), startRow = 1, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", t(row2), startRow = 2, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", t(row3), startRow = 3, colNames = FALSE)
openxlsx::writeData(wb, "Sheet1", out, startRow = 4, colNames = FALSE)

for (rng in pop_merge_ranges)    openxlsx::mergeCells(wb, "Sheet1", cols = rng, rows = 1)
for (rng in effect_merge_ranges) openxlsx::mergeCells(wb, "Sheet1", cols = rng, rows = 2)

fs::dir_create(dir_modified)
openxlsx::saveWorkbook(wb, fs::path(dir_modified, "Table S8. PFAS and SNP_081026_M.xlsx"), overwrite = TRUE)
