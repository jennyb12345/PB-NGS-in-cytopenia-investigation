# ============================================================
# Table 4 PPV of recurrent genes (>=6 detections) for Myeloid outcome
# with 95% CI
# Exclude VUS and VCS-germline only
# ============================================================

# 1) Install/load required packages
required_packages <- c("readxl", "dplyr", "stringr", "tidyr", "purrr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# 2) File + sheet settings
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
clinical_tab <- "Clinic information + FBC"
pb_tab       <- "PB NGS"
id_col       <- "Patient ID"  

# 3) Read sheets
clinical_df <- read_excel(file_path, sheet = clinical_tab)
pb_df       <- read_excel(file_path, sheet = pb_tab)

# 4) Clean column names
names(clinical_df) <- str_squish(names(clinical_df))
names(pb_df)       <- str_squish(names(pb_df))

# 5) Clinical outcome dataset (patient-level)
clin_outcome <- clinical_df %>%
  transmute(
    patient_id = .data[[id_col]],
    final_dx = as.character(`Final Dx Classification`)
  ) %>%
  mutate(
    myeloid = if_else(str_squish(final_dx) == "Myeloid", 1L, 0L)
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

# 6) Detect slot columns in PB NGS
gene_cols <- names(pb_df)[str_detect(names(pb_df), regex("^Gene\\s*\\d+$", ignore_case = TRUE))]
var_cols  <- names(pb_df)[str_detect(names(pb_df), regex("^Variant\\s*\\d+$", ignore_case = TRUE))]
cls_cols  <- names(pb_df)[str_detect(names(pb_df), regex("^Classification\\s*\\d+$", ignore_case = TRUE))]

if (length(gene_cols) == 0 || length(var_cols) == 0 || length(cls_cols) == 0) {
  stop("Could not detect Gene/Variant/Classification slot columns (e.g., Gene 1, Variant 1, Classification 1).")
}

# 7) Sample-level exclusion in PB NGS
pb_filt <- pb_df %>%
  mutate(
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

# 8) Long format and slot-align Gene N + Variant N + Classification N
pb_long_gene <- pb_filt %>%
  select(all_of(c(id_col, gene_cols))) %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "gene_col", values_to = "gene") %>%
  mutate(slot = str_extract(gene_col, "\\d+"))

pb_long_var <- pb_filt %>%
  select(all_of(c(id_col, var_cols))) %>%
  pivot_longer(cols = all_of(var_cols), names_to = "var_col", values_to = "variant") %>%
  mutate(slot = str_extract(var_col, "\\d+")) %>%
  select(all_of(id_col), slot, variant)

pb_long_cls <- pb_filt %>%
  select(all_of(c(id_col, cls_cols))) %>%
  pivot_longer(cols = all_of(cls_cols), names_to = "cls_col", values_to = "classification") %>%
  mutate(
    slot = str_extract(cls_col, "\\d+"),
    classification = str_to_upper(str_squish(as.character(classification)))
  ) %>%
  select(all_of(id_col), slot, classification)

pb_long <- pb_long_gene %>%
  left_join(pb_long_var, by = c(id_col, "slot")) %>%
  left_join(pb_long_cls, by = c(id_col, "slot")) %>%
  transmute(
    patient_id = .data[[id_col]],
    gene = str_to_upper(str_squish(as.character(gene))),
    variant = str_squish(as.character(variant)),
    classification = classification
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(variant), variant != "") %>%
  # Keep only Path variants
  filter(classification == "PATH")

# 9) Patient-gene presence (one row per patient per gene)
patient_gene <- pb_long %>%
  distinct(patient_id, gene)

# 10) Identify recurrent genes (>=6 unique patients)
recurrent_genes <- patient_gene %>%
  count(gene, name = "n_detected", sort = TRUE) %>%
  filter(n_detected >= 6)

if (nrow(recurrent_genes) == 0) {
  stop("No recurrent genes found with detection count >= 6 after exclusions/filters.")
}

# 11) Build patient universe for PPV calculations
patients <- clin_outcome %>%
  select(patient_id, myeloid)

# 12) Function: PPV + exact 95% CI per gene
calc_ppv <- function(g) {
  # Patients positive for this gene
  pos_ids <- patient_gene %>%
    filter(gene == g) %>%
    distinct(patient_id)
  
  # Merge onto all patients
  df <- patients %>%
    left_join(pos_ids %>% mutate(test_pos = 1L), by = "patient_id") %>%
    mutate(test_pos = if_else(is.na(test_pos), 0L, test_pos))
  
  TP <- sum(df$test_pos == 1 & df$myeloid == 1, na.rm = TRUE)
  FP <- sum(df$test_pos == 1 & df$myeloid == 0, na.rm = TRUE)
  
  n_test_pos <- TP + FP
  
  # PPV = TP / (TP + FP)
  ppv <- ifelse(n_test_pos > 0, TP / n_test_pos, NA_real_)
  
  # Exact binomial CI for PPV among test positives
  ci <- if (n_test_pos > 0) binom.test(TP, n_test_pos, conf.level = 0.95)$conf.int else c(NA, NA)
  
  tibble(
    gene = g,
    n_detected = n_test_pos,
    TP_myeloid = TP,
    FP_non_myeloid = FP,
    PPV = ppv,
    PPV_95CI_low = ci[1],
    PPV_95CI_high = ci[2]
  )
}

# 13) Compute PPV table for all recurrent genes
ppv_results <- map_dfr(recurrent_genes$gene, calc_ppv) %>%
  left_join(recurrent_genes, by = "gene", suffix = c("", "_check")) %>%
  mutate(
    PPV_percent = round(100 * PPV, 1),
    PPV_95CI = paste0(
      round(100 * PPV_95CI_low, 1), "% to ",
      round(100 * PPV_95CI_high, 1), "%"
    )
  ) %>%
  arrange(desc(PPV), desc(n_detected))

# 14) Print results
print(ppv_results)

# 15) Export
write.csv(ppv_results, "PPV_recurrent_genes_myeloid_with95CI.csv", row.names = FALSE)
write.csv(recurrent_genes, "Recurrent_genes_detected_ge4.csv", row.names = FALSE)

# ============================================================
# 16) Additional PPV analysis by gene group
# ============================================================

# Define gene groups
gene_group_map <- tibble::tibble(
  gene = c("SRSF2","SF3B1","U2AF1",
           "TET2","ASXL1","DNMT3A","IDH2",
           "CBL","KRAS",
           "CUX1","TP53","STAG2"),
  gene_group = c(rep("Spliceosome", 3),
                 rep("Epigenetic", 4),
                 rep("Signalling", 2),
                 rep("Other", 3))
)

# Keep only patient-gene calls that are in grouped genes
patient_gene_grouped <- patient_gene %>%
  inner_join(gene_group_map, by = "gene")

# Collapse to patient-group presence (patient positive for a group if any gene in that group detected)
patient_group <- patient_gene_grouped %>%
  distinct(patient_id, gene_group)

# Function to compute PPV for a gene group
calc_ppv_group <- function(ggrp) {
  pos_ids <- patient_group %>%
    filter(gene_group == ggrp) %>%
    distinct(patient_id)
  
  df <- patients %>%
    left_join(pos_ids %>% mutate(test_pos = 1L), by = "patient_id") %>%
    mutate(test_pos = if_else(is.na(test_pos), 0L, test_pos))
  
  TP <- sum(df$test_pos == 1 & df$myeloid == 1, na.rm = TRUE)
  FP <- sum(df$test_pos == 1 & df$myeloid == 0, na.rm = TRUE)
  n_test_pos <- TP + FP
  
  ppv <- ifelse(n_test_pos > 0, TP / n_test_pos, NA_real_)
  ci <- if (n_test_pos > 0) binom.test(TP, n_test_pos, conf.level = 0.95)$conf.int else c(NA, NA)
  
  tibble::tibble(
    analysis_level = "Gene_group",
    item = ggrp,
    n_detected = n_test_pos,
    TP_myeloid = TP,
    FP_non_myeloid = FP,
    PPV = ppv,
    PPV_95CI_low = ci[1],
    PPV_95CI_high = ci[2],
    PPV_percent = round(100 * ppv, 1),
    PPV_95CI = paste0(round(100 * ci[1], 1), "% to ", round(100 * ci[2], 1), "%")
  )
}

# Run grouped PPV
ppv_group_results <- purrr::map_dfr(unique(gene_group_map$gene_group), calc_ppv_group)

# ============================================================
# 17) Combine per-gene + gene-group into one CSV
# ============================================================

# Harmonize per-gene table to same columns
ppv_genes_export <- ppv_results %>%
  transmute(
    analysis_level = "Gene",
    item = gene,
    n_detected = n_detected,
    TP_myeloid = TP_myeloid,
    FP_non_myeloid = FP_non_myeloid,
    PPV = PPV,
    PPV_95CI_low = PPV_95CI_low,
    PPV_95CI_high = PPV_95CI_high,
    PPV_percent = PPV_percent,
    PPV_95CI = PPV_95CI
  )

ppv_combined <- bind_rows(ppv_genes_export, ppv_group_results) %>%
  arrange(analysis_level, desc(PPV), desc(n_detected))

# Save one combined CSV
write.csv(ppv_combined, "PPV_myeloid_genes_and_groups_with95CI_ge5.csv", row.names = FALSE)

cat("\nDone.\nAdditional file created:\n- PPV_myeloid_genes_and_groups_with95CI_ge5.csv\n")
