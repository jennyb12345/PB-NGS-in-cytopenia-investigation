# ============================================================
# Table 5- Recurrent specific variants (>=3) + PPV (95% CI) for Myeloid
# Exclude VUS and VCS germline only= "path"/VCS variants only.
# ============================================================

# 1) Packages
required_packages <- c("readxl", "dplyr", "stringr", "tidyr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

# 2) File + tabs
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
pb_tab       <- "PB NGS"
clinical_tab <- "Clinic information + FBC"
id_col       <- "Patient ID"  

# If path issue:
# file_path <- file.choose()

# 3) Read sheets
pb <- read_excel(file_path, sheet = pb_tab)
cl <- read_excel(file_path, sheet = clinical_tab)

# 4) Clean names
names(pb) <- str_squish(names(pb))
names(cl) <- str_squish(names(cl))

# 5) Detect slot columns
gene_cols <- names(pb)[str_detect(names(pb), regex("^Gene\\s*\\d+$", ignore_case = TRUE))]
var_cols  <- names(pb)[str_detect(names(pb), regex("^Variant\\s*\\d+$", ignore_case = TRUE))]
cls_cols  <- names(pb)[str_detect(names(pb), regex("^Classification\\s*\\d+$", ignore_case = TRUE))]

if (length(gene_cols) == 0 || length(var_cols) == 0 || length(cls_cols) == 0) {
  stop("Could not detect Gene/Variant/Classification columns (e.g., Gene 1, Variant 1, Classification 1).")
}

# 6) Clinical outcome table (patient-level)
clin_outcome <- cl %>%
  transmute(
    patient_id = .data[[id_col]],
    final_dx = as.character(`Final Dx Classification`)
  ) %>%
  mutate(myeloid = if_else(str_squish(final_dx) == "Myeloid", 1L, 0L)) %>%
  distinct(patient_id, .keep_all = TRUE)

# 7) PB sample-level exclusions
pb_filt <- pb %>%
  mutate(
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

# 8) Long format: align Gene N + Variant N + Classification N
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

# 9) Keep only path variants and build specific variant label
specific_variants <- pb_long_gene %>%
  left_join(pb_long_var, by = c(id_col, "slot")) %>%
  left_join(pb_long_cls, by = c(id_col, "slot")) %>%
  transmute(
    patient_id = .data[[id_col]],
    gene = str_to_upper(str_squish(as.character(gene))),
    variant = str_squish(as.character(variant)),
    classification = classification
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(variant), variant != "") %>%
  filter(classification == "PATH") %>%
  mutate(variant_id = paste(gene, variant, sep = " | "))

# 10) Patient-level unique specific variant calls
patient_variant <- specific_variants %>%
  distinct(patient_id, variant_id)

# 11) Recurrent specific variants (>=3 patients)
recurrent_variants <- patient_variant %>%
  count(variant_id, name = "n_detected", sort = TRUE) %>%
  filter(n_detected >= 3)

if (nrow(recurrent_variants) == 0) {
  stop("No recurrent specific variants with n>=3 after exclusions/filters.")
}


# 12) PPV + exact 95% CI against Myeloid for each recurrent specific variant
patients <- clin_outcome %>%
  select(patient_id, myeloid)

calc_ppv_variant <- function(v) {
  pos_ids <- patient_variant %>%
    filter(variant_id == v) %>%
    distinct(patient_id)
  
  df <- patients %>%
    left_join(pos_ids %>% mutate(test_pos = 1L), by = "patient_id") %>%
    mutate(test_pos = if_else(is.na(test_pos), 0L, test_pos))
  
  TP <- sum(df$test_pos == 1 & df$myeloid == 1, na.rm = TRUE)
  FP <- sum(df$test_pos == 1 & df$myeloid == 0, na.rm = TRUE)
  n_test_pos <- TP + FP
  
  ppv <- ifelse(n_test_pos > 0, TP / n_test_pos, NA_real_)
  ci <- if (n_test_pos > 0) binom.test(TP, n_test_pos, conf.level = 0.95)$conf.int else c(NA, NA)
  
  data.frame(
    variant_id = v,
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

ppv_results <- do.call(
  rbind,
  lapply(recurrent_variants$variant_id, calc_ppv_variant)
) %>%
  as_tibble() %>%
  arrange(desc(PPV), desc(n_detected))

print(ppv_results)

# 14) Export CSV outputs
write.csv(recurrent_variants, "recurrent_specific_variants_counts.csv", row.names = FALSE)
write.csv(ppv_results, "PPV_recurrent_specific_variants_myeloid_95CI.csv", row.names = FALSE)

cat("\nDone.\nFiles created:\n",
    "- recurrent_specific_variants_counts.csv\n",
    "- PPV_recurrent_specific_variants_myeloid_95CI.csv\n")
