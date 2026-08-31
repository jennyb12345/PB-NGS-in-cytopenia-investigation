# ============================================================
# GENE-SET MODELS A-E ON ANALYSIS COHORT
# Condition+ = Myeloid
# Test+ = >=1 PATH variant in model genes
# Exclusions: Variant Classifications = "VUS only" or "VCS germline only"
# Cohort: Whole eligible cohort
# Outputs: PPV, NPV, LR+, LR- with 95% CI
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# 1) File + tabs
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
pb_tab       <- "PB NGS"
clinical_tab <- "Clinic information + FBC"
id_col       <- "Patient ID"

# 2) Read data
pb <- read_excel(file_path, sheet = pb_tab)
cl <- read_excel(file_path, sheet = clinical_tab)

# 3) Clean names
names(pb) <- str_squish(names(pb))
names(cl) <- str_squish(names(cl))

# 4) Detect slot columns
gene_cols <- names(pb)[str_detect(names(pb), regex("^Gene\\s*\\d+$", ignore_case = TRUE))]
cls_cols  <- names(pb)[str_detect(names(pb), regex("^Classification\\s*\\d+$", ignore_case = TRUE))]
if (length(gene_cols) == 0 || length(cls_cols) == 0) {
  stop("Could not detect Gene/Classification slot columns.")
}

# 5) PB exclusions + valid PB results (same cohort logic as original analysis)
pb_filt <- pb %>%
  mutate(
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`))),
    pb_result = str_to_upper(str_squish(as.character(`PB NGS results`))),
    patient_id = .data[[id_col]]
  ) %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean)) %>%
  filter(pb_result %in% c("POS", "NEG")) %>%
  distinct(patient_id, .keep_all = TRUE)

# 6) Long format align Gene N with Classification N; keep PATH calls
pb_long_gene <- pb_filt %>%
  select(all_of(c("patient_id", gene_cols))) %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "gene_col", values_to = "gene") %>%
  mutate(slot = str_extract(gene_col, "\\d+"))

pb_long_cls <- pb_filt %>%
  select(all_of(c("patient_id", cls_cols))) %>%
  pivot_longer(cols = all_of(cls_cols), names_to = "cls_col", values_to = "classification") %>%
  mutate(
    slot = str_extract(cls_col, "\\d+"),
    classification = str_to_upper(str_squish(as.character(classification)))
  ) %>%
  select(patient_id, slot, classification)

pb_long <- pb_long_gene %>%
  left_join(pb_long_cls, by = c("patient_id", "slot")) %>%
  transmute(
    patient_id,
    gene = str_to_upper(str_squish(as.character(gene))),
    classification
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  filter(classification == "PATH") %>%
  distinct(patient_id, gene)

# 7) Clinical outcomes; restrict to same PB cohort IDs
clin_keep <- cl %>%
  transmute(
    patient_id = .data[[id_col]],
    final_dx = str_squish(as.character(`Final Dx Classification`))
  ) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  inner_join(pb_filt %>% select(patient_id), by = "patient_id") %>%
  mutate(
    condition_pos = if_else(str_to_lower(final_dx) == "myeloid", 1L, 0L)
  )

# 8) Model definitions
model_genes <- list(
  A = c("TET2", "SRSF2", "ASXL1", "DNMT3A", "SF3B1"),
  B = c("TET2", "SRSF2", "ASXL1", "DNMT3A", "SF3B1", "U2AF1", "CUX1", "TP53", "CBL", "KRAS"),
  C = c("TET2", "SRSF2", "ASXL1", "DNMT3A", "SF3B1", "U2AF1", "CUX1", "TP53", "CBL", "KRAS",
        "STAG2", "IDH2", "PTPN11", "ZRSR2", "DDX41"),
  D = c("TET2", "SRSF2", "ASXL1", "DNMT3A", "SF3B1", "U2AF1", "CUX1", "TP53", "CBL", "KRAS",
        "STAG2", "IDH2", "PTPN11", "ZRSR2", "DDX41", "EZH2", "GNAS", "JAK2", "PPM1D", "NPM1", "RUNX1"),
  E = NULL # all genes
)

# 9) Metric function
calc_model_metrics <- function(model_name, genes_vec) {
  
  # Test-positive IDs for this gene model
  if (is.null(genes_vec)) {
    test_pos_ids <- pb_long %>% distinct(patient_id)
    genes_used <- "ALL GENES"
  } else {
    test_pos_ids <- pb_long %>%
      filter(gene %in% genes_vec) %>%
      distinct(patient_id)
    genes_used <- paste(genes_vec, collapse = ", ")
  }
  
  # Full eligible cohort base
  df <- clin_keep %>%
    left_join(test_pos_ids %>% mutate(test_pos = 1L), by = "patient_id") %>%
    mutate(
      test_pos = if_else(is.na(test_pos), 0L, test_pos),
      test_neg = 1L - test_pos
    )
  
  # 2x2
  TP <- sum(df$test_pos == 1 & df$condition_pos == 1, na.rm = TRUE)
  FP <- sum(df$test_pos == 1 & df$condition_pos == 0, na.rm = TRUE)
  FN <- sum(df$test_pos == 0 & df$condition_pos == 1, na.rm = TRUE)
  TN <- sum(df$test_pos == 0 & df$condition_pos == 0, na.rm = TRUE)
  
  # PPV / NPV
  ppv <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  npv <- ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)
  
  ppv_ci <- if ((TP + FP) > 0) binom.test(TP, TP + FP)$conf.int else c(NA, NA)
  npv_ci <- if ((TN + FN) > 0) binom.test(TN, TN + FN)$conf.int else c(NA, NA)
  
  # LR
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  
  lr_pos <- ifelse(!is.na(sens) & !is.na(spec) & (1 - spec) > 0, sens / (1 - spec), NA_real_)
  lr_neg <- ifelse(!is.na(sens) & !is.na(spec) & spec > 0, (1 - sens) / spec, NA_real_)
  
  # LR CI (log method + 0.5 correction if any zero cell)
  a <- TP; b <- FP; c <- FN; d <- TN
  if (any(c(a,b,c,d) == 0)) { a <- a+0.5; b <- b+0.5; c <- c+0.5; d <- d+0.5 }
  
  se_log_lrpos <- sqrt((1/a) - (1/(a+c)) + (1/b) - (1/(b+d)))
  se_log_lrneg <- sqrt((1/c) - (1/(a+c)) + (1/d) - (1/(b+d)))
  
  sens_c <- a / (a + c)
  spec_c <- d / (b + d)
  lr_pos_c <- sens_c / (1 - spec_c)
  lr_neg_c <- (1 - sens_c) / spec_c
  
  lr_pos_ci <- c(exp(log(lr_pos_c) - 1.96 * se_log_lrpos),
                 exp(log(lr_pos_c) + 1.96 * se_log_lrpos))
  lr_neg_ci <- c(exp(log(lr_neg_c) - 1.96 * se_log_lrneg),
                 exp(log(lr_neg_c) + 1.96 * se_log_lrneg))
  
  data.frame(
    Model = model_name,
    Genes_in_model = genes_used,
    N_total_patients = nrow(df),
    N_test_positive = TP + FP,
    N_test_negative = TN + FN,
    TP_myeloid = TP,
    FP_non_myeloid = FP,
    FN_myeloid = FN,
    TN_non_myeloid = TN,
    
    PPV_raw = ppv,
    PPV_95CI_low_raw = ppv_ci[1],
    PPV_95CI_high_raw = ppv_ci[2],
    
    NPV_raw = npv,
    NPV_95CI_low_raw = npv_ci[1],
    NPV_95CI_high_raw = npv_ci[2],
    
    LR_pos_raw = lr_pos,
    LR_pos_95CI_low_raw = lr_pos_ci[1],
    LR_pos_95CI_high_raw = lr_pos_ci[2],
    
    LR_neg_raw = lr_neg,
    LR_neg_95CI_low_raw = lr_neg_ci[1],
    LR_neg_95CI_high_raw = lr_neg_ci[2],
    
    PPV_display = ifelse(is.na(ppv), NA,
                         paste0(round(100 * ppv, 1), "% (", round(100 * ppv_ci[1], 1), "–", round(100 * ppv_ci[2], 1), ")")),
    NPV_display = ifelse(is.na(npv), NA,
                         paste0(round(100 * npv, 1), "% (", round(100 * npv_ci[1], 1), "–", round(100 * npv_ci[2], 1), ")")),
    LR_pos_display = ifelse(is.na(lr_pos), NA,
                            paste0(round(lr_pos, 2), " (", round(lr_pos_ci[1], 2), "–", round(lr_pos_ci[2], 2), ")")),
    LR_neg_display = ifelse(is.na(lr_neg), NA,
                            paste0(round(lr_neg, 2), " (", round(lr_neg_ci[1], 2), "–", round(lr_neg_ci[2], 2), ")")),
    
    stringsAsFactors = FALSE
  )
}

# 10) Run models A-E
results <- imap_dfr(model_genes, ~calc_model_metrics(.y, .x))

# 11) Print concise summary
print(results %>%
        select(Model, N_total_patients, N_test_positive, N_test_negative,
               PPV_display, NPV_display, LR_pos_display, LR_neg_display))

# 12) Save output
write.csv(results, "Gene_set_Models_Myeloid_PPV_NPV_LR.csv", row.names = FALSE)

cat("\nDone.\nFile created:\n- Gene_set_Models_Myeloid_PPV_NPV_LR.csv\n")
