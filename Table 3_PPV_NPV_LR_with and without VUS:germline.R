# ============================================================
# Table 3 PB NGS POS predicts Myeloid
# Outputs: PPV, NPV, LR+, LR- with 95% CIs
# Cohort A: EXCLUDING "VUS only" + "VCS germline only"
# Cohort B: INCLUDING "VUS only" + "VCS germline only"
# ============================================================

# 1) Packages
required_packages <- c("readxl", "dplyr", "stringr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(stringr)

# 2) File + tabs
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
clinical_tab <- "Clinic information + FBC"
pb_tab       <- "PB NGS"
id_col       <- "Patient ID"

# 3) Read data
clinical_df <- read_excel(file_path, sheet = clinical_tab)
pb_df       <- read_excel(file_path, sheet = pb_tab)

# 4) Clean column names
names(clinical_df) <- str_squish(names(clinical_df))
names(pb_df)       <- str_squish(names(pb_df))

# 5) Build PB base table (no exclusion yet)
pb_base <- pb_df %>%
  transmute(
    patient_id = .data[[id_col]],
    pb_result = str_to_upper(str_squish(as.character(`PB NGS results`))),
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(pb_result %in% c("POS", "NEG")) %>%
  distinct(patient_id, .keep_all = TRUE)

# 6) Clinical table
clin_keep <- clinical_df %>%
  transmute(
    patient_id = .data[[id_col]],
    final_dx = str_squish(as.character(`Final Dx Classification`))
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

# ---------- helper function ----------
compute_metrics <- function(pb_data, clin_data, cohort_label) {
  df <- pb_data %>%
    inner_join(clin_data, by = "patient_id") %>%
    mutate(
      disease_pos = final_dx == "Myeloid",
      disease_neg = !disease_pos,
      test_pos = pb_result == "POS",
      test_neg = pb_result == "NEG"
    )
  
  TP <- sum(df$test_pos & df$disease_pos, na.rm = TRUE)
  FP <- sum(df$test_pos & df$disease_neg, na.rm = TRUE)
  FN <- sum(df$test_neg & df$disease_pos, na.rm = TRUE)
  TN <- sum(df$test_neg & df$disease_neg, na.rm = TRUE)
  
  # Core metrics
  PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)
  
  ppv_ci <- if ((TP + FP) > 0) binom.test(TP, TP + FP, conf.level = 0.95)$conf.int else c(NA_real_, NA_real_)
  npv_ci <- if ((TN + FN) > 0) binom.test(TN, TN + FN, conf.level = 0.95)$conf.int else c(NA_real_, NA_real_)
  
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  
  LR_pos <- ifelse(!is.na(sens) & !is.na(spec) & (1 - spec) > 0, sens / (1 - spec), NA_real_)
  LR_neg <- ifelse(!is.na(sens) & !is.na(spec) & spec > 0, (1 - sens) / spec, NA_real_)
  
  # LR CIs (log method with 0.5 correction if any zero cell)
  a <- TP; b <- FP; c <- FN; d <- TN
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  
  sens_c <- a / (a + c)
  spec_c <- d / (b + d)
  LR_pos_c <- sens_c / (1 - spec_c)
  LR_neg_c <- (1 - sens_c) / spec_c
  
  se_log_lrpos <- sqrt((1/a) - (1/(a+c)) + (1/b) - (1/(b+d)))
  se_log_lrneg <- sqrt((1/c) - (1/(a+c)) + (1/d) - (1/(b+d)))
  
  lrpos_ci <- c(
    exp(log(LR_pos_c) - 1.96 * se_log_lrpos),
    exp(log(LR_pos_c) + 1.96 * se_log_lrpos)
  )
  lrneg_ci <- c(
    exp(log(LR_neg_c) - 1.96 * se_log_lrneg),
    exp(log(LR_neg_c) + 1.96 * se_log_lrneg)
  )
  
  counts <- data.frame(
    Cohort = cohort_label,
    N = nrow(df),
    TP = TP, FP = FP, FN = FN, TN = TN
  )
  
  results <- data.frame(
    Cohort = cohort_label,
    Metric = c("PPV", "NPV", "LR+", "LR-"),
    Estimate_raw = c(PPV, NPV, LR_pos, LR_neg),
    CI_low_raw = c(ppv_ci[1], npv_ci[1], lrpos_ci[1], lrneg_ci[1]),
    CI_high_raw = c(ppv_ci[2], npv_ci[2], lrpos_ci[2], lrneg_ci[2]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Summary_display = case_when(
        Metric %in% c("PPV", "NPV") ~
          paste0(round(100 * Estimate_raw, 2), "% (", round(100 * CI_low_raw, 2), "–", round(100 * CI_high_raw, 2), ")"),
        TRUE ~
          paste0(round(Estimate_raw, 3), " (", round(CI_low_raw, 3), "–", round(CI_high_raw, 3), ")")
      )
    )
  
  list(counts = counts, results = results)
}

# 7) Cohort A: Excluding VUS/VCS germline only
pb_excluding <- pb_base %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

out_excluding <- compute_metrics(pb_excluding, clin_keep, "Excluding VUS only + VCS germline only")

# 8) Cohort B: Including VUS/VCS germline only
pb_including <- pb_base
out_including <- compute_metrics(pb_including, clin_keep, "Including VUS only + VCS germline only")

# 9) Combine outputs
counts_all <- bind_rows(out_excluding$counts, out_including$counts)
results_all <- bind_rows(out_excluding$results, out_including$results)

cat("\n=== COUNTS USED ===\n")
print(counts_all)

cat("\n=== FINAL METRICS ===\n")
print(results_all[, c("Cohort", "Metric", "Summary_display")])

cat("\n=== RAW METRICS ===\n")
print(results_all)

# 10) Save files
write.csv(results_all, "PB_NGS_Myeloid_metrics_compare_excluding_vs_including_VUS_VCS.csv", row.names = FALSE)
write.csv(counts_all, "PB_NGS_Myeloid_counts_compare_excluding_vs_including_VUS_VCS.csv", row.names = FALSE)

cat("\nDone.\nFiles created:\n",
    "- PB_NGS_Myeloid_metrics_compare_excluding_vs_including_VUS_VCS.csv\n",
    "- PB_NGS_Myeloid_counts_compare_excluding_vs_including_VUS_VCS.csv\n")