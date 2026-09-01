# ============================================================
# MODELS (all outcome = Myeloid):
# 1) >=2 "path" variants (existing model)
# 2) Any "path" variant with VAF > 10
# 3) Any "path" variant with VAF > 20
#
# Exclusion first: Variant Classifications = "VUS only" or "VCS germline only"
# Outputs: PPV, NPV, LR+, LR- with 95% CIs for each model
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

# 4) Clean names
names(clinical_df) <- str_squish(names(clinical_df))
names(pb_df)       <- str_squish(names(pb_df))

# 5) Exclude VUS-only / VCS germline-only
pb_df <- pb_df %>%
  mutate(
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

# 6) Find classification and VAF columns
class_cols <- grep("^Classification\\s*\\d+$", names(pb_df), value = TRUE, ignore.case = TRUE)
vaf_cols   <- grep("^VAF\\s*\\d+$", names(pb_df), value = TRUE, ignore.case = TRUE)

if (length(class_cols) == 0) stop("No 'Classification 1/2/..' columns found.")
if (length(vaf_cols) == 0) stop("No 'VAF 1/2/..' columns found.")

# Match by numeric suffix so Classification k pairs with VAF k
get_suffix_num <- function(x) as.integer(str_extract(x, "\\d+"))
class_map <- data.frame(class_col = class_cols, k = get_suffix_num(class_cols))
vaf_map   <- data.frame(vaf_col   = vaf_cols,   k = get_suffix_num(vaf_cols))

pairs <- inner_join(class_map, vaf_map, by = "k") %>% arrange(k)

if (nrow(pairs) == 0) stop("No matching Classification/VAF pairs found by numeric suffix.")

# 7) Build PB model flags
pb_keep <- pb_df %>%
  mutate(patient_id = .data[[id_col]]) %>%
  rowwise() %>%
  mutate(
    # Existing model: >=2 path
    n_path = sum(
      str_to_lower(str_squish(as.character(c_across(all_of(class_cols))))) == "path",
      na.rm = TRUE
    ),
    test_ge2_path = n_path >= 2,
    
    # New model: any path with VAF > 10
    test_any_path_vaf_gt10 = {
      .row <- pick(everything())
      any_hit <- FALSE
      for (i in seq_len(nrow(pairs))) {
        cval <- str_to_lower(str_squish(as.character(.row[[pairs$class_col[i]]])))
        vraw <- as.character(.row[[pairs$vaf_col[i]]])
        vnum <- suppressWarnings(as.numeric(str_replace_all(vraw, "%", "")))
        if (!is.na(cval) && cval == "path" && !is.na(vnum) && vnum > 10) {
          any_hit <- TRUE
          break
        }
      }
      any_hit
    },
    
    # New model: any path with VAF > 20
    test_any_path_vaf_gt20 = {
      .row <- pick(everything())
      any_hit <- FALSE
      for (i in seq_len(nrow(pairs))) {
        cval <- str_to_lower(str_squish(as.character(.row[[pairs$class_col[i]]])))
        vraw <- as.character(.row[[pairs$vaf_col[i]]])
        vnum <- suppressWarnings(as.numeric(str_replace_all(vraw, "%", "")))
        if (!is.na(cval) && cval == "path" && !is.na(vnum) && vnum > 20) {
          any_hit <- TRUE
          break
        }
      }
      any_hit
    }
  ) %>%
  ungroup() %>%
  select(patient_id, n_path, test_ge2_path, test_any_path_vaf_gt10, test_any_path_vaf_gt20) %>%
  distinct(patient_id, .keep_all = TRUE)

# 8) Clinical outcome
clin_keep <- clinical_df %>%
  transmute(
    patient_id = .data[[id_col]],
    final_dx = str_squish(as.character(`Final Dx Classification`)),
    disease_pos = final_dx == "Myeloid"
  ) %>%
  mutate(disease_neg = !disease_pos) %>%
  distinct(patient_id, .keep_all = TRUE)

# 9) Merge
df <- pb_keep %>%
  inner_join(clin_keep, by = "patient_id")

# ---------- helper to compute one model ----------
compute_one_model <- function(data, test_col, model_name, small_n_threshold = 20) {
  test_pos <- data[[test_col]]
  test_neg <- !test_pos
  
  TP <- sum(test_pos & data$disease_pos, na.rm = TRUE)
  FP <- sum(test_pos & data$disease_neg, na.rm = TRUE)
  FN <- sum(test_neg & data$disease_pos, na.rm = TRUE)
  TN <- sum(test_neg & data$disease_neg, na.rm = TRUE)
  N  <- nrow(data)
  
  PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)
  
  ppv_ci <- if ((TP + FP) > 0) binom.test(TP, TP + FP, conf.level = 0.95)$conf.int else c(NA_real_, NA_real_)
  npv_ci <- if ((TN + FN) > 0) binom.test(TN, TN + FN, conf.level = 0.95)$conf.int else c(NA_real_, NA_real_)
  
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  
  LR_pos <- ifelse(!is.na(sens) & !is.na(spec) & (1 - spec) > 0, sens / (1 - spec), NA_real_)
  LR_neg <- ifelse(!is.na(sens) & !is.na(spec) & spec > 0, (1 - sens) / spec, NA_real_)
  
  # LR CI (log method; 0.5 correction if any zero cells)
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
  
  lrpos_ci <- c(exp(log(LR_pos_c) - 1.96 * se_log_lrpos),
                exp(log(LR_pos_c) + 1.96 * se_log_lrpos))
  lrneg_ci <- c(exp(log(LR_neg_c) - 1.96 * se_log_lrneg),
                exp(log(LR_neg_c) + 1.96 * se_log_lrneg))
  
  warning_msgs <- c()
  if (N < small_n_threshold) warning_msgs <- c(warning_msgs, paste0("Small sample size (N=", N, ")"))
  if (any(c(TP, FP, FN, TN) == 0)) warning_msgs <- c(warning_msgs, "Zero cell present; 0.5 correction applied for LR CI")
  warning_text <- ifelse(length(warning_msgs) == 0, "None", paste(warning_msgs, collapse = " | "))
  
  counts <- data.frame(
    Model = model_name, N = N, TP = TP, FP = FP, FN = FN, TN = TN, Warning = warning_text
  )
  
  results <- data.frame(
    Model = model_name,
    Metric = c("PPV", "NPV", "LR+", "LR-"),
    Estimate_raw = c(PPV, NPV, LR_pos, LR_neg),
    CI_low_raw = c(ppv_ci[1], npv_ci[1], lrpos_ci[1], lrneg_ci[1]),
    CI_high_raw = c(ppv_ci[2], npv_ci[2], lrpos_ci[2], lrneg_ci[2]),
    Warning = warning_text,
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

# 10) Run all 3 models
out1 <- compute_one_model(df, "test_ge2_path", ">=2 path variants")
out2 <- compute_one_model(df, "test_any_path_vaf_gt10", "Any path variant with VAF > 10")
out3 <- compute_one_model(df, "test_any_path_vaf_gt20", "Any path variant with VAF > 20")

counts_all  <- bind_rows(out1$counts, out2$counts, out3$counts)
results_all <- bind_rows(out1$results, out2$results, out3$results)

cat("\n=== COUNTS USED ===\n")
print(counts_all)

cat("\n=== FINAL METRICS ===\n")
print(results_all[, c("Model", "Metric", "Summary_display", "Warning")])

cat("\n=== RAW METRICS ===\n")
print(results_all)

write.csv(results_all, "PPV_NPV_2+ and VAF_Myeloid_metrics_with_CI.csv", row.names = FALSE)
write.csv(counts_all,  "PPV_NPV_2+ and VAF_Myeloid_counts.csv", row.names = FALSE)

cat("\nDone.\nFiles created:\n",
    "- PPV_NPV_2+ and VAF_Myeloid_metrics_with_CI.csv\n",
    "- PPV_NPV_2+ and VAF_Myeloid_3models_counts.csv\n")
