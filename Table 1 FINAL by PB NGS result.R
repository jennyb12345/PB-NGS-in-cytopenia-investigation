# ============================================================
# TABLE 1- BY PB NGS OUTCOME (VCS vs Neg), Exclude VUS and VCS germline only
# ============================================================
# 1) Install/load required packages
required_packages <- c("readxl", "dplyr", "gtsummary", "flextable", "stringr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(gtsummary)
library(flextable)
library(stringr)

# 2) File + tab settings
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
clinical_tab <- "Clinic information + FBC"
pbngs_tab    <- "PB NGS"
id_col       <- "Patient ID"

# 3) Read data
clinical_df <- read_excel(file_path, sheet = clinical_tab)
pbngs_df    <- read_excel(file_path, sheet = pbngs_tab)

# 4) Clean column names
names(clinical_df) <- str_squish(names(clinical_df))
names(pbngs_df)    <- str_squish(names(pbngs_df))

# 5) Keep required clinical fields
clinical_keep <- clinical_df %>%
  transmute(
    join_id      = .data[[id_col]],
    age          = as.numeric(`Age at referral`),
    sex_raw      = as.character(`Sex`),
    bmbx_raw     = as.character(`Proceed to BMBx?`),
    wbc          = as.numeric(`WBC x10^9/L`),
    hb           = as.numeric(`Hb g/L`),
    plts         = as.numeric(`Plts x10^9/L`),
    neuts        = as.numeric(`Neuts x10^9/L`),
    ferritin     = as.numeric(`Ferritin (serum) microgram/L`),
    b12          = as.numeric(`B12 ng/L`),
    folate       = as.numeric(`Folate microgram/L`)
  )

# 6) Keep PB NGS fields (including Final Dx classification from PB NGS tab)
pbngs_keep <- pbngs_df %>%
  transmute(
    join_id             = .data[[id_col]],
    final_dx_pbngs_raw  = str_to_lower(str_squish(as.character(`Final Dx Classification`))),
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`))),
    pbngs_result_raw    = str_to_upper(str_squish(as.character(`PB NGS results`)))
  )

# QC check: inspect classes before filtering
cat("\nVariant classifications (cleaned):\n")
print(table(pbngs_keep$variant_class_clean, useNA = "ifany"))

# QC check: inspect PB NGS Final Dx classification values
cat("\nPB NGS Final Dx Classification (raw cleaned):\n")
print(table(pbngs_keep$final_dx_pbngs_raw, useNA = "ifany"))

# 7) Merge and exclude EXACT unwanted variant categories
df <- clinical_keep %>%
  left_join(pbngs_keep, by = "join_id") %>%
  filter(
    is.na(variant_class_clean) |
      !(variant_class_clean %in% c("vus only", "vcs germline only"))
  )

# QC: show rows kept/excluded
cat("\nRows after exclusions:", nrow(df), "\n")
cat(
  "Excluded rows count (in PB NGS table):",
  sum(pbngs_keep$variant_class_clean %in% c("vus only", "vcs germline only"), na.rm = TRUE),
  "\n"
)

# 8) Create PB NGS grouping + derived characteristics
df <- df %>%
  mutate(
    # Outcome groups for table columns
    ngs_group = case_when(
      pbngs_result_raw == "POS" ~ "Variants of Clin.Signif",
      pbngs_result_raw == "NEG" ~ "No variants of Clin.Signif",
      TRUE ~ NA_character_ ),
    ngs_group = factor(
      ngs_group,
      levels = c("Variants of Clin.Signif", "No variants of Clin.Signif")),
    
    # Standardize sex
    sex = case_when(
      str_to_upper(str_trim(sex_raw)) %in% c("F", "FEMALE") ~ "Female",
      str_to_upper(str_trim(sex_raw)) %in% c("M", "MALE") ~ "Male",
      TRUE ~ NA_character_ ),
    sex = factor(sex, levels = c("Female", "Male")),
    
    # BMBx
    bmbx_yes = if_else(str_to_upper(str_trim(bmbx_raw)) == "Y", "Y", "No"),
    bmbx_yes = factor(bmbx_yes, levels = c("Y", "No")),
    
    # PB NGS Final Dx classification mapped to 3 categories
    diagnostic_category = case_when(
      final_dx_pbngs_raw == "myeloid" ~ "Myeloid",
      final_dx_pbngs_raw == "lymphoid" ~ "Lymphoid",
      final_dx_pbngs_raw %in% c(
        "non-haematological", "non haematological",
        "nonhematological", "non-hematological"
      ) ~ "Non-Haematological",
      TRUE ~ NA_character_ ),
    diagnostic_category = factor(
      diagnostic_category,
      levels = c("Myeloid", "Lymphoid", "Non-Haematological")),
    
    # Three separate diagnostic lines
    dx_myeloid = factor(if_else(diagnostic_category == "Myeloid", "Yes", "No", missing = "No"),
                        levels = c("Yes", "No")),
    dx_lymphoid = factor(if_else(diagnostic_category == "Lymphoid", "Yes", "No", missing = "No"),
                         levels = c("Yes", "No")),
    dx_non_haem = factor(if_else(diagnostic_category == "Non-Haematological", "Yes", "No", missing = "No"),
                         levels = c("Yes", "No")),
    
    # Cytopenia definitions
    anaemia_log = !is.na(hb) & hb < 120,
    thrombocytopenia_log = !is.na(plts) & plts < 150,
    leucopenia_log = !is.na(wbc) & wbc < 3.7,
    neutropenia_log = !is.na(neuts) & neuts < 1.8,
    
    # Third lineage counted as ONE if either leucopenia OR neutropenia is true
    WBC_lineage_low_log = leucopenia_log | neutropenia_log,
    cytopenia_n = anaemia_log + thrombocytopenia_log + WBC_lineage_low_log,
    
    anaemia = factor(if_else(anaemia_log, "Yes", "No"), levels = c("Yes", "No")),
    thrombocytopenia = factor(if_else(thrombocytopenia_log, "Yes", "No"), levels = c("Yes", "No")),
    leucopenia = factor(if_else(leucopenia_log, "Yes", "No"), levels = c("Yes", "No")),
    neutropenia = factor(if_else(neutropenia_log, "Yes", "No"), levels = c("Yes", "No")),
    isolated_cytopenia = factor(if_else(cytopenia_n == 1, "Yes", "No"), levels = c("Yes", "No")),
    bicytopenia = factor(if_else(cytopenia_n == 2, "Yes", "No"), levels = c("Yes", "No")),
    pancytopenia = factor(if_else(cytopenia_n == 3, "Yes", "No"), levels = c("Yes", "No"))) %>%
  filter(!is.na(ngs_group))

# 9) Build table
tbl_ngs <- df %>%
  select(
    ngs_group,
    age, sex, bmbx_yes,
    dx_myeloid, dx_lymphoid, dx_non_haem,
    wbc, hb, plts, neuts, ferritin, b12, folate,
    anaemia, thrombocytopenia, leucopenia, neutropenia,
    isolated_cytopenia, bicytopenia, pancytopenia
  ) %>%
  tbl_summary(
    by = ngs_group,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    value = list(
      bmbx_yes ~ "Y",
      dx_myeloid ~ "Yes",
      dx_lymphoid ~ "Yes",
      dx_non_haem ~ "Yes",
      anaemia ~ "Yes",
      thrombocytopenia ~ "Yes",
      leucopenia ~ "Yes",
      neutropenia ~ "Yes",
      isolated_cytopenia ~ "Yes",
      bicytopenia ~ "Yes",
      pancytopenia ~ "Yes"
    ),
    label = list(
      age ~ "Age at referral",
      sex ~ "Sex",
      bmbx_yes ~ "Underwent BMBx",
      dx_myeloid ~ "Myeloid diagnosis",
      dx_lymphoid ~ "Lymphoid diagnosis",
      dx_non_haem ~ "Non-Haem diagnosis",
      wbc ~ "White blood cells",
      hb ~ "Haemoglobin",
      plts ~ "Platelets",
      neuts ~ "Neutrophils",
      ferritin ~ "Ferritin",
      b12 ~ "B12",
      folate ~ "Folate",
      anaemia ~ "Anaemia",
      thrombocytopenia ~ "Thrombocytopenia",
      leucopenia ~ "Leucopenia",
      neutropenia ~ "Neutropenia",
      isolated_cytopenia ~ "Isolated cytopenia",
      bicytopenia ~ "Bicytopenia",
      pancytopenia ~ "Pancytopenia"
    ),
    missing = "ifany",
    digits = all_continuous() ~ 1
  ) %>%
  add_p(
    test = list(
      age ~ "wilcox.test",
      wbc ~ "wilcox.test",
      hb ~ "wilcox.test",
      plts ~ "wilcox.test",
      neuts ~ "wilcox.test",
      ferritin ~ "wilcox.test",
      b12 ~ "wilcox.test",
      folate ~ "wilcox.test",
      
      sex ~ "chisq.test",
      bmbx_yes ~ "chisq.test",
      anaemia ~ "chisq.test",
      thrombocytopenia ~ "chisq.test",
      leucopenia ~ "chisq.test",
      neutropenia ~ "chisq.test",
      isolated_cytopenia ~ "chisq.test",
      bicytopenia ~ "chisq.test",
      pancytopenia ~ "chisq.test",
      
      dx_myeloid ~ "fisher.test",
      dx_lymphoid ~ "fisher.test",
      dx_non_haem ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall(last = FALSE, col_label = "**Overall**") %>%
  bold_labels()

# 10) View and export
tbl_ngs
save_as_docx(
  as_flex_table(tbl_ngs), path = "Table1_by_PB_NGS_outcome.docx")

cat("\nDone.\nFile created:\n- Table1_by_PB_NGS_outcome.docx\n")