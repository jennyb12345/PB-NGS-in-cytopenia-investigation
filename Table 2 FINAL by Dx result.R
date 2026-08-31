# ============================================================
# # TABLE 2- BY DIAGNOSTIC OUTCOME
# Includes:
# - Myeloid, Lymphoid, No evidence of haem cause (+ Overall)
# - Excludes PB NGS classes: "VUS only", "VCS germline only"
# ============================================================

# 1) Install/load required packages
required_packages <- c("readxl", "dplyr", "stringr", "gtsummary", "flextable")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(stringr)
library(gtsummary)
library(flextable)

# 2) File + tab settings
file_path    <- "MP_PB_audit_cleandata_V2.xlsx"
clinical_tab <- "Clinic information + FBC"
pbngs_tab    <- "PB NGS"
id_col       <- "Patient ID"

# If path issue, use:
# file_path <- file.choose()

# 3) Read sheets
clinical_df <- read_excel(file_path, sheet = clinical_tab)
pbngs_df    <- read_excel(file_path, sheet = pbngs_tab)

# 4) Clean column names
names(clinical_df) <- str_squish(names(clinical_df))
names(pbngs_df)    <- str_squish(names(pbngs_df))

# 5) PB NGS fields + exclusions
pbngs_keep <- pbngs_df %>%
  transmute(
    join_id             = .data[[id_col]],
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`))),
    pbngs_result_raw    = str_to_upper(str_squish(as.character(`PB NGS results`)))
  ) %>%
  filter(
    is.na(variant_class_clean) |
      !(variant_class_clean %in% c("vus only", "vcs germline only"))
  ) %>%
  mutate(
    var_clinsig_pos = if_else(pbngs_result_raw == "POS", "Yes", "No"),
    var_clinsig_neg = if_else(pbngs_result_raw == "NEG", "Yes", "No"),
    var_clinsig_pos = factor(var_clinsig_pos, levels = c("Yes", "No")),
    var_clinsig_neg = factor(var_clinsig_neg, levels = c("Yes", "No"))
  ) %>%
  distinct(join_id, .keep_all = TRUE)

# 6) Clinical fields
clinical_keep <- clinical_df %>%
  transmute(
    join_id = .data[[id_col]],
    final_dx_raw = as.character(`Final Dx Classification`),
    age = as.numeric(`Age at referral`),
    sex_raw = as.character(`Sex`),
    bmbx_raw = as.character(`Proceed to BMBx?`),
    wbc = as.numeric(`WBC x10^9/L`),
    hb = as.numeric(`Hb g/L`),
    plts = as.numeric(`Plts x10^9/L`),
    neuts = as.numeric(`Neuts x10^9/L`),
    ferritin = as.numeric(`Ferritin (serum) microgram/L`),
    b12 = as.numeric(`B12 ng/L`),
    folate = as.numeric(`Folate microgram/L`)
  )

# 7) Merge + derive variables
df <- clinical_keep %>%
  inner_join(pbngs_keep, by = "join_id") %>%
  mutate(
    dx_group = case_when(
      str_squish(final_dx_raw) == "Myeloid" ~ "Myeloid",
      str_squish(final_dx_raw) == "Lymphoid" ~ "Lymphoid",
      str_squish(final_dx_raw) == "Non-Haematological" ~ "Non-haem diagnosis",
      TRUE ~ NA_character_
    ),
    dx_group = factor(dx_group, levels = c("Myeloid", "Lymphoid", "Non-haem diagnosis")),
    
    sex = case_when(
      str_to_upper(str_trim(sex_raw)) %in% c("F", "FEMALE") ~ "Female",
      str_to_upper(str_trim(sex_raw)) %in% c("M", "MALE") ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    
    bmbx_yes = if_else(str_to_upper(str_trim(bmbx_raw)) == "Y", "Y", "No"),
    bmbx_yes = factor(bmbx_yes, levels = c("Y", "No")),
    
    # Cytopenia definitions
    anaemia_log = !is.na(hb) & hb < 120,
    thrombocytopenia_log = !is.na(plts) & plts < 150,
    leucopenia_log = !is.na(wbc) & wbc < 3.7,
    neutropenia_log = !is.na(neuts) & neuts < 1.8,
    
    # Third lineage = leucopenia OR neutropenia
    wbc_lineage_low_log = leucopenia_log | neutropenia_log,
    cytopenia_n = anaemia_log + thrombocytopenia_log + wbc_lineage_low_log,
    
    anaemia = factor(if_else(anaemia_log, "Yes", "No"), levels = c("Yes", "No")),
    thrombocytopenia = factor(if_else(thrombocytopenia_log, "Yes", "No"), levels = c("Yes", "No")),
    leucopenia = factor(if_else(leucopenia_log, "Yes", "No"), levels = c("Yes", "No")),
    neutropenia = factor(if_else(neutropenia_log, "Yes", "No"), levels = c("Yes", "No")),
    isolated_cytopenia = factor(if_else(cytopenia_n == 1, "Yes", "No"), levels = c("Yes", "No")),
    bicytopenia = factor(if_else(cytopenia_n == 2, "Yes", "No"), levels = c("Yes", "No")),
    pancytopenia = factor(if_else(cytopenia_n == 3, "Yes", "No"), levels = c("Yes", "No"))
  ) %>%
  filter(!is.na(dx_group))

# Optional QC
cat("\nGroup counts:\n")
print(table(df$dx_group, useNA = "ifany"))

# 8) Build summary table
tbl_dx <- df %>%
  select(
    dx_group,
    age, sex, bmbx_yes,
    wbc, hb, plts, neuts, ferritin, b12, folate,
    anaemia, thrombocytopenia, leucopenia, neutropenia,
    isolated_cytopenia, bicytopenia, pancytopenia,
    var_clinsig_pos, var_clinsig_neg
  ) %>%
  tbl_summary(
    by = dx_group,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    value = list(
      bmbx_yes ~ "Y",
      anaemia ~ "Yes",
      thrombocytopenia ~ "Yes",
      leucopenia ~ "Yes",
      neutropenia ~ "Yes",
      isolated_cytopenia ~ "Yes",
      bicytopenia ~ "Yes",
      pancytopenia ~ "Yes",
      var_clinsig_pos ~ "Yes",
      var_clinsig_neg ~ "Yes"
    ),
    label = list(
      age ~ "Age at referral",
      sex ~ "Sex",
      bmbx_yes ~ "Underwent BMBx",
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
      pancytopenia ~ "Pancytopenia",
      var_clinsig_pos ~ "Variant of Clin.Signif ",
      var_clinsig_neg ~ "No variants of Clin.Signif"
    ),
    missing = "ifany",
    digits = all_continuous() ~ 1
  ) %>%
  add_p(
    test = list(
      age ~ "kruskal.test",
      wbc ~ "kruskal.test",
      hb ~ "kruskal.test",
      plts ~ "kruskal.test",
      neuts ~ "kruskal.test",
      ferritin ~ "kruskal.test",
      b12 ~ "kruskal.test",
      folate ~ "kruskal.test",
      sex ~ "fisher.test",
      bmbx_yes ~ "fisher.test",
      anaemia ~ "fisher.test",
      thrombocytopenia ~ "fisher.test",
      leucopenia ~ "fisher.test",
      neutropenia ~ "fisher.test",
      isolated_cytopenia ~ "fisher.test",
      bicytopenia ~ "fisher.test",
      pancytopenia ~ "fisher.test",
      var_clinsig_pos ~ "fisher.test",
      var_clinsig_neg ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall(last = FALSE, col_label = "**Overall**") %>%
  bold_labels()

# 9) View & export table
tbl_dx
save_as_docx(as_flex_table(tbl_dx), path = "Table2_by_Diagnostic_Outcome.docx")

cat("\nDone.\nFile created:\n- Table2_by_Diagnostic_Outcome.docx\n")