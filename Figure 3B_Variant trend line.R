# ============================================================
# geom_smooth trend plot:
# all specific Path variants vs number of patients
# ============================================================

# 1) Install/load packages
required_packages <- c("readxl", "dplyr", "stringr", "tidyr", "ggplot2")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 2) File + tab settings
file_path <- "MP_PB_audit_cleandata_V2.xlsx"
pb_tab <- "PB NGS"
id_col <- "Patient ID"

# If path issue:
# file_path <- file.choose()

# 3) Read PB NGS tab
pb <- read_excel(file_path, sheet = pb_tab)

# 4) Clean names
names(pb) <- str_squish(names(pb))

# 5) Detect slot columns
gene_cols <- names(pb)[str_detect(names(pb), regex("^Gene\\s*\\d+$", ignore_case = TRUE))]
var_cols  <- names(pb)[str_detect(names(pb), regex("^Variant\\s*\\d+$", ignore_case = TRUE))]
cls_cols  <- names(pb)[str_detect(names(pb), regex("^Classification\\s*\\d+$", ignore_case = TRUE))]

if (length(gene_cols) == 0 || length(var_cols) == 0 || length(cls_cols) == 0) {
  stop("Could not detect Gene/Variant/Classification slot columns.")
}

# 6) Sample-level exclusions
pb_filt <- pb %>%
  mutate(
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

# 7) Long format with slot alignment
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

# 8) Keep only Path variants and create specific variant ID
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

# 9) Count unique patients per specific variant
variant_counts <- specific_variants %>%
  distinct(patient_id, variant_id) %>%
  count(variant_id, name = "n_patients", sort = TRUE) %>%
  mutate(variant_rank = row_number())

# 10) Plot geom_smooth trend
p <- ggplot(variant_counts, aes(x = variant_rank, y = n_patients)) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 2) +
  labs(
    title = "Variant Detection trend",
    x = "Individual variants",
    y = "Number of patients"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )

print(p)

# 11) Save outputs
ggsave("Figure 3B_variant_trend_geom_smooth.png", p, width = 5, height = 3.5, dpi = 600)
write.csv(variant_counts, "specific_variant_patient_counts.csv", row.names = FALSE)

cat("\nDone.\nFiles created:\n",
    "- Figure 3B_variant_trend_geom_smooth.png\n",
    "- specific_variant_patient_counts.csv\n")