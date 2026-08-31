# ============================================================
# FINAL ONCOPLOT
# 'path'/VCS variants only. Categorise by diagnosis, gender and frequency of genes
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 1) File + sheet settings
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
var_cols  <- names(pb)[str_detect(names(pb), regex("^Variant\\s*\\d+$", ignore_case = TRUE))]
cls_cols  <- names(pb)[str_detect(names(pb), regex("^Classification\\s*\\d+$", ignore_case = TRUE))]

if (length(gene_cols) == 0 || length(var_cols) == 0 || length(cls_cols) == 0) {
  stop("Could not find Gene/Variant/Classification slot columns.")
}

# 5) Clinical annotation (unchanged)
clin_annot <- cl %>%
  transmute(
    patient_id = .data[[id_col]],
    diagnosis_raw = as.character(`Final Dx Classification`),
    sex_raw = as.character(`Sex`)
  ) %>%
  mutate(
    diagnosis = case_when(
      str_squish(diagnosis_raw) == "Myeloid" ~ "Myeloid",
      str_squish(diagnosis_raw) == "Lymphoid" ~ "Lymphoid",
      str_squish(diagnosis_raw) %in% c("Non-Haematological", "No evidence of haem cause") ~ "No haem cause",
      TRUE ~ "Other"
    ),
    sex = case_when(
      str_to_upper(str_trim(sex_raw)) %in% c("F", "FEMALE") ~ "Female",
      str_to_upper(str_trim(sex_raw)) %in% c("M", "MALE") ~ "Male",
      TRUE ~ "Unknown"
    )
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

# 6) PB filter (sample-level exclusion first)
pb_filt <- pb %>%
  mutate(
    pb_result = str_to_upper(str_squish(as.character(`PB NGS results`))),
    variant_class_clean = str_to_lower(str_squish(as.character(`Variant Classifications`)))
  ) %>%
  filter(pb_result == "POS") %>%
  filter(!(variant_class_clean %in% c("vus only", "vcs germline only")) | is.na(variant_class_clean))

# 7) Long format align Gene N / Variant N / Classification N
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
  mutate(slot = str_extract(cls_col, "\\d+"),
         classification = str_to_upper(str_squish(as.character(classification)))) %>%
  select(all_of(id_col), slot, classification)

# 8) Keep only Path variants at slot-level
mut_long <- pb_long_gene %>%
  left_join(pb_long_var, by = c(id_col, "slot")) %>%
  left_join(pb_long_cls, by = c(id_col, "slot")) %>%
  transmute(
    patient_id = .data[[id_col]],
    gene = str_squish(as.character(gene)),
    variant = str_squish(as.character(variant)),
    classification = classification
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(variant), variant != "") %>%
  filter(classification == "PATH")

# 8) Count variants per patient-gene + annotation
oncodata <- mut_long %>%
  group_by(patient_id, gene) %>%
  summarise(n_variants = n(), .groups = "drop") %>%
  left_join(clin_annot, by = "patient_id") %>%
  filter(diagnosis %in% c("Myeloid", "Lymphoid", "No haem cause")) %>%
  mutate(
    variant_multiplicity = if_else(n_variants > 1, "Multiple variants", "Single variant")
  )

# 9) Gene order on Y axis (most frequent top)
gene_order <- oncodata %>%
  distinct(patient_id, gene) %>%
  count(gene, sort = TRUE) %>%
  pull(gene)

oncodata <- oncodata %>%
  mutate(gene = factor(gene, levels = rev(gene_order)))

# 10) Patient order by diagnosis, sex, then gene-frequency

# Create global gene frequency ranking (most frequent = rank 1)
gene_rank_tbl <- oncodata %>%
  distinct(patient_id, gene) %>%
  count(gene, sort = TRUE) %>%
  mutate(
    gene_rank = row_number(),
    gene_upper = str_to_upper(gene)
  )

# Patient-level priority score based on ranked genes
# Higher score = has more high-frequency genes
patient_priority <- oncodata %>%
  mutate(gene_upper = str_to_upper(gene)) %>%
  left_join(gene_rank_tbl %>% select(gene_upper, gene_rank), by = "gene_upper") %>%
  group_by(patient_id, diagnosis, sex) %>%
  summarise(
    n_genes_mutated = n_distinct(gene_upper),
    total_variant_events = sum(n_variants),
    # inverse rank weighting: rank 1 contributes most
    freq_priority_score = sum(1 / gene_rank, na.rm = TRUE),
    .groups = "drop"
  )

# Final patient order:
# diagnosis block -> sex block -> highest frequency-priority first
patient_order <- patient_priority %>%
  mutate(
    dx_rank = case_when(
      diagnosis == "Myeloid" ~ 1L,
      diagnosis == "Lymphoid" ~ 2L,
      diagnosis == "No haem cause" ~ 3L,
      TRUE ~ 4L
    ),
    sex_rank = case_when(
      sex == "Male" ~ 1L,
      sex == "Female" ~ 2L,
      TRUE ~ 3L
    ),
    patient_id_num = suppressWarnings(as.numeric(as.character(patient_id)))
  ) %>%
  arrange(
    dx_rank,
    sex_rank,
    desc(freq_priority_score),
    desc(n_genes_mutated),
    desc(total_variant_events),
    patient_id_num,
    patient_id
  ) %>%
  pull(patient_id)

oncodata <- oncodata %>%
  mutate(patient_id = factor(patient_id, levels = patient_order))

patient_annot <- oncodata %>%
  distinct(patient_id, diagnosis, sex)

# Export patient ordering scores for audit
write.csv(patient_priority, "oncoplot_patient_priority_scores.csv", row.names = FALSE)

patient_priority_ordered <- patient_priority %>%
  mutate(patient_id = factor(patient_id, levels = patient_order)) %>%
  arrange(patient_id)

write.csv(patient_priority_ordered, "oncoplot_patient_priority_scores_ordered.csv", row.names = FALSE)

# 11) Numeric coordinates
plot_df <- oncodata %>%
  mutate(x = as.numeric(patient_id), y = as.numeric(gene))

annot_df <- patient_annot %>%
  mutate(x = as.numeric(factor(patient_id, levels = levels(oncodata$patient_id))))

n_genes <- length(levels(oncodata$gene))
n_patients <- length(levels(oncodata$patient_id))

y_diag <- 0
y_sex  <- -1

# Diagnosis group separators
group_bounds <- annot_df %>%
  arrange(x) %>%
  group_by(diagnosis) %>%
  summarise(xmin = min(x), xmax = max(x), .groups = "drop") %>%
  arrange(xmin)

sep_x <- head(group_bounds$xmax + 0.5, -1)

# 12) Build plot
p <- ggplot() +
  # Main mutation tiles
  geom_tile(
    data = plot_df,
    aes(x = x, y = y),
    width = 0.95, height = 0.95,
    fill = "grey35", color = "white", linewidth = 0.2
  ) +
  # White diagonal for multiple variants
  geom_segment(
    data = plot_df %>% filter(n_variants > 1),
    aes(x = x - 0.45, xend = x + 0.45, y = y - 0.45, yend = y + 0.45),
    color = "white", linewidth = 0.5
  ) +
  # Diagnosis row
  geom_tile(
    data = annot_df,
    aes(x = x, y = y_diag, fill = diagnosis),
    width = 0.95, height = 0.95,
    color = "white", linewidth = 0.2
  ) +
  # Sex row
  geom_tile(
    data = annot_df,
    aes(x = x, y = y_sex, fill = sex),
    width = 0.95, height = 0.95,
    color = "white", linewidth = 0.2
  ) +

  # Boundary gridlines
  geom_vline(xintercept = seq(0.5, n_patients + 0.5, by = 1), color = "grey88", linewidth = 0.2) +
  geom_hline(yintercept = seq(-1.5, n_genes + 0.5, by = 1), color = "grey88", linewidth = 0.2) +
  
  # Single thick white separator between diagnosis groups
  geom_vline(xintercept = sep_x, color = "white", linewidth = 1) +
  
  scale_fill_manual(
    values = c(
      "Female" = "#8e44ad",
      "Male" = "#27ae60",
      "Myeloid" = "#1f77b4",
      "Lymphoid" = "#f1c40f",
      "No haem cause" = "#9e9e9e"
    )
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_patients + 0.5)) +
  scale_y_continuous(
    breaks = c(y_sex, y_diag, 1:n_genes),
    labels = c("Sex", "Diagnosis", levels(oncodata$gene)),
    expand = c(0, 0),
    limits = c(-1.5, n_genes + 0.5)
  ) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),
    panel.background = element_rect(fill = "grey95", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "plain", size = 11, family = "Arial"),
    legend.position = "none",
    plot.margin = margin(10, 20, 20, 10)
  )

print(p)

# 13) Save outputs
ggsave("oncoplot_VCS_final.png", p, width = 16, height = 10, dpi = 300)
write.csv(oncodata, "Oncoplot_matrix_data.csv", row.names = FALSE)
