# ==============================================================================
# 02_preprocess.R
# Build a clean analysis object: expression of the 4 genes/miRNAs of interest,
# a defined TNBC cohort, and survival variables.
#
# THIS IS THE MOST IMPORTANT SCRIPT IN THE PROJECT.
# If the TNBC definition is wrong, everything downstream is wrong.
# ==============================================================================

library(SummarizedExperiment)
library(tidyverse)

rna_se      <- readRDS("data/processed/tcga_brca_rna_se.rds")
mirna_df    <- readRDS("data/processed/tcga_brca_mirna.rds")
clinical    <- readRDS("data/processed/tcga_brca_clinical.rds")
patient_tbl <- readRDS("data/processed/tcga_brca_patient_biotab.rds")

# ------------------------------------------------------------------------------
# 1. Genes and miRNAs of interest
# ------------------------------------------------------------------------------
genes_of_interest  <- c("IGF1R", "AKT3")
mirnas_of_interest <- c("hsa-mir-497", "hsa-mir-424")

# Negative controls: miRNAs with no reported link to this axis.
# Including these is what separates a real analysis from a fishing expedition.
mirna_controls <- c("hsa-mir-16-1", "hsa-let-7a-1")

# ------------------------------------------------------------------------------
# 2. Extract mRNA expression (TPM, log2-transformed)
# ------------------------------------------------------------------------------
tpm <- assay(rna_se, "tpm_unstrand")
gene_symbols <- rowData(rna_se)$gene_name

idx <- which(gene_symbols %in% genes_of_interest)
mrna_mat <- log2(tpm[idx, , drop = FALSE] + 1)
rownames(mrna_mat) <- gene_symbols[idx]

# Collapse duplicate symbols by taking the highest-expressed row
if (any(duplicated(rownames(mrna_mat)))) {
  mrna_mat <- do.call(rbind, lapply(split(seq_len(nrow(mrna_mat)),
                                          rownames(mrna_mat)), function(i) {
    if (length(i) == 1) mrna_mat[i, ] else mrna_mat[i[which.max(rowMeans(mrna_mat[i, ]))], ]
  }))
}

mrna_long <- as.data.frame(t(mrna_mat)) %>%
  rownames_to_column("barcode") %>%
  mutate(
    patient     = substr(barcode, 1, 12),
    sample_type = ifelse(as.numeric(substr(barcode, 14, 15)) < 10,
                         "Tumor", "Normal")
  )

cat("mRNA samples - Tumor:", sum(mrna_long$sample_type == "Tumor"),
    "| Normal:", sum(mrna_long$sample_type == "Normal"), "\n")

# ------------------------------------------------------------------------------
# 3. Extract miRNA expression (reads_per_million, log2-transformed)
# ------------------------------------------------------------------------------
rpm_cols <- grep("^reads_per_million", names(mirna_df), value = TRUE)

mirna_mat <- mirna_df %>%
  filter(miRNA_ID %in% c(mirnas_of_interest, mirna_controls)) %>%
  column_to_rownames("miRNA_ID") %>%
  select(all_of(rpm_cols)) %>%
  as.matrix()

mirna_mat <- log2(mirna_mat + 1)
colnames(mirna_mat) <- sub("^reads_per_million_miRNA_mapped_", "", colnames(mirna_mat))

mirna_long <- as.data.frame(t(mirna_mat)) %>%
  rownames_to_column("barcode") %>%
  mutate(
    patient     = substr(barcode, 1, 12),
    sample_type = ifelse(as.numeric(substr(barcode, 14, 15)) < 10,
                         "Tumor", "Normal")
  )

# ------------------------------------------------------------------------------
# 4. Define the TNBC cohort - TWO independent definitions
# ------------------------------------------------------------------------------
# Definition A: IHC receptor status (the clinical definition)
# The exact column names shift between TCGA releases. Check them before trusting.
receptor_cols <- grep("estrogen|progesterone|her2", names(patient_tbl),
                      ignore.case = TRUE, value = TRUE)
cat("\nReceptor columns found:\n"); print(receptor_cols)

# EDIT THESE THREE LINES to match the real column names printed above.
er_col   <- "er_status_by_ihc"
pr_col   <- "pr_status_by_ihc"
her2_col <- "her2_status_by_ihc"

subtype_tbl <- patient_tbl %>%
  select(patient = bcr_patient_barcode,
         ER = all_of(er_col), PR = all_of(pr_col), HER2 = all_of(her2_col)) %>%
  filter(!patient %in% c("bcr_patient_barcode", "CDE_ID:2003301")) %>%  # header rows
  mutate(
    subtype_ihc = case_when(
      ER == "Negative" & PR == "Negative" & HER2 == "Negative" ~ "TNBC",
      ER %in% c("Positive", "Negative") &
        PR %in% c("Positive", "Negative") &
        HER2 %in% c("Positive", "Negative")                    ~ "non-TNBC",
      TRUE                                                     ~ NA_character_
    )
  )

cat("\nTNBC cohort by IHC:\n")
print(table(subtype_tbl$subtype_ihc, useNA = "ifany"))

# Definition B: PAM50 Basal-like (the molecular definition)
# Use TCGAbiolinks::PanCancerAtlas_subtypes() or TCGAquery_subtype("BRCA").
# Report results under BOTH definitions. Reviewers always ask.

# ------------------------------------------------------------------------------
# 5. Survival variables
# ------------------------------------------------------------------------------
surv_tbl <- clinical %>%
  transmute(
    patient   = submitter_id,
    os_status = ifelse(vital_status == "Dead", 1, 0),
    os_time   = ifelse(vital_status == "Dead",
                       days_to_death, days_to_last_follow_up),
    age       = age_at_index,
    stage     = ajcc_pathologic_stage
  ) %>%
  filter(!is.na(os_time), os_time > 0) %>%
  mutate(
    os_months  = os_time / 30.44,
    stage_simple = case_when(
      str_detect(stage, "Stage IV")              ~ "III-IV",
      str_detect(stage, "Stage III")             ~ "III-IV",
      str_detect(stage, "Stage II")              ~ "II",
      str_detect(stage, "Stage I")               ~ "I",
      TRUE                                       ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# 6. Merge everything into one analysis table
# ------------------------------------------------------------------------------
analysis_tbl <- mrna_long %>%
  filter(sample_type == "Tumor") %>%
  distinct(patient, .keep_all = TRUE) %>%
  inner_join(mirna_long %>% filter(sample_type == "Tumor") %>%
               distinct(patient, .keep_all = TRUE) %>%
               select(-barcode, -sample_type),
             by = "patient") %>%
  left_join(subtype_tbl, by = "patient") %>%
  left_join(surv_tbl,    by = "patient")

saveRDS(analysis_tbl, "data/processed/analysis_table.rds")

cat("\n--- Final analysis table ---\n")
cat("Patients with matched mRNA + miRNA:", nrow(analysis_tbl), "\n")
cat("TNBC patients:", sum(analysis_tbl$subtype_ihc == "TNBC", na.rm = TRUE), "\n")
cat("TNBC with survival data:",
    sum(analysis_tbl$subtype_ihc == "TNBC" & !is.na(analysis_tbl$os_time),
        na.rm = TRUE), "\n")

# ==============================================================================
# SANITY CHECKS - do not skip these
# ------------------------------------------------------------------------------
# 1. Is the TNBC count roughly 150-200? If it is 12 or 900, the receptor
#    column names are wrong.
# 2. Do the expression values look like log2 values (roughly 0-15)?
# 3. Does every patient appear exactly once?
# Write the answers in your lab notebook. Week 10 you will need them.
# ==============================================================================

