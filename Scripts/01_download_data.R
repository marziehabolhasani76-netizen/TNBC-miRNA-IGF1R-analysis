# ==============================================================================
# 01_download_data.R
# Download TCGA-BRCA RNA-seq, miRNA-seq and clinical data.
# Expect this to take 30-90 minutes and a few GB of disk. Run it once,
# save the result, and never download again.
# ==============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. mRNA expression (STAR - Counts)
# ------------------------------------------------------------------------------
query_rna <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_rna, directory = "data/raw", files.per.chunk = 50)
rna_se <- GDCprepare(query_rna, directory = "data/raw")
saveRDS(rna_se, "data/processed/tcga_brca_rna_se.rds")

cat("mRNA data:", dim(rna_se)[1], "genes x", dim(rna_se)[2], "samples\n")

# ------------------------------------------------------------------------------
# 2. miRNA expression
# ------------------------------------------------------------------------------
query_mirna <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "miRNA Expression Quantification"
)

GDCdownload(query_mirna, directory = "data/raw", files.per.chunk = 50)
mirna_df <- GDCprepare(query_mirna, directory = "data/raw")
saveRDS(mirna_df, "data/processed/tcga_brca_mirna.rds")

cat("miRNA data:", dim(mirna_df)[1], "miRNAs x", ncol(mirna_df), "columns\n")

# ------------------------------------------------------------------------------
# 3. Clinical data (indexed - includes vital status and follow-up time)
# ------------------------------------------------------------------------------
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
saveRDS(clinical, "data/processed/tcga_brca_clinical.rds")

cat("Clinical data:", nrow(clinical), "patients\n")

# ------------------------------------------------------------------------------
# 4. Receptor status (ER / PR / HER2) - needed to define TNBC
# ------------------------------------------------------------------------------
# The IHC receptor calls live in the BCR biotab supplementary clinical files,
# not in the indexed clinical table above.
query_biotab <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)
GDCdownload(query_biotab, directory = "data/raw")
biotab <- GDCprepare(query_biotab, directory = "data/raw")

# The patient table is usually named like "clinical_patient_brca"
patient_tbl <- biotab[[grep("patient", names(biotab))[1]]]
saveRDS(patient_tbl, "data/processed/tcga_brca_patient_biotab.rds")

cat("\nColumn names containing 'receptor' or 'her2':\n")
print(grep("receptor|her2|estrogen|progesterone", names(patient_tbl),
           ignore.case = TRUE, value = TRUE))

# ==============================================================================
# FALLBACK: manual download from UCSC Xena
# ------------------------------------------------------------------------------
# If GDC downloads keep failing, go to https://xenabrowser.net, choose the
# GDC TCGA Breast Cancer (BRCA) cohort, and download these as .tsv.gz:
#   - HTSeq/STAR gene expression (log2(count+1) or TPM)
#   - miRNA expression
#   - phenotype / clinical matrix
#   - survival data
# Then read them in with data.table::fread() and skip straight to 02_preprocess.R.
# Same analysis, same conclusions - only the download route differs.
# ==============================================================================

cat("\n--- Download step complete. Move on to 02_preprocess.R ---\n")
