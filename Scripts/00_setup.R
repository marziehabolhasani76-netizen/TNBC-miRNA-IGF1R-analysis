# ==============================================================================
# 00_setup.R
# Install and load every package the project needs. Run this ONCE.
# Project: miR-497/miR-424 - IGF1R/AKT3 axis in TNBC
# ==============================================================================

# --- CRAN packages ------------------------------------------------------------
cran_pkgs <- c(
  "tidyverse",    # data wrangling + ggplot2
  "data.table",
  "survival",     # Cox, Kaplan-Meier
  "survminer",    # KM plots
  "maxstat",      # optimal cutpoint (use with caution - see notes)
  "ggpubr",       # publication-ready plots + stat comparisons
  "patchwork",    # combining panels into one figure
  "here"          # sane file paths
)

new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran)) install.packages(new_cran)

# --- Bioconductor packages ----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "TCGAbiolinks",         # download TCGA data
  "SummarizedExperiment",
  "DESeq2",
  "limma",
  "edgeR",
  "clusterProfiler",      # GO / KEGG enrichment
  "org.Hs.eg.db",         # human gene annotation
  "enrichplot",
  "ComplexHeatmap"
)

new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc)) BiocManager::install(new_bioc, ask = FALSE, update = FALSE)

# --- Load ---------------------------------------------------------------------
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

cat("\n--- Setup complete ---\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

# --- Record versions for the Methods section ----------------------------------
# Journals ask for software versions. Save them now, thank yourself in week 10.
writeLines(capture.output(sessionInfo()), "session_info.txt")
cat("Software versions saved to session_info.txt\n")

# ==============================================================================
# NOTE ON INSTALLATION PROBLEMS
# ------------------------------------------------------------------------------
# If Bioconductor installs fail from Iran, try:
#   1. A different CRAN mirror:
#        options(repos = c(CRAN = "https://cloud.r-project.org"))
#   2. Installing packages one at a time to see which one actually fails.
#   3. On Linux, missing system libraries are a common cause. Errors mentioning
#      libcurl, libxml2, or openssl mean you need the system dev packages first.
# TCGAbiolinks is the heaviest dependency. If it refuses to install, you can
# still do the whole project by downloading TCGA data manually from UCSC Xena
# (https://xenabrowser.net) as plain .tsv files - see 01_download_data.R.
# ==============================================================================

