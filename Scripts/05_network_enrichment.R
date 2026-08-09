# ==============================================================================
# 05_network_enrichment.R
# Research question 4: what network are IGF1R and AKT3 embedded in, in TNBC?
# Produces Figure 4 (enrichment) and the input file for Figure 5 (Cytoscape)
# ==============================================================================

library(SummarizedExperiment)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

rna_se       <- readRDS("data/processed/tcga_brca_rna_se.rds")
analysis_tbl <- readRDS("data/processed/analysis_table.rds")

tnbc_patients <- analysis_tbl %>%
  filter(subtype_ihc == "TNBC") %>% pull(patient)

# ------------------------------------------------------------------------------
# 1. Build the TNBC-only expression matrix
# ------------------------------------------------------------------------------
tpm <- assay(rna_se, "tpm_unstrand")
colnames_patient <- substr(colnames(tpm), 1, 12)
is_tumor <- as.numeric(substr(colnames(tpm), 14, 15)) < 10

keep <- colnames_patient %in% tnbc_patients & is_tumor
expr <- log2(tpm[, keep, drop = FALSE] + 1)
rownames(expr) <- rowData(rna_se)$gene_name

# Drop low-expressed genes: keep genes detected in at least half the samples
expr <- expr[rowMeans(expr > 1) > 0.5, ]
cat("Genes retained:", nrow(expr), "| TNBC samples:", ncol(expr), "\n")

# ------------------------------------------------------------------------------
# 2. Co-expression with IGF1R and AKT3
# ------------------------------------------------------------------------------
get_coexpressed <- function(gene, rho_cut = 0.4, fdr_cut = 0.05) {
  v <- expr[gene, ]
  res <- apply(expr, 1, function(x) {
    ct <- suppressWarnings(cor.test(x, v, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  })
  tibble(
    gene = colnames(res) %||% rownames(expr),
    rho  = res["rho", ],
    p    = res["p", ]
  ) %>%
    mutate(fdr = p.adjust(p, method = "BH")) %>%
    filter(abs(rho) > rho_cut, fdr < fdr_cut, gene != !!gene) %>%
    arrange(desc(abs(rho)))
}

coexp_igf1r <- get_coexpressed("IGF1R")
coexp_akt3  <- get_coexpressed("AKT3")

cat("Co-expressed with IGF1R:", nrow(coexp_igf1r), "genes\n")
cat("Co-expressed with AKT3 :", nrow(coexp_akt3),  "genes\n")

write_csv(coexp_igf1r, "tables/Table_S3_coexpression_IGF1R.csv")
write_csv(coexp_akt3,  "tables/Table_S4_coexpression_AKT3.csv")

# ------------------------------------------------------------------------------
# 3. GO and KEGG enrichment
# ------------------------------------------------------------------------------
run_enrichment <- function(gene_vec, label) {
  entrez <- bitr(gene_vec, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)$ENTREZID
  universe <- bitr(rownames(expr), fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)$ENTREZID

  # A proper universe (all expressed genes) matters. Using the whole genome as
  # background inflates significance and reviewers notice.
  go <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                 universe = universe, pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
  kegg <- enrichKEGG(entrez, organism = "hsa", universe = universe,
                     pvalueCutoff = 0.05)

  write_csv(as.data.frame(go),   paste0("tables/GO_", label, ".csv"))
  write_csv(as.data.frame(kegg), paste0("tables/KEGG_", label, ".csv"))
  list(go = go, kegg = kegg)
}

enr_igf1r <- run_enrichment(coexp_igf1r$gene, "IGF1R")
enr_akt3  <- run_enrichment(coexp_akt3$gene,  "AKT3")

# --- Figure 4 -----------------------------------------------------------------
pdf("figures/Figure4_enrichment.pdf", width = 10, height = 8)
print(dotplot(enr_igf1r$go, showCategory = 15,
              title = "GO BP - genes co-expressed with IGF1R (TNBC)"))
print(dotplot(enr_akt3$go, showCategory = 15,
              title = "GO BP - genes co-expressed with AKT3 (TNBC)"))
dev.off()

# ------------------------------------------------------------------------------
# 4. Export node list for STRING / Cytoscape (Figure 5)
# ------------------------------------------------------------------------------
top_nodes <- unique(c("IGF1R", "AKT3",
                      head(coexp_igf1r$gene, 50),
                      head(coexp_akt3$gene, 50)))

writeLines(top_nodes, "tables/string_input_genes.txt")

cat("\n--- Next step for Figure 5 ---\n")
cat("1. Paste tables/string_input_genes.txt into https://string-db.org\n")
cat("2. Set confidence to 0.7 (high), export as TSV\n")
cat("3. Import into Cytoscape, run cytoHubba (MCC method) for hub genes\n")
cat("4. Export the figure at 300 dpi into figures/Figure5_network\n")
