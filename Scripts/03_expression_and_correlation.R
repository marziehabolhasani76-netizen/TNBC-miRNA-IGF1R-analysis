# ==============================================================================
# 03_expression_and_correlation.R
# Research questions 1 and 2.
# Produces Figure 1 (differential expression) and Figure 2 (miRNA-mRNA correlation)
# ==============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)

analysis_tbl <- readRDS("data/processed/analysis_table.rds")
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

theme_set(theme_classic(base_size = 11))

# ------------------------------------------------------------------------------
# Q1: Is IGF1R / AKT3 / miR-497 / miR-424 differentially expressed in TNBC?
# ------------------------------------------------------------------------------
targets <- c("IGF1R", "AKT3", "hsa-mir-497", "hsa-mir-424")

de_tbl <- analysis_tbl %>%
  filter(!is.na(subtype_ihc)) %>%
  select(patient, subtype_ihc, all_of(targets)) %>%
  pivot_longer(all_of(targets), names_to = "feature", values_to = "expr")

# Wilcoxon rank-sum, NOT t-test. Expression data are not normally distributed.
de_results <- de_tbl %>%
  group_by(feature) %>%
  summarise(
    median_TNBC     = median(expr[subtype_ihc == "TNBC"], na.rm = TRUE),
    median_nonTNBC  = median(expr[subtype_ihc == "non-TNBC"], na.rm = TRUE),
    log2FC          = median_TNBC - median_nonTNBC,   # already log2 scale
    p_value         = wilcox.test(expr ~ subtype_ihc)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_csv(de_results, "tables/Table_S1_differential_expression.csv")
print(de_results)

# --- Figure 1 -----------------------------------------------------------------
fig1 <- de_tbl %>%
  mutate(feature = factor(feature, levels = targets)) %>%
  ggplot(aes(x = subtype_ihc, y = expr, fill = subtype_ihc)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  facet_wrap(~ feature, scales = "free_y", nrow = 1) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  scale_fill_manual(values = c("TNBC" = "#C0392B", "non-TNBC" = "#5D8AA8")) +
  labs(x = NULL, y = expression(log[2]*"(expression + 1)"), fill = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/Figure1_differential_expression.pdf", fig1,
       width = 9, height = 3.2)
ggsave("figures/Figure1_differential_expression.png", fig1,
       width = 9, height = 3.2, dpi = 300)

# ------------------------------------------------------------------------------
# Q2: Do miR-497 and miR-424 correlate inversely with IGF1R and AKT3?
#     IMPORTANT: correlate WITHIN the TNBC cohort only.
#     Correlating across all breast cancers gives a subtype artefact, not biology.
# ------------------------------------------------------------------------------
tnbc <- analysis_tbl %>% filter(subtype_ihc == "TNBC")
cat("\nTNBC patients used for correlation:", nrow(tnbc), "\n")

pairs_to_test <- tribble(
  ~mirna,          ~gene,    ~type,
  "hsa-mir-497",   "IGF1R",  "hypothesised",
  "hsa-mir-424",   "AKT3",   "hypothesised",
  "hsa-mir-497",   "AKT3",   "cross-check",
  "hsa-mir-424",   "IGF1R",  "cross-check",
  "hsa-mir-16-1",  "IGF1R",  "negative control",
  "hsa-let-7a-1",  "AKT3",   "negative control"
)

cor_results <- pairs_to_test %>%
  rowwise() %>%
  mutate(
    test = list(cor.test(tnbc[[mirna]], tnbc[[gene]], method = "spearman",
                         exact = FALSE)),
    rho     = test$estimate,
    p_value = test$p.value
  ) %>%
  ungroup() %>%
  select(-test) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write_csv(cor_results, "tables/Table_S2_correlations.csv")
print(cor_results)

# --- Figure 2 -----------------------------------------------------------------
make_scatter <- function(mirna, gene) {
  ggplot(tnbc, aes(x = .data[[mirna]], y = .data[[gene]])) +
    geom_point(alpha = 0.55, size = 1.4, colour = "#34495E") +
    geom_smooth(method = "lm", se = TRUE, colour = "#C0392B", linewidth = 0.7) +
    stat_cor(method = "spearman", size = 3.2) +
    labs(x = paste0(mirna, "  (log2 RPM)"),
         y = paste0(gene, "  (log2 TPM)"))
}

fig2 <- make_scatter("hsa-mir-497", "IGF1R") +
        make_scatter("hsa-mir-424", "AKT3") +
        plot_annotation(tag_levels = "A")

ggsave("figures/Figure2_correlation.pdf", fig2, width = 7.5, height = 3.4)
ggsave("figures/Figure2_correlation.png", fig2, width = 7.5, height = 3.4, dpi = 300)

cat("\n--- Figures 1 and 2 saved to figures/ ---\n")

# ==============================================================================
# HOW TO READ YOUR OWN RESULT - read this before you feel disappointed
# ------------------------------------------------------------------------------
# If the inverse correlation is NOT there, or is weak: that is a real finding,
# not a failure. Bulk tumour tissue mixes many cell types, and miRNA repression
# often acts on protein rather than transcript level, so a modest or absent
# mRNA-level correlation is entirely expected and well documented.
#
# Report exactly what you find. Discuss why it might be so. A paper that says
# "the hypothesised inverse relationship was not detectable at transcript level
# in bulk tumours, for these reasons" is publishable and honest.
# A paper that quietly drops the inconvenient panel is neither.
# ==============================================================================
