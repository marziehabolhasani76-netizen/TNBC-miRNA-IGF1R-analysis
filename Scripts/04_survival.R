# ==============================================================================
# 04_survival.R
# Research question 3: does expression predict survival in TNBC?
# Produces Figure 3 (Kaplan-Meier) and Table 1 (multivariable Cox)
# ==============================================================================

library(tidyverse)
library(survival)
library(survminer)

analysis_tbl <- readRDS("data/processed/analysis_table.rds")

surv_data <- analysis_tbl %>%
  filter(subtype_ihc == "TNBC", !is.na(os_time), os_time > 0)

cat("TNBC patients with survival data:", nrow(surv_data), "\n")
cat("Number of death events:", sum(surv_data$os_status), "\n\n")

# ------------------------------------------------------------------------------
# POWER CHECK - do this before you interpret anything
# ------------------------------------------------------------------------------
# TCGA TNBC has few deaths (often under 30). With that few events, survival
# analysis is badly underpowered and a non-significant result means almost
# nothing. Report the event count honestly in the manuscript, and lean on
# METABRIC (many more events) as the primary survival cohort if TCGA is thin.
if (sum(surv_data$os_status) < 20) {
  warning("Fewer than 20 events. Survival analysis here is exploratory only. ",
          "Say so explicitly in the paper.")
}

targets <- c("IGF1R", "AKT3", "hsa-mir-497", "hsa-mir-424")

# ------------------------------------------------------------------------------
# 1. Kaplan-Meier, median split (PRIMARY analysis)
# ------------------------------------------------------------------------------
# The median is chosen a priori. This matters: optimising the cutpoint and then
# reporting the resulting p-value is one of the most common reasons these
# papers get rejected.
km_data <- surv_data %>%
  mutate(across(all_of(targets),
                ~ factor(ifelse(.x > median(.x, na.rm = TRUE), "High", "Low"),
                         levels = c("Low", "High")),
                .names = "{.col}_grp"))

km_plots <- lapply(targets, function(g) {
  grp <- paste0(g, "_grp")
  fit <- survfit(as.formula(paste0("Surv(os_months, os_status) ~ ", "`", grp, "`")),
                 data = km_data)
  ggsurvplot(fit, data = km_data,
             pval = TRUE, pval.size = 3.5,
             risk.table = TRUE, risk.table.height = 0.28,
             legend.title = g, legend.labs = c("Low", "High"),
             palette = c("#5D8AA8", "#C0392B"),
             xlab = "Time (months)", ylab = "Overall survival",
             ggtheme = theme_classic(base_size = 10))
})

pdf("figures/Figure3_survival.pdf", width = 10, height = 9)
print(arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2, print = FALSE))
dev.off()

# ------------------------------------------------------------------------------
# 2. Univariable and multivariable Cox regression
# ------------------------------------------------------------------------------
# The multivariable model adjusted for age and stage is what lifts this above
# a GEPIA2 screenshot. Do not skip it.
cox_results <- map_dfr(targets, function(g) {
  d <- surv_data %>% rename(expr = all_of(g))

  uni <- coxph(Surv(os_months, os_status) ~ expr, data = d)
  s_uni <- summary(uni)

  multi <- coxph(Surv(os_months, os_status) ~ expr + age + stage_simple, data = d)
  s_multi <- summary(multi)

  tibble(
    feature      = g,
    HR_uni       = s_uni$coefficients["expr", "exp(coef)"],
    CI_uni       = paste0(round(s_uni$conf.int["expr", "lower .95"], 2), "-",
                          round(s_uni$conf.int["expr", "upper .95"], 2)),
    p_uni        = s_uni$coefficients["expr", "Pr(>|z|)"],
    HR_multi     = s_multi$coefficients["expr", "exp(coef)"],
    CI_multi     = paste0(round(s_multi$conf.int["expr", "lower .95"], 2), "-",
                          round(s_multi$conf.int["expr", "upper .95"], 2)),
    p_multi      = s_multi$coefficients["expr", "Pr(>|z|)"]
  )
}) %>%
  mutate(p_multi_adj = p.adjust(p_multi, method = "BH"))

write_csv(cox_results, "tables/Table1_cox_regression.csv")
print(cox_results)

# ------------------------------------------------------------------------------
# 3. Proportional hazards assumption - reviewers check this
# ------------------------------------------------------------------------------
cat("\n--- Proportional hazards test (Schoenfeld residuals) ---\n")
for (g in targets) {
  d <- surv_data %>% rename(expr = all_of(g))
  m <- coxph(Surv(os_months, os_status) ~ expr + age + stage_simple, data = d)
  cat("\n", g, ":\n"); print(cox.zph(m))
}
# A global p < 0.05 means the PH assumption is violated. Say so and use a
# time-stratified model or restricted mean survival time instead.

cat("\n--- Figure 3 and Table 1 saved ---\n")
