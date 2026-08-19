# =========================================================
# Table: Alpha diversity summary with mixed-model significance
# (Shannon, Pielou's evenness, Faith's PD)
#
# UPDATED after diagnostic review:
#   - Shannon diversity: Gaussian LMM (diagnostics were fine)
#       Shannon ~ HostType * Season + (1 | Site)
#   - Faith's PD: Gaussian LMM on log-transformed response
#     (raw-scale diagnostics showed heteroscedasticity + right
#     skew; Faith's PD is strictly positive and naturally skewed)
#       log(Faith_PD) ~ HostType * Season + (1 | Site)
#   - Pielou's evenness: BETA regression mixed model, NOT
#     Gaussian LMM (raw diagnostics showed strong non-normality
#     -- evenness is a proportion bounded in (0,1) and values
#     cluster near the upper bound, which a Gaussian LMM cannot
#     handle correctly)
#       Evenness ~ HostType * Season + (1 | Site), family = beta_family()
#     Fit with glmmTMB (supports beta family + random effects +
#     works cleanly with emmeans for post-hoc comparisons).
#
# All three models: Site as random intercept (repeated sampling
# of the same locations across seasons), HostType * Season as
# fixed effects, Type II Wald chi-square tests for fixed-effect
# significance, and full pairwise post-hoc comparisons via
# emmeans with Benjamini-Hochberg (FDR) correction.
#
# Requires: results/phylo_raw.rds
# Packages: lme4, lmerTest, car, emmeans, glmmTMB (new)
# =========================================================

library(phyloseq)
library(microbiome)
library(picante)
library(ape)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(glmmTMB)

results_dir <- "results"
phylo_raw <- readRDS(file.path(results_dir, "phylo_raw.rds"))

stopifnot(all(c("HostType", "Site") %in% colnames(sample_data(phylo_raw))))

season_levels <- c("Spring", "Summer", "Autumn", "Winter")
host_levels <- c("Water sample", "Arbacia lixula", "Paracentrotus lividus")

# ---------------------------------------------------------
# Rarefy at 5000 reads (identical to Figure 3)
# ---------------------------------------------------------
set.seed(123)
rare_depth <- 5000

phylo_rarefied <- rarefy_even_depth(
  phylo_raw, sample.size = rare_depth, rngseed = 123,
  replace = FALSE, trimOTUs = TRUE, verbose = TRUE
)

shannon_vals  <- microbiome::diversity(phylo_rarefied, index = "shannon")
richness_vals <- microbiome::richness(phylo_rarefied, index = "observed")

shannon_named  <- setNames(shannon_vals$shannon, rownames(shannon_vals))
richness_named <- setNames(richness_vals$observed, rownames(richness_vals))
evenness_named <- shannon_named / log(richness_named)

otu_mat <- as(otu_table(phylo_rarefied), "matrix")
if (taxa_are_rows(phylo_rarefied)) otu_mat <- t(otu_mat)

tree <- phy_tree(phylo_rarefied)
common_taxa <- intersect(colnames(otu_mat), tree$tip.label)
otu_mat_pd <- otu_mat[, common_taxa, drop = FALSE]
tree_pd <- ape::keep.tip(tree, common_taxa)

faith_pd <- picante::pd(otu_mat_pd, tree_pd, include.root = TRUE)

meta_df <- data.frame(sample_data(phylo_rarefied), check.names = FALSE, stringsAsFactors = FALSE)
meta_df$Sample <- rownames(meta_df)

div_df <- data.frame(
  Sample = rownames(otu_mat),
  Shannon = shannon_named[rownames(otu_mat)],
  Evenness = evenness_named[rownames(otu_mat)],
  Faith_PD = faith_pd$PD[match(rownames(otu_mat), rownames(faith_pd))]
) %>%
  left_join(meta_df, by = "Sample") %>%
  mutate(
    Season = factor(Season, levels = season_levels),
    HostType = factor(HostType, levels = host_levels),
    Site = factor(Site),
    log_Faith_PD = log(Faith_PD),
    # glmmTMB's beta family requires values strictly inside (0,1);
    # evenness of exactly 1 (perfectly even, richness-limited edge
    # case) or 0 would break the model, so nudge any boundary
    # values inward by a tiny epsilon (standard practice, e.g.
    # Smithson & Verkuilen 2006).
    Evenness_beta = pmin(pmax(Evenness, 1e-4), 1 - 1e-4)
  )

metric_labels <- c(Shannon = "Shannon diversity", Evenness = "Pielou's evenness", Faith_PD = "Faith's PD")

sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# ---------------------------------------------------------
# Descriptive stats: mean +/- SD per HostType (RAW scale,
# not transformed -- these are purely descriptive, unaffected
# by which model is used for significance testing)
# ---------------------------------------------------------
desc_df <- div_df %>%
  tidyr::pivot_longer(cols = c(Shannon, Evenness, Faith_PD), names_to = "Metric", values_to = "Value") %>%
  dplyr::group_by(Metric, HostType) %>%
  dplyr::summarise(
    mean_sd = sprintf("%.2f \u00b1 %.2f", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = HostType, values_from = mean_sd)

# ---------------------------------------------------------
# Helper: extract fixed-effect + pairwise results from any
# model object that emmeans/car::Anova can handle (lmer or
# glmmTMB), given the metric label to attach
# ---------------------------------------------------------
extract_results <- function(mod, metric_label, type2_test = "Chisq") {

  aov_tab <- car::Anova(mod, type = "II", test.statistic = type2_test)

  fixed_p <- data.frame(
    Metric = metric_label,
    HostType_p = aov_tab["HostType", grep("^Pr", colnames(aov_tab))],
    Season_p = aov_tab["Season", grep("^Pr", colnames(aov_tab))],
    Interaction_p = aov_tab["HostType:Season", grep("^Pr", colnames(aov_tab))]
  )

  emm_host <- emmeans(mod, pairwise ~ HostType, adjust = "fdr", type = "response")
  emm_season <- emmeans(mod, pairwise ~ Season, adjust = "fdr", type = "response")

  host_pairs <- as.data.frame(emm_host$contrasts) %>%
    dplyr::mutate(Metric = metric_label, Comparison = contrast) %>%
    dplyr::select(Metric, Comparison, p.value)

  season_pairs <- as.data.frame(emm_season$contrasts) %>%
    dplyr::mutate(Metric = metric_label, Comparison = contrast) %>%
    dplyr::select(Metric, Comparison, p.value)

  list(fixed_p = fixed_p, host_pairs = host_pairs, season_pairs = season_pairs, model = mod)
}

# ---------------------------------------------------------
# Shannon: Gaussian LMM (diagnostics OK)
# ---------------------------------------------------------
mod_shannon <- lmerTest::lmer(Shannon ~ HostType * Season + (1 | Site), data = div_df, REML = TRUE)
res_shannon <- extract_results(mod_shannon, "Shannon")

# ---------------------------------------------------------
# Faith's PD: Gaussian LMM on log-transformed response
# ---------------------------------------------------------
mod_faith <- lmerTest::lmer(log_Faith_PD ~ HostType * Season + (1 | Site), data = div_df, REML = TRUE)
res_faith <- extract_results(mod_faith, "Faith_PD")
# NOTE: emmeans contrasts above are on the log scale (differences
# in log(Faith_PD)); p-values are unaffected by back-transformation,
# but if you report effect sizes, back-transform as needed.

# ---------------------------------------------------------
# Evenness: Beta regression mixed model (glmmTMB)
# ---------------------------------------------------------
mod_evenness <- glmmTMB(
  Evenness_beta ~ HostType * Season + (1 | Site),
  family = beta_family(link = "logit"),
  data = div_df
)
res_evenness <- extract_results(mod_evenness, "Evenness", type2_test = "Chisq")

results_list <- list(Shannon = res_shannon, Faith_PD = res_faith, Evenness = res_evenness)

fixed_p_df <- dplyr::bind_rows(lapply(results_list, `[[`, "fixed_p"))
host_pairs_df <- dplyr::bind_rows(lapply(results_list, `[[`, "host_pairs"))
season_pairs_df <- dplyr::bind_rows(lapply(results_list, `[[`, "season_pairs"))

# ---------------------------------------------------------
# Format fixed-effect p-values with significance stars
# ---------------------------------------------------------
fixed_p_fmt <- fixed_p_df %>%
  dplyr::mutate(
    `HostType (p)` = sprintf("%.3f %s", HostType_p, sig_stars(HostType_p)),
    `Season (p)` = sprintf("%.3f %s", Season_p, sig_stars(Season_p)),
    `HostType x Season (p)` = sprintf("%.3f %s", Interaction_p, sig_stars(Interaction_p))
  ) %>%
  dplyr::select(Metric, `HostType (p)`, `Season (p)`, `HostType x Season (p)`)

# ---------------------------------------------------------
# Pairwise tables (wide format, one column per comparison)
# ---------------------------------------------------------
host_pairs_wide <- host_pairs_df %>%
  dplyr::mutate(p.value = sprintf("%.3f", p.value)) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = p.value) %>%
  dplyr::rename_with(~ paste0("Host: ", .x), -Metric)

season_pairs_wide <- season_pairs_df %>%
  dplyr::mutate(p.value = sprintf("%.3f", p.value)) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = p.value) %>%
  dplyr::rename_with(~ paste0("Season: ", .x), -Metric)

# ---------------------------------------------------------
# Combine everything into one master table
# ---------------------------------------------------------
table_diversity <- desc_df %>%
  dplyr::left_join(fixed_p_fmt, by = "Metric") %>%
  dplyr::left_join(host_pairs_wide, by = "Metric") %>%
  dplyr::left_join(season_pairs_wide, by = "Metric") %>%
  dplyr::mutate(Metric = metric_labels[Metric]) %>%
  dplyr::rename(
    `Water sample (mean ± SD)` = `Water sample`,
    `A. lixula (mean ± SD)` = `Arbacia lixula`,
    `P. lividus (mean ± SD)` = `Paracentrotus lividus`
  )

# ---------------------------------------------------------
# Export
# ---------------------------------------------------------
dir.create(results_dir, showWarnings = FALSE)
write.csv(table_diversity, file.path(results_dir, "TableSXX_alpha_diversity_summary.csv"), row.names = FALSE)
print(table_diversity)

cat("\nModels used:\n")
cat("  Shannon diversity : Shannon ~ HostType * Season + (1 | Site)              [Gaussian LMM]\n")
cat("  Faith's PD        : log(Faith_PD) ~ HostType * Season + (1 | Site)        [Gaussian LMM, log-transformed]\n")
cat("  Pielou's evenness : Evenness ~ HostType * Season + (1 | Site)             [Beta regression, glmmTMB]\n")
cat("Fixed-effect p-values: car::Anova(type = 'II', test.statistic = 'Chisq').\n")
cat("Pairwise post-hoc: emmeans, Benjamini-Hochberg (FDR)-adjusted p-values.\n")
cat("Significance codes: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant.\n")

# ---------------------------------------------------------
# Diagnostics -- re-check all three models with the corrected
# specifications (Faith's PD now log-transformed, Evenness now
# beta regression)
# ---------------------------------------------------------
dir.create("figures/diagnostics", showWarnings = FALSE, recursive = TRUE)

# Shannon and Faith's PD (lmer): base-R residual plots as before
for (nm in c("Shannon", "Faith_PD")) {
  mod <- results_list[[nm]]$model
  png(file.path("figures/diagnostics", paste0("diagnostics_", nm, "_v2.png")), width = 900, height = 450)
  par(mfrow = c(1, 2))
  plot(fitted(mod), resid(mod),
       xlab = "Fitted", ylab = "Residuals", main = paste(nm, "- Residuals vs Fitted (v2)"))
  abline(h = 0, lty = 2)
  qqnorm(resid(mod), main = paste(nm, "- QQ plot (v2)"))
  qqline(resid(mod))
  dev.off()
}

# Evenness (glmmTMB beta): use DHARMa for correct simulated
# residuals, since raw Pearson residuals from a beta GLMM are
# not expected to look Gaussian and base-R qqnorm would be
# misleading here.
if (requireNamespace("DHARMa", quietly = TRUE)) {
  library(DHARMa)
  sim_res <- simulateResiduals(mod_evenness)
  png("figures/diagnostics/diagnostics_Evenness_v2_DHARMa.png", width = 900, height = 450)
  plot(sim_res)
  dev.off()
  cat("\nEvenness diagnostics (DHARMa simulated residuals) saved.\n")
} else {
  cat("\nPackage 'DHARMa' not installed -- install it (install.packages('DHARMa'))\n")
  cat("to properly check the beta-regression Evenness model; base-R qqnorm()\n")
  cat("is not appropriate for non-Gaussian GLMM residuals.\n")
}

cat("\nCheck figures/diagnostics/*_v2* before finalising p-values.\n")
