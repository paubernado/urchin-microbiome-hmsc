# =========================================================
# NULL MODELS (intercept-only)
# MCMC settings are now read from R/00_mcmc_config.R, the
# single source of truth shared with 04_hmsc_models_with_phylogeny.R
# =========================================================

library(Hmsc)
library(coda)

source("R/00_mcmc_config.R")

results_dir <- "results"

stopifnot(exists("HMSC_PA_with_phylogeny_logDepth_siteFixed"))
stopifnot(exists("HMSC_hurdle_abundance_with_phylogeny_parallel"))

# ---------------------------------------------------------
# 1) NULL MODEL — PRESENCE/ABSENCE
# ---------------------------------------------------------

fitted_pa <- HMSC_PA_with_phylogeny_logDepth_siteFixed

m_null_pa <- Hmsc(
  Y = fitted_pa$Y,
  XData = fitted_pa$XData,
  XFormula = ~1,
  phyloTree = fitted_pa$phyloTree,
  distr = "probit"
)

m_null_pa <- sampleMcmc(
  m_null_pa,
  samples   = mcmc_config$pa$samples,
  thin      = mcmc_config$pa$thin,
  transient = mcmc_config$pa$transient,
  nChains   = mcmc_config$pa$nChains,
  nParallel = mcmc_config$n_parallel,
  useSocket = FALSE,
  verbose = 500
)

mpost_null_pa <- convertToCodaObject(m_null_pa)

cat("Null PA model - PSRF Beta:\n")
print(gelman.diag(mpost_null_pa$Beta, multivariate = TRUE)$mpsrf)
cat("Null PA model - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_null_pa$Beta)))

preds_null_pa <- computePredictedValues(m_null_pa, nParallel = mcmc_config$n_parallel)
MF_null_pa <- evaluateModelFit(m_null_pa, preds_null_pa)
cat("Null PA model - Mean AUC:\n")
print(mean(MF_null_pa$AUC, na.rm = TRUE))
cat("Null PA model - Mean Tjur R2:\n")
print(mean(MF_null_pa$TjurR2, na.rm = TRUE))

saveRDS(m_null_pa, file.path(results_dir, "HMSC_null_PA_with_phylogeny.rds"))
saveRDS(mpost_null_pa, file.path(results_dir, "mpost_null_PA_with_phylogeny.rds"))
saveRDS(preds_null_pa, file.path(results_dir, "preds_null_PA_with_phylogeny.rds"))
saveRDS(MF_null_pa, file.path(results_dir, "MF_null_PA_with_phylogeny.rds"))

# Cross-validated version (4-fold, matching Table SXX structure)
set.seed(123)
partition_null_pa <- createPartition(m_null_pa, nfolds = 4)
preds_null_pa_cv <- computePredictedValues(
  m_null_pa, partition = partition_null_pa, nParallel = min(mcmc_config$n_parallel, 5)
)
MF_null_pa_cv <- evaluateModelFit(hM = m_null_pa, predY = preds_null_pa_cv)
cat("Null PA model CROSS-VALIDATED - Mean AUC:\n")
print(mean(MF_null_pa_cv$AUC, na.rm = TRUE))
cat("Null PA model CROSS-VALIDATED - Mean Tjur R2:\n")
print(mean(MF_null_pa_cv$TjurR2, na.rm = TRUE))
saveRDS(MF_null_pa_cv, file.path(results_dir, "MF_null_PA_with_phylogeny_cv.rds"))

# ---------------------------------------------------------
# 2) NULL MODEL — CONDITIONAL ABUNDANCE
# ---------------------------------------------------------

fitted_abund <- HMSC_hurdle_abundance_with_phylogeny_parallel

m_null_abund <- Hmsc(
  Y = fitted_abund$Y,
  XData = fitted_abund$XData,
  XFormula = ~1,
  phyloTree = fitted_abund$phyloTree,
  distr = "normal"
)

m_null_abund <- sampleMcmc(
  m_null_abund,
  samples   = mcmc_config$abund$samples,
  thin      = mcmc_config$abund$thin,
  transient = mcmc_config$abund$transient,
  nChains   = mcmc_config$abund$nChains,
  nParallel = mcmc_config$n_parallel,
  useSocket = FALSE,
  verbose = 500
)

mpost_null_abund <- convertToCodaObject(m_null_abund)

cat("Null abundance model - PSRF Beta:\n")
print(gelman.diag(mpost_null_abund$Beta, multivariate = TRUE)$mpsrf)
cat("Null abundance model - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_null_abund$Beta)))

preds_null_abund <- computePredictedValues(m_null_abund, nParallel = mcmc_config$n_parallel)
MF_null_abund <- evaluateModelFit(m_null_abund, preds_null_abund)
cat("Null abundance model - Mean R2:\n")
print(mean(MF_null_abund$R2, na.rm = TRUE))

saveRDS(m_null_abund, file.path(results_dir, "HMSC_null_abundance_with_phylogeny.rds"))
saveRDS(mpost_null_abund, file.path(results_dir, "mpost_null_abundance_with_phylogeny.rds"))
saveRDS(preds_null_abund, file.path(results_dir, "preds_null_abundance_with_phylogeny.rds"))
saveRDS(MF_null_abund, file.path(results_dir, "MF_null_abundance_with_phylogeny.rds"))

# Cross-validated version
set.seed(123)
partition_null_abund <- createPartition(m_null_abund, nfolds = 4)
preds_null_abund_cv <- computePredictedValues(
  m_null_abund, partition = partition_null_abund, nParallel = min(mcmc_config$n_parallel, 5)
)
MF_null_abund_cv <- evaluateModelFit(hM = m_null_abund, predY = preds_null_abund_cv)
cat("Null abundance model CROSS-VALIDATED - Mean R2:\n")
print(mean(MF_null_abund_cv$R2, na.rm = TRUE))
saveRDS(MF_null_abund_cv, file.path(results_dir, "MF_null_abundance_with_phylogeny_cv.rds"))

# ---------------------------------------------------------
# 3) MODEL vs NULL COMPARISON TABLE (matching Table SXX)
# ---------------------------------------------------------

comparison_table <- data.frame(
  Model = c(
    "PA (full: HostType * Season + Site)",
    "PA (null: intercept only)",
    "Abundance (full: HostType * Season + Site)",
    "Abundance (null: intercept only)"
  ),
  Mean_AUC = c(
    mean(MF_pa_phy$AUC, na.rm = TRUE),
    mean(MF_null_pa$AUC, na.rm = TRUE),
    NA,
    NA
  ),
  Mean_AUC_CV = c(
    mean(MF_pa_phy_cv$AUC, na.rm = TRUE),
    mean(MF_null_pa_cv$AUC, na.rm = TRUE),
    NA,
    NA
  ),
  Mean_R2 = c(
    NA,
    NA,
    mean(MF_abundance_with_phylogeny$R2, na.rm = TRUE),
    mean(MF_null_abund$R2, na.rm = TRUE)
  ),
  Mean_R2_CV = c(
    NA,
    NA,
    mean(MF_abundance_with_phylogeny_cv$R2, na.rm = TRUE),
    mean(MF_null_abund_cv$R2, na.rm = TRUE)
  )
)

print(comparison_table)
write.csv(comparison_table, file.path(results_dir, "TableSXX_model_vs_null_comparison.csv"), row.names = FALSE)

cat("\nNull models complete. Comparison table saved to TableSXX_model_vs_null_comparison.csv\n")
