# =========================================================
# HMSC MODELS WITHOUT PHYLOGENY
# (corresponds to Blocks 6-9 of the original analysis notebook)
#
# NOTE ON MCMC SETTINGS: only the *phylogeny* PA and abundance
# models were independently confirmed from fitted objects (see
# R/00_mcmc_config.R). The no-phylogeny models below use the same
# settings for consistency across the four HMSC models reported
# in the paper. If your no-phylogeny fitted objects used different
# settings, update mcmc_config in R/00_mcmc_config.R accordingly
# rather than editing the numbers here.
#
# NOTE ON logDepth_z: the confirmed production scripts only cover
# the *phylogeny* PA and abundance models, both of which use
# XFormula = ~ HostType * Season + Site + logDepth_z (see
# R/04_hmsc_models_with_phylogeny.R). The no-phylogeny models below
# include the same sequencing-depth covariate for consistency, but
# this has not been independently confirmed from a fitted
# no-phylogeny object. If your no-phylogeny models were fit with a
# different formula, update XFormula here to match.
# =========================================================

library(Hmsc)
library(coda)

source("R/00_mcmc_config.R")

results_dir <- "results"

Y_counts <- readRDS(file.path(results_dir, "Y_counts.rds"))
XData    <- readRDS(file.path(results_dir, "XData.rds"))

# ---------------------------------------------------------
# BLOCK 6. Presence/absence, no phylogeny
# ---------------------------------------------------------
Y_pa <- 1 * (Y_counts > 0)

m_pa <- Hmsc(
  Y = Y_pa,
  XData = XData,
  XFormula = ~ HostType * Season + Site + logDepth_z,
  distr = "probit"
)

m_pa <- sampleMcmc(
  m_pa,
  samples   = mcmc_config$pa$samples,
  thin      = mcmc_config$pa$thin,
  transient = mcmc_config$pa$transient,
  nChains   = mcmc_config$pa$nChains,
  nParallel = mcmc_config$n_parallel,
  verbose   = 500
)

mpost_pa <- convertToCodaObject(m_pa)

cat("\nPA model without phylogeny - PSRF Beta:\n")
print(gelman.diag(mpost_pa$Beta, multivariate = TRUE)$mpsrf)
cat("\nPA model without phylogeny - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_pa$Beta)))

preds_pa <- computePredictedValues(m_pa, nParallel = mcmc_config$n_parallel)
MF_pa <- evaluateModelFit(m_pa, preds_pa)

cat("\nPA model without phylogeny - Mean AUC:\n")
print(mean(MF_pa$AUC, na.rm = TRUE))
cat("\nPA model without phylogeny - Mean Tjur R2:\n")
print(mean(MF_pa$TjurR2, na.rm = TRUE))

saveRDS(m_pa, file.path(results_dir, "HMSC_PA_no_phylogeny.rds"))
saveRDS(mpost_pa, file.path(results_dir, "mpost_PA_no_phylogeny.rds"))
saveRDS(preds_pa, file.path(results_dir, "preds_PA_no_phylogeny.rds"))
saveRDS(MF_pa, file.path(results_dir, "MF_PA_no_phylogeny_insample.rds"))

# 4-fold CV
set.seed(123)
partition_pa <- createPartition(m_pa, nfolds = 4)
preds_pa_cv <- computePredictedValues(
  m_pa, partition = partition_pa, nParallel = min(mcmc_config$n_parallel, 5)
)
MF_pa_cv <- evaluateModelFit(hM = m_pa, predY = preds_pa_cv)

cat("\nPA model CROSS-VALIDATED - Mean AUC:\n")
print(mean(MF_pa_cv$AUC, na.rm = TRUE))
cat("\nPA model CROSS-VALIDATED - Mean Tjur R2:\n")
print(mean(MF_pa_cv$TjurR2, na.rm = TRUE))
saveRDS(MF_pa_cv, file.path(results_dir, "MF_PA_no_phylogeny_cv.rds"))

# Variance partitioning
group_pa <- rep(NA_integer_, ncol(m_pa$X))
group_pa[colnames(m_pa$X) == "(Intercept)"] <- 1
group_pa[grepl("^HostType", colnames(m_pa$X)) & !grepl(":", colnames(m_pa$X))] <- 2
group_pa[grepl("^Season", colnames(m_pa$X)) & !grepl(":", colnames(m_pa$X))]   <- 3
group_pa[grepl("HostType.*:Season|Season.*:HostType", colnames(m_pa$X))]      <- 4
group_pa[grepl("^Site", colnames(m_pa$X))]     <- 5
group_pa[grepl("^logDepth_z$", colnames(m_pa$X))] <- 6
groupnames_pa <- c("(Intercept)", "HostType", "Season", "HostType:Season", "Site", "Depth")

if (any(is.na(group_pa))) {
  stop("Check the column names of the PA (no phylogeny) design matrix.")
}

VP_pa <- computeVariancePartitioning(m_pa, group = group_pa, groupnames = groupnames_pa)
saveRDS(VP_pa, file.path(results_dir, "VP_PA_no_phylogeny.rds"))

# ---------------------------------------------------------
# BLOCK 8. Conditional abundance, no phylogeny
# ---------------------------------------------------------
Y_abund <- log1p(Y_counts)
Y_abund[Y_counts == 0] <- NA

keep_taxa_abund <- colSums(!is.na(Y_abund)) >= 5
Y_abund <- Y_abund[, keep_taxa_abund, drop = FALSE]

studyDesign_abund <- data.frame(
  Sample = factor(rownames(Y_abund)),
  row.names = rownames(Y_abund)
)

rL_sample_abund <- HmscRandomLevel(units = levels(studyDesign_abund$Sample))
rL_sample_abund <- setPriors(rL_sample_abund, nfMax = 5)

m_abund <- Hmsc(
  Y = Y_abund,
  XData = XData[rownames(Y_abund), , drop = FALSE],
  XFormula = ~ HostType * Season + Site + logDepth_z,
  studyDesign = studyDesign_abund,
  ranLevels = list(Sample = rL_sample_abund),
  distr = "normal"
)

m_abund <- sampleMcmc(
  m_abund,
  samples   = mcmc_config$abund$samples,
  thin      = mcmc_config$abund$thin,
  transient = mcmc_config$abund$transient,
  nChains   = mcmc_config$abund$nChains,
  nParallel = mcmc_config$n_parallel,
  verbose   = 500
)

mpost_abund <- convertToCodaObject(m_abund)

cat("\nAbundance model without phylogeny - PSRF Beta:\n")
print(gelman.diag(mpost_abund$Beta, multivariate = TRUE)$mpsrf)
cat("\nAbundance model without phylogeny - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_abund$Beta)))

preds_abund <- computePredictedValues(m_abund, nParallel = mcmc_config$n_parallel)
MF_abund <- evaluateModelFit(m_abund, preds_abund)

cat("\nAbundance model without phylogeny - R2 summary:\n")
print(summary(MF_abund$R2))

saveRDS(m_abund, file.path(results_dir, "HMSC_abundance_no_phylogeny.rds"))
saveRDS(mpost_abund, file.path(results_dir, "mpost_abundance_no_phylogeny.rds"))
saveRDS(preds_abund, file.path(results_dir, "preds_abundance_no_phylogeny.rds"))
saveRDS(MF_abund, file.path(results_dir, "MF_abundance_no_phylogeny_insample.rds"))

# 4-fold CV
set.seed(123)
partition_abund <- createPartition(m_abund, nfolds = 4)
preds_abund_cv <- computePredictedValues(
  m_abund, partition = partition_abund, nParallel = min(mcmc_config$n_parallel, 5)
)
MF_abund_cv <- evaluateModelFit(hM = m_abund, predY = preds_abund_cv)

cat("\nAbundance model CROSS-VALIDATED - R2 summary:\n")
print(summary(MF_abund_cv$R2))
saveRDS(MF_abund_cv, file.path(results_dir, "MF_abundance_no_phylogeny_cv.rds"))

# Variance partitioning
group_abund <- rep(NA_integer_, ncol(m_abund$X))
group_abund[colnames(m_abund$X) == "(Intercept)"] <- 1
group_abund[grepl("^HostType", colnames(m_abund$X)) & !grepl(":", colnames(m_abund$X))] <- 2
group_abund[grepl("^Season", colnames(m_abund$X)) & !grepl(":", colnames(m_abund$X))]   <- 3
group_abund[grepl("HostType.*:Season|Season.*:HostType", colnames(m_abund$X))]          <- 4
group_abund[grepl("^Site", colnames(m_abund$X))]     <- 5
group_abund[grepl("^logDepth_z$", colnames(m_abund$X))] <- 6
groupnames_abund <- c("(Intercept)", "HostType", "Season", "HostType:Season", "Site", "Depth")

if (any(is.na(group_abund))) {
  stop("Check the column names of the abundance (no phylogeny) design matrix.")
}

VP_abund <- computeVariancePartitioning(m_abund, group = group_abund, groupnames = groupnames_abund)
saveRDS(VP_abund, file.path(results_dir, "VP_abundance_no_phylogeny.rds"))

cat("\nBlock 6-9 (HMSC models without phylogeny) complete.\n")
