# =========================================================
# HMSC MODELS WITH PHYLOGENY
# Corrected version of Blocks 12 and 14 from Microbiome_HMSC.Rmd
#
# CHANGE LOG vs. the original .Rmd:
#   - nChains was 4 in the .Rmd draft; the actual fitted objects
#     used nChains = 8. Corrected here.
#   - Abundance model samples/thin/transient were 12000/50/4000
#     in the .Rmd draft; the actual fitted object used
#     10000/20/6000. Corrected here.
#   - All values now come from R/00_mcmc_config.R so this script
#     and any future refit stay in sync automatically.
# =========================================================

library(Hmsc)
library(coda)

source("R/00_mcmc_config.R")

results_dir <- "results"

Y_counts   <- readRDS(file.path(results_dir, "Y_counts.rds"))
XData      <- readRDS(file.path(results_dir, "XData.rds"))
tree_genus <- readRDS(file.path(results_dir, "tree_genus.rds"))

# Restrict to genera present in the collapsed tree
common_genus <- intersect(tree_genus$tip.label, colnames(Y_counts))
tree_genus <- ape::keep.tip(tree_genus, common_genus)
Y_counts_phy <- Y_counts[, tree_genus$tip.label, drop = FALSE]
stopifnot(identical(colnames(Y_counts_phy), tree_genus$tip.label))

Y_pa_phy <- 1 * (Y_counts_phy > 0)

Y_abund_phy <- log1p(Y_counts_phy)
Y_abund_phy[Y_counts_phy == 0] <- NA
keep_taxa_abund_phy <- colSums(!is.na(Y_abund_phy)) >= 5
Y_abund_phy <- Y_abund_phy[, keep_taxa_abund_phy, drop = FALSE]

tree_genus_abund <- ape::keep.tip(tree_genus, colnames(Y_abund_phy))
Y_abund_phy <- Y_abund_phy[, tree_genus_abund$tip.label, drop = FALSE]

# ---------------------------------------------------------
# BLOCK 12. PRESENCE/ABSENCE MODEL WITH PHYLOGENY
# ---------------------------------------------------------

m_pa_phy <- Hmsc(
  Y = Y_pa_phy,
  XData = XData[rownames(Y_pa_phy), , drop = FALSE],
  XFormula = ~ HostType * Season + Site + logDepth_z,
  phyloTree = tree_genus,
  distr = "probit"
)

m_pa_phy <- sampleMcmc(
  m_pa_phy,
  samples   = mcmc_config$pa$samples,
  thin      = mcmc_config$pa$thin,
  transient = mcmc_config$pa$transient,
  nChains   = mcmc_config$pa$nChains,
  nParallel = mcmc_config$n_parallel,
  verbose   = 500
)

mpost_pa_phy <- convertToCodaObject(m_pa_phy)

cat("\nPA model with phylogeny - PSRF Beta:\n")
print(gelman.diag(mpost_pa_phy$Beta, multivariate = TRUE)$mpsrf)

cat("\nPA model with phylogeny - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_pa_phy$Beta)))

if (!is.null(mpost_pa_phy$Rho)) {
  cat("\nPA model with phylogeny - PSRF Rho:\n")
  print(gelman.diag(mpost_pa_phy$Rho, autoburnin = FALSE))

  cat("\nPA model with phylogeny - ESS Rho summary:\n")
  print(summary(effectiveSize(mpost_pa_phy$Rho)))
}

preds_pa_phy <- computePredictedValues(m_pa_phy, nParallel = mcmc_config$n_parallel)
MF_pa_phy <- evaluateModelFit(m_pa_phy, preds_pa_phy)

cat("\nPA model with phylogeny - Mean AUC:\n")
print(mean(MF_pa_phy$AUC, na.rm = TRUE))

cat("\nPA model with phylogeny - Mean Tjur R2:\n")
print(mean(MF_pa_phy$TjurR2, na.rm = TRUE))

saveRDS(m_pa_phy, file.path(results_dir, "HMSC_PA_with_phylogeny_logDepth_siteFixed.rds"))
saveRDS(mpost_pa_phy, file.path(results_dir, "mpost_PA_with_phylogeny.rds"))
saveRDS(preds_pa_phy, file.path(results_dir, "preds_PA_with_phylogeny.rds"))
saveRDS(MF_pa_phy, file.path(results_dir, "MF_PA_with_phylogeny_insample.rds"))

# Cross-validated version (4-fold)
set.seed(123)
partition_pa_phy <- createPartition(m_pa_phy, nfolds = 4)

preds_pa_phy_cv <- computePredictedValues(
  m_pa_phy,
  partition = partition_pa_phy,
  nParallel = min(mcmc_config$n_parallel, 5)
)

MF_pa_phy_cv <- evaluateModelFit(hM = m_pa_phy, predY = preds_pa_phy_cv)

cat("\nPA model with phylogeny CROSS-VALIDATED - Mean AUC:\n")
print(mean(MF_pa_phy_cv$AUC, na.rm = TRUE))

cat("\nPA model with phylogeny CROSS-VALIDATED - Mean Tjur R2:\n")
print(mean(MF_pa_phy_cv$TjurR2, na.rm = TRUE))

saveRDS(MF_pa_phy_cv, file.path(results_dir, "MF_PA_with_phylogeny_cv.rds"))

# Variance partitioning
group_pa_phy <- rep(NA_integer_, ncol(m_pa_phy$X))
group_pa_phy[colnames(m_pa_phy$X) == "(Intercept)"] <- 1
group_pa_phy[grepl("^HostType", colnames(m_pa_phy$X)) &
               !grepl(":", colnames(m_pa_phy$X))] <- 2
group_pa_phy[grepl("^Season", colnames(m_pa_phy$X)) &
               !grepl(":", colnames(m_pa_phy$X))] <- 3
group_pa_phy[grepl("HostType.*:Season|Season.*:HostType", colnames(m_pa_phy$X))] <- 4
group_pa_phy[grepl("^Site", colnames(m_pa_phy$X))] <- 5
group_pa_phy[grepl("^logDepth_z$", colnames(m_pa_phy$X))] <- 6

groupnames_pa_phy <- c(
  "(Intercept)", "HostType", "Season", "HostType:Season", "Site", "Depth"
)

if (any(is.na(group_pa_phy))) {
  cat("Columns without group assignment in m_pa_phy:\n")
  print(colnames(m_pa_phy$X)[is.na(group_pa_phy)])
  stop("Check the names of the columns in m_pa_phy$X.")
}

VP_PA_with_phylogeny <- computeVariancePartitioning(
  m_pa_phy,
  group = group_pa_phy,
  groupnames = groupnames_pa_phy
)

saveRDS(VP_PA_with_phylogeny, file.path(results_dir, "VP_PA_with_phylogeny.rds"))


# ---------------------------------------------------------
# BLOCK 14. CONDITIONAL ABUNDANCE MODEL WITH PHYLOGENY
# ---------------------------------------------------------

HMSC_hurdle_abundance_with_phylogeny_parallel <- Hmsc(
  Y = Y_abund_phy,
  XData = XData[rownames(Y_abund_phy), , drop = FALSE],
  XFormula = ~ HostType * Season + Site + logDepth_z,
  phyloTree = tree_genus_abund,
  distr = "normal"
)

HMSC_hurdle_abundance_with_phylogeny_parallel <- sampleMcmc(
  HMSC_hurdle_abundance_with_phylogeny_parallel,
  samples   = mcmc_config$abund$samples,
  thin      = mcmc_config$abund$thin,
  transient = mcmc_config$abund$transient,
  nChains   = mcmc_config$abund$nChains,
  nParallel = mcmc_config$n_parallel,
  verbose   = 500
)

mpost_abundance_with_phylogeny <- convertToCodaObject(HMSC_hurdle_abundance_with_phylogeny_parallel)

cat("\nAbundance model with phylogeny - PSRF Beta:\n")
print(gelman.diag(mpost_abundance_with_phylogeny$Beta, multivariate = TRUE)$mpsrf)

cat("\nAbundance model with phylogeny - ESS Beta summary:\n")
print(summary(effectiveSize(mpost_abundance_with_phylogeny$Beta)))

if (!is.null(mpost_abundance_with_phylogeny$Rho)) {
  cat("\nAbundance model with phylogeny - PSRF Rho:\n")
  print(gelman.diag(mpost_abundance_with_phylogeny$Rho, autoburnin = FALSE))

  cat("\nAbundance model with phylogeny - ESS Rho summary:\n")
  print(summary(effectiveSize(mpost_abundance_with_phylogeny$Rho)))
}

preds_abundance_with_phylogeny <- computePredictedValues(
  HMSC_hurdle_abundance_with_phylogeny_parallel,
  nParallel = mcmc_config$n_parallel
)

MF_abundance_with_phylogeny <- evaluateModelFit(
  HMSC_hurdle_abundance_with_phylogeny_parallel,
  preds_abundance_with_phylogeny
)

cat("\nAbundance model with phylogeny - R2 summary:\n")
print(summary(MF_abundance_with_phylogeny$R2))

cat("\nAbundance model with phylogeny - RMSE summary:\n")
print(summary(MF_abundance_with_phylogeny$RMSE))

cat("\nAbundance model with phylogeny - Mean R2:\n")
print(mean(MF_abundance_with_phylogeny$R2, na.rm = TRUE))

saveRDS(HMSC_hurdle_abundance_with_phylogeny_parallel, file.path(results_dir, "HMSC_hurdle_abundance_with_phylogeny_parallel.rds"))
saveRDS(mpost_abundance_with_phylogeny, file.path(results_dir, "mpost_abundance_with_phylogeny.rds"))
saveRDS(preds_abundance_with_phylogeny, file.path(results_dir, "preds_abundance_with_phylogeny.rds"))
saveRDS(MF_abundance_with_phylogeny, file.path(results_dir, "MF_abundance_with_phylogeny_insample.rds"))

# Cross-validated version (4-fold)
set.seed(123)
cv_partition_abundance <- createPartition(HMSC_hurdle_abundance_with_phylogeny_parallel, nfolds = 4)

preds_abundance_with_phylogeny_cv <- computePredictedValues(
  HMSC_hurdle_abundance_with_phylogeny_parallel,
  partition = cv_partition_abundance,
  nParallel = min(mcmc_config$n_parallel, 5)
)

MF_abundance_with_phylogeny_cv <- evaluateModelFit(
  hM = HMSC_hurdle_abundance_with_phylogeny_parallel,
  predY = preds_abundance_with_phylogeny_cv
)

cat("\nAbundance model with phylogeny CROSS-VALIDATED - R2 summary:\n")
print(summary(MF_abundance_with_phylogeny_cv$R2))

cat("\nAbundance model with phylogeny CROSS-VALIDATED - RMSE summary:\n")
print(summary(MF_abundance_with_phylogeny_cv$RMSE))

cat("\nAbundance model with phylogeny CROSS-VALIDATED - Mean R2:\n")
print(mean(MF_abundance_with_phylogeny_cv$R2, na.rm = TRUE))

saveRDS(MF_abundance_with_phylogeny_cv, file.path(results_dir, "MF_abundance_with_phylogeny_cv.rds"))

# Variance partitioning
group_abund_phy <- rep(NA_integer_, ncol(HMSC_hurdle_abundance_with_phylogeny_parallel$X))
group_abund_phy[colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X) == "(Intercept)"] <- 1
group_abund_phy[grepl("^HostType", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X)) &
                   !grepl(":", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X))] <- 2
group_abund_phy[grepl("^Season", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X)) &
                   !grepl(":", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X))] <- 3
group_abund_phy[grepl("HostType.*:Season|Season.*:HostType",
                       colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X))] <- 4
group_abund_phy[grepl("^Site", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X))] <- 5
group_abund_phy[grepl("^logDepth_z$", colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X))] <- 6

groupnames_abund_phy <- c(
  "(Intercept)", "HostType", "Season", "HostType:Season", "Site", "Depth"
)

if (any(is.na(group_abund_phy))) {
  cat("Columns without group assignment in HMSC_hurdle_abundance_with_phylogeny_parallel:\n")
  print(colnames(HMSC_hurdle_abundance_with_phylogeny_parallel$X)[is.na(group_abund_phy)])
  stop("Check the column names of the abundance model design matrix.")
}

VP_abundance_with_phylogeny <- computeVariancePartitioning(
  HMSC_hurdle_abundance_with_phylogeny_parallel,
  group = group_abund_phy,
  groupnames = groupnames_abund_phy
)

saveRDS(VP_abundance_with_phylogeny, file.path(results_dir, "VP_abundance_with_phylogeny.rds"))

# ---------------------------------------------------------
# Save XData_final: the exact covariate table used to fit both
# phylogeny-informed models. Figure scripts (06_figures/*) read
# this object by name, so it must be saved here rather than
# relying on a variable that only exists in an interactive
# session.
# ---------------------------------------------------------
XData_final <- XData[union(rownames(Y_pa_phy), rownames(Y_abund_phy)), , drop = FALSE]
saveRDS(XData_final, file.path(results_dir, "XData_final.rds"))

cat("\nBlock 12 (PA + phylogeny) and Block 14 (Abundance + phylogeny) complete.\n")
