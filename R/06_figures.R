# =========================================================
# MASTER FIGURE-GENERATION SCRIPT
# Loads every object the figure sub-scripts need directly from
# the .rds files saved by the earlier pipeline steps, then
# sources the four figure scripts in the paper's figure order:
#
#   1) R/figures/fig2_composition_diversity.R
#      -> Figure2_composition_and_diversity.{pdf,png}
#   2) R/figures/fig_variance_partitioning.R
#      -> FIGURE_VP_both_models.{pdf,png}
#   3) R/figures/fig_seasonal_deviations.R
#      -> FIGURE_seasonal_deviations_FINAL.{pdf,png}
#   4) R/figures/fig_host_contrasts.R
#      -> FIGURE_block4_ABUNDANCE.{pdf,png}, FIGURE_block4_PA.{pdf,png}
#
# Each sub-script is self-contained (own library() calls, own
# stopifnot() checks) so it can also be run/sourced on its own,
# as long as the objects below have been loaded first.
# =========================================================

# ---------------------------------------------------------
# 0) Paths — adjust if your saved .rds files live elsewhere
# ---------------------------------------------------------
results_dir <- "results"   # where 01-05 scripts saved their .rds output

# ---------------------------------------------------------
# 1) Load objects shared across figure scripts
# ---------------------------------------------------------

# From 01_import_decontam.R
phylo_raw <- readRDS(file.path(results_dir, "phylo_raw.rds"))

# From 04_hmsc_models_with_phylogeny.R
XData_final <- readRDS(file.path(results_dir, "XData_final.rds"))

VP_abundance_with_phylogeny <- readRDS(file.path(results_dir, "VP_abundance_with_phylogeny.rds"))
VP_PA_with_phylogeny        <- readRDS(file.path(results_dir, "VP_PA_with_phylogeny.rds"))

preds_abundance_with_phylogeny <- readRDS(file.path(results_dir, "preds_abundance_with_phylogeny.rds"))
preds_PA_with_phylogeny        <- readRDS(file.path(results_dir, "preds_PA_with_phylogeny.rds"))

mpost_abundance_with_phylogeny <- readRDS(file.path(results_dir, "mpost_abundance_with_phylogeny.rds"))
mpost_PA_with_phylogeny        <- readRDS(file.path(results_dir, "mpost_PA_with_phylogeny.rds"))

cat("All figure input objects loaded.\n")

# ---------------------------------------------------------
# 2) Source figure scripts in order
# ---------------------------------------------------------
source("R/figures/fig2_composition_diversity.R")
source("R/figures/fig_variance_partitioning.R")
source("R/figures/fig_seasonal_deviations.R")
source("R/figures/fig_host_contrasts.R")

cat("\nAll figures generated.\n")
