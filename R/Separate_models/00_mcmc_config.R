# =========================================================
# MCMC CONFIGURATION — SINGLE SOURCE OF TRUTH
# =========================================================
# These values are confirmed directly from the fitted model
# objects (not from earlier drafts of the analysis code):
#
#   HMSC_PA_with_phylogeny_logDepth_siteFixed$samples    -> 12000
#   HMSC_PA_with_phylogeny_logDepth_siteFixed$thin        -> 50
#   HMSC_PA_with_phylogeny_logDepth_siteFixed$transient   -> 4000
#   length(HMSC_PA_with_phylogeny_logDepth_siteFixed$postList) -> 8
#
#   HMSC_hurdle_abundance_with_phylogeny_parallel$samples    -> 10000
#   HMSC_hurdle_abundance_with_phylogeny_parallel$thin        -> 20
#   HMSC_hurdle_abundance_with_phylogeny_parallel$transient   -> 6000
#   length(HMSC_hurdle_abundance_with_phylogeny_parallel$postList) -> 8
#
# Every script that fits, refits, or null-models these HMSC
# objects (with or without phylogeny) should source this file
# and read values from `mcmc_config` rather than hardcoding
# numbers. This guarantees that the manuscript's Methods /
# Table SXX and the actual code can never drift apart again.
# =========================================================

mcmc_config <- list(

  # Presence/absence models (probit)
  pa = list(
    samples   = 12000,
    thin      = 50,
    transient = 4000,
    nChains   = 8
  ),

  # Conditional abundance models (normal, hurdle)
  abund = list(
    samples   = 10000,
    thin      = 20,
    transient = 6000,
    nChains   = 8
  ),

  # Parallel workers for sampleMcmc() / computePredictedValues().
  # Adjust to the number of cores actually available on the
  # machine/cluster running the script; this does NOT change
  # the statistical results, only wall-clock time.
  n_parallel = 8
)

# Convenience helper: total post-burn-in draws per chain and
# across all chains, useful when reporting ESS relative to the
# maximum possible number of samples.
mcmc_total_draws <- function(cfg) {
  cfg$samples * cfg$nChains
}

cat("MCMC config loaded.\n")
cat("  PA model:       samples =", mcmc_config$pa$samples,
    "| thin =", mcmc_config$pa$thin,
    "| transient =", mcmc_config$pa$transient,
    "| nChains =", mcmc_config$pa$nChains,
    "| total draws =", mcmc_total_draws(mcmc_config$pa), "\n")
cat("  Abundance model: samples =", mcmc_config$abund$samples,
    "| thin =", mcmc_config$abund$thin,
    "| transient =", mcmc_config$abund$transient,
    "| nChains =", mcmc_config$abund$nChains,
    "| total draws =", mcmc_total_draws(mcmc_config$abund), "\n")
