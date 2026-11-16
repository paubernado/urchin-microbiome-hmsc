#!/usr/bin/env Rscript
# =============================================================================
# HMSC PIPELINE - STEP 1: MODEL DEFINITION + FITTING (+ convergence, in-sample)
# =============================================================================
# One script for the three models, chosen with --model:
#   pa     presence/absence (probit)
#   abund  abundance conditional on presence (normal on log1p reads, NA if absent)
#   joint  both in a single model (each genus twice: <genus>_PA and <genus>_ABUC),
#          as in Abrego et al. (2025)
#
# COMMON TO THE THREE MODELS
#   - data from STEP 0 (same samples, genera, covariates and rooted tree)
#   - presence = >= 5 reads (the same threshold defines PA and ABUC)
#   - fixed effects: HostType * Season + logDepth_z (water = reference level)
#   - random effects: sample (latent factors, number chosen by Hmsc) and site (1 factor)
#   - phylogeny from the rooted genus tree
#
# JOINT MODEL ONLY
#   - response-type "trait" (PA / ABUC): community responses (Gamma) estimated
#     separately for occurrence and abundance
#   - phylogeny as a block-diagonal correlation matrix C (Brownian motion,
#     ape::vcv(corr = TRUE), what Hmsc does internally with a tree), because
#     each genus is duplicated: phylogenetic correlation within PA and within
#     ABUC columns, none between them; a single Rho
#   - ABUC columns standardised inside the model (YScale = TRUE)
#
# BLOCKS
#   1  Load STEP 0 data
#   2  Response matrix
#   3  Phylogeny
#   4  Random levels
#   5  Define and fit the model (MCMC); save it before anything else
#   6  Convergence diagnostics (PSRF, ESS)
#   7  In-sample predictions and model fit
#
# MCMC DEFAULTS (change with --samples, --thin, --transient, --n-chains)
#   joint     : transient 50,000 + 250 samples x thin 400 = 150,000 iter/chain
#               (Abrego et al. 2025 used transient 12,500 + 250 x thin 100;
#               thinning increased for this larger model) -> 2,000 posterior samples
#   pa, abund : transient 6,000 + 10,000 samples x thin 20 = 206,000 iter/chain
#               (settings of the separate models already fitted)
#   8 chains, one core per chain (run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1)
#
# USAGE (hmsc_env active, inside screen):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 \
#     Rscript HMSC_STEP1_fit.R --model joint --datadir results_data --outdir results_joint
#   TEST RUN (minutes; 2 chains x 100 samples, output in test_<model>/):
#   Rscript HMSC_STEP1_fit.R --model joint --datadir results_data --test
#   then, to test STEP 2 on the same folder:
#   Rscript HMSC_STEP2_cv_results.R --model joint --outdir test_joint --task results --n-parallel 2
# =============================================================================

options(stringsAsFactors = FALSE, width = 200)

# -----------------------------
# User-tunable parameters
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg_value <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] + 1 > length(args)) stop(sprintf("Missing value for argument %s", flag))
  args[idx[1] + 1]
}
get_int_arg <- function(flag, default) as.integer(get_arg_value(flag, as.character(default)))

model_type <- get_arg_value("--model", "joint")
stopifnot(model_type %in% c("pa", "abund", "joint"))
datadir <- get_arg_value("--datadir", "results_data")

# TEST RUN: "--test" runs the whole model with very few iterations (minutes),
# only to check that the code works from start to end. Results are NOT valid.
# Outputs go to "test_<model>" so they never mix with the real results.
# STEP 2 run on that folder also refits with these short settings.
test_run <- "--test" %in% args

outdir  <- get_arg_value("--outdir", paste0(if (test_run) "test_" else "results_", model_type))

mcmc_default <- if (test_run) list(samples = 100, thin = 1, transient = 50, chains = 2) else
                if (model_type == "joint") list(samples = 250, thin = 400, transient = 50000, chains = 8) else
                                           list(samples = 10000, thin = 20, transient = 6000, chains = 8)
samples    <- get_int_arg("--samples",   mcmc_default$samples)
thin       <- get_int_arg("--thin",      mcmc_default$thin)
transient  <- get_int_arg("--transient", mcmc_default$transient)
n_chains   <- get_int_arg("--n-chains",  mcmc_default$chains)
n_parallel <- get_int_arg("--n-parallel", n_chains)
nf_min_sample <- get_int_arg("--nf-min-sample", 2)   # min latent factors at sample level (Hmsc default)
min_abuc_obs  <- get_int_arg("--min-abuc-obs", 5)    # min presences for an ABUC column
pred_draws    <- get_int_arg("--pred-draws", 1000)   # posterior draws used for predictions

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== HMSC pipeline - STEP 1:", model_type, "model ===\n")
if (test_run) cat("*** TEST RUN: very few iterations, results NOT valid (only checks the code) ***\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Data:", datadir, "| output:", outdir, "\n")
cat("MCMC: chains =", n_chains, "| samples =", samples, "| thin =", thin, "| transient =", transient,
    "| iterations/chain =", transient + samples * thin, "| nfMin(sample) =", nf_min_sample, "\n")

required_pkgs <- c("Hmsc", "ape", "coda")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({ library(Hmsc); library(ape); library(coda) })
cat("Hmsc version:", as.character(packageVersion("Hmsc")), "\n")

save_rds <- function(obj, filename) saveRDS(obj, file.path(outdir, filename))
model_file <- file.path(outdir, paste0("HMSC_", model_type, "_phylogeny_siteRandom.rds"))

# =============================================================================
# Bloc 1: Load STEP 0 data
# =============================================================================
Y_counts   <- readRDS(file.path(datadir, "Y_counts.rds"))
XData      <- readRDS(file.path(datadir, "XData.rds"))
tree_genus <- readRDS(file.path(datadir, "tree_genus_rooted.rds"))
data_set   <- readRDS(file.path(datadir, "data_settings.rds"))
min_reads  <- data_set$min_reads
stopifnot(all(rownames(Y_counts) == rownames(XData)),
          all(colnames(Y_counts) == tree_genus$tip.label), is.rooted(tree_genus))
cat("Samples:", nrow(Y_counts), "| genera:", ncol(Y_counts), "| tree rooted by:", data_set$root_method, "\n")

XD <- XData[, c("HostType", "Season", "Site", "logDepth_z")]

# =============================================================================
# Bloc 2: Response matrix
# =============================================================================
# PA   : 1 if >= 5 reads, 0 otherwise. Genera present in ALL samples have no
#        variation and are dropped from PA (they keep their ABUC column).
# ABUC : log(1 + reads) if >= 5 reads, NA otherwise (abundance conditional on
#        presence); genera need >= min_abuc_obs presences.
cat("\n--- Bloc 2: Response matrix ---\n")
present <- Y_counts >= min_reads
Y_pa    <- 1 * present
Y_pa    <- Y_pa[, apply(Y_pa, 2, var) > 0, drop = FALSE]
Y_abuc  <- log1p(Y_counts); Y_abuc[!present] <- NA
Y_abuc  <- Y_abuc[, colSums(!is.na(Y_abuc)) >= min_abuc_obs, drop = FALSE]

if (model_type == "pa") {
  Y <- Y_pa;   type <- factor(rep("PA",   ncol(Y)), levels = c("PA", "ABUC")); genus <- colnames(Y)
} else if (model_type == "abund") {
  Y <- Y_abuc; type <- factor(rep("ABUC", ncol(Y)), levels = c("PA", "ABUC")); genus <- colnames(Y)
} else {
  colnames(Y_pa)   <- paste0(colnames(Y_pa),   "_PA")
  colnames(Y_abuc) <- paste0(colnames(Y_abuc), "_ABUC")
  Y <- cbind(Y_pa, Y_abuc)
  type  <- factor(sub(".*_(PA|ABUC)$", "\\1", colnames(Y)), levels = c("PA", "ABUC"))
  genus <- sub("_(PA|ABUC)$", "", colnames(Y))
}
distr <- ifelse(type == "PA", "probit", "normal")
cat("Columns:", ncol(Y), "| PA:", sum(type == "PA"), "| ABUC:", sum(type == "ABUC"), "\n")

# =============================================================================
# Bloc 3: Phylogeny
# =============================================================================
if (model_type == "joint") {
  # Block-diagonal phylogenetic correlation matrix (each genus appears twice)
  C_genus <- vcv(tree_genus, corr = TRUE)
  C <- C_genus[genus, genus]
  C[type == "PA", type == "ABUC"] <- 0
  C[type == "ABUC", type == "PA"] <- 0
  dimnames(C) <- list(colnames(Y), colnames(Y))
  stopifnot(isSymmetric(C))
  tree_model <- NULL
} else {
  # Separate models: the rooted tree is passed directly (Hmsc builds C itself)
  tree_model <- keep.tip(tree_genus, colnames(Y))
  Y <- Y[, tree_model$tip.label, drop = FALSE]
  genus <- colnames(Y)
  C <- NULL
}

# =============================================================================
# Bloc 4: Random levels
# =============================================================================
# sample: one unit per sample; latent factors give the residual associations
#         Omega. Their number is chosen adaptively by Hmsc during the transient
#         (minimum nf_min_sample = 2, the Hmsc default; the separate abundance
#         model ended up with 6 factors in every chain).
# site:   the 4 sites; only 1 latent factor (4 units are too few for more).
studyDesign <- data.frame(sample = factor(rownames(XD)), site = factor(XD$Site))
# nfMax is given explicitly (Inf = Hmsc default, adaptive): some Hmsc versions
# (e.g. 3.3.7) fail when nfMin is set without nfMax
rL_sample <- setPriors(HmscRandomLevel(units = levels(studyDesign$sample)),
                       nfMax = Inf, nfMin = nf_min_sample)
rL_site   <- setPriors(HmscRandomLevel(units = levels(studyDesign$site)), nfMax = 1, nfMin = 1)

# =============================================================================
# Bloc 5: Define and fit the model
# =============================================================================
cat("\n--- Bloc 5: Model definition and MCMC ---\n")
hmsc_args <- list(Y = Y, XData = XD, XFormula = ~ HostType * Season + logDepth_z,
                  studyDesign = studyDesign, ranLevels = list(sample = rL_sample, site = rL_site),
                  distr = if (model_type == "joint") distr else distr[1])
if (model_type == "joint") {
  hmsc_args$TrData    <- data.frame(type = type, row.names = colnames(Y))  # PA / ABUC
  hmsc_args$TrFormula <- ~ type
  hmsc_args$C         <- C
  hmsc_args$YScale    <- TRUE   # standardise ABUC columns; predictions back on log1p scale
} else {
  hmsc_args$phyloTree <- tree_model
}
m <- do.call(Hmsc, hmsc_args)

save_rds(Y, paste0("Y_", model_type, ".rds"))
save_rds(XD, paste0("XData_", model_type, ".rds"))
if (!is.null(tree_model)) save_rds(tree_model, paste0("tree_", model_type, ".rds"))

t0 <- Sys.time()
m <- sampleMcmc(m, samples = samples, thin = thin, transient = transient,
                nChains = n_chains, nParallel = n_parallel, useSocket = FALSE,
                verbose = max(1, round((transient + samples * thin) / 100)))
cat("MCMC run time:", format(Sys.time() - t0), "\n")

# Save the model with all chains BEFORE any other computation
saveRDS(m, model_file, compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, paste0("sessionInfo_step1_", model_type, ".txt")))
cat("Model saved:", model_file, "\n")

# =============================================================================
# Bloc 6: Convergence diagnostics
# =============================================================================
# PSRF (Gelman-Rubin) ~1 = chains agree; criterion: median ~1.00-1.05 and few
# parameters > 1.1. ESS = effectively independent samples.
# Only identifiable parameters are converted (Eta, Lambda, Omega would take
# many GB); Omega is checked on random pairs of columns, computed from Lambda.
cat("\n--- Bloc 6: Convergence ---\n")
cc_args <- list(hM = m, Beta = TRUE, Gamma = TRUE, V = TRUE, Sigma = TRUE,
                Eta = FALSE, Lambda = FALSE, Omega = FALSE, Psi = FALSE, Delta = FALSE)
f <- names(formals(convertToCodaObject))             # "Rho" (Hmsc <= 3.3) or "rho" (newer)
if ("Rho" %in% f) cc_args$Rho <- TRUE else cc_args$rho <- TRUE
if ("Alpha" %in% f) cc_args$Alpha <- FALSE else cc_args$alpha <- FALSE
mpost <- do.call(convertToCodaObject, cc_args)
rho_post <- if (!is.null(mpost$Rho)) mpost$Rho else mpost$rho
save_rds(mpost, paste0("mpost_", model_type, ".rds"))

conv <- function(x, label, max_par = 1000) {
  if (is.null(x)) return(NULL)
  set.seed(1)
  np <- ncol(as.matrix(x[[1]]))
  if (np > max_par) x <- x[, sort(sample(np, max_par)), drop = FALSE]
  x <- x[, apply(as.matrix(x[[1]]), 2, var) > 0, drop = FALSE]     # drop fixed parameters
  if (ncol(as.matrix(x[[1]])) == 0) return(NULL)
  psrf <- gelman.diag(x, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
  ess  <- effectiveSize(x)
  data.frame(Parameter = label, n_checked = length(psrf),
             PSRF_median = round(median(psrf), 3), PSRF_q95 = round(quantile(psrf, 0.95), 3),
             Pct_PSRF_gt_1.1 = round(100 * mean(psrf > 1.1), 1),
             ESS_median = round(median(ess)), ESS_min = round(min(ess)),
             ESS_max_possible = length(x) * nrow(as.matrix(x[[1]])))
}
omega_pairs <- function(r, i_set, j_set, n = 300) {
  set.seed(1)
  pr <- cbind(sample(i_set, n, TRUE), sample(j_set, n, TRUE))
  pr <- unique(pr[pr[, 1] != pr[, 2], , drop = FALSE])
  mcmc.list(lapply(m$postList, function(ch) mcmc(t(sapply(ch, function(s) {
    L <- s$Lambda[[r]]
    if (length(dim(L)) == 3) L <- matrix(L[, , 1], nrow = dim(L)[1])
    colSums(L[, pr[, 1], drop = FALSE] * L[, pr[, 2], drop = FALSE])
  })))))
}

beta_col <- rep(seq_len(m$ns), each = m$nc)            # column of each Beta parameter
pa_i <- which(type == "PA"); ab_i <- which(type == "ABUC")
conv_tab <- rbind(
  if (length(pa_i)) conv(mpost$Beta[, type[beta_col] == "PA",   drop = FALSE], "Beta (PA)"),
  if (length(ab_i)) conv(mpost$Beta[, type[beta_col] == "ABUC", drop = FALSE], "Beta (ABUC)"),
  conv(mpost$Gamma, "Gamma"), conv(mpost$V, "V"), conv(rho_post, "Rho"),
  if (length(ab_i)) conv(mpost$Sigma[, ab_i, drop = FALSE], "Sigma (ABUC)"))
for (r in seq_along(m$rLNames)) {
  lv <- m$rLNames[r]
  if (length(pa_i) > 1) conv_tab <- rbind(conv_tab, conv(omega_pairs(r, pa_i, pa_i), paste0("Omega ", lv, " (PA-PA)")))
  if (length(ab_i) > 1) conv_tab <- rbind(conv_tab, conv(omega_pairs(r, ab_i, ab_i), paste0("Omega ", lv, " (ABUC-ABUC)")))
  if (length(pa_i) && length(ab_i)) conv_tab <- rbind(conv_tab, conv(omega_pairs(r, pa_i, ab_i), paste0("Omega ", lv, " (PA-ABUC)")))
}
print(conv_tab, row.names = FALSE)
write.csv(conv_tab, file.path(outdir, paste0("convergence_", model_type, ".csv")), row.names = FALSE)

# Chains stuck in different solutions? Mean sigma per chain (ABUC columns)
if (length(ab_i)) {
  sig_chain <- sapply(mpost$Sigma, function(ch) colMeans(as.matrix(ch)[, ab_i, drop = FALSE]))
  ratio <- apply(sig_chain, 1, max) / apply(sig_chain, 1, min)
  cat("ABUC columns whose chains differ > 2x in mean sigma:", sum(ratio > 2), "\n")
  write.csv(data.frame(Column = colnames(Y)[ab_i], round(sig_chain, 3), max_min_ratio = round(ratio, 2)),
            file.path(outdir, paste0("sigma_by_chain_", model_type, ".csv")), row.names = FALSE)
}

# =============================================================================
# Bloc 7: In-sample predictions and model fit
# =============================================================================
# ~1,000 posterior draws are enough for predictions and keep files manageable
cat("\n--- Bloc 7: In-sample predictions and model fit ---\n")
pthin <- max(1, floor(samples * n_chains / pred_draws))
preds <- computePredictedValues(m, thin = pthin, nParallel = n_parallel)
saveRDS(preds, file.path(outdir, paste0("preds_", model_type, "_insample.rds")), compress = FALSE)
MF <- evaluateModelFit(m, preds)
save_rds(MF, paste0("MF_", model_type, "_insample.rds"))
if (length(pa_i)) cat("PA   - mean AUC:", round(mean(MF$AUC[pa_i], na.rm = TRUE), 3),
                      "| mean Tjur R2:", round(mean(MF$TjurR2[pa_i], na.rm = TRUE), 3), "\n")
if (length(ab_i)) cat("ABUC - mean R2:", round(mean(MF$R2[ab_i], na.rm = TRUE), 3),
                      "| mean RMSE:", round(mean(MF$RMSE[ab_i], na.rm = TRUE), 3), "\n")

if (test_run) cat("\n*** TEST RUN finished OK. Run again without --test for the real model. ***\n")
cat("\n=== STEP 1 (", model_type, ") finished ===\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
