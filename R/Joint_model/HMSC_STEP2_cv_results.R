#!/usr/bin/env Rscript
# =============================================================================
# HMSC PIPELINE - STEP 2: CROSS-VALIDATION + RESULTS
# Run after HMSC_STEP1_fit.R, for the same --model (pa, abund or joint).
# =============================================================================
#
# TASKS (--task); each can be sent as a separate cluster job:
#   insample        in-sample predictions and fit (only needed for models not
#                   fitted with HMSC_STEP1_fit.R)
#   cv_random       random 4-fold cross-validation across samples, stratified
#                   by HostType x Season (every fold has every combination)
#   cv_site         leave-one-site-out cross-validation: each site is predicted
#                   by a model fitted to the other 3 (transferability to new
#                   sites, the reviewer's point)
#   cv_conditional  (joint model only) same folds as cv_random, but PA columns
#                   are predicted given the observed ABUC of the test samples and
#                   vice versa (adapted from Abrego et al. 2025). If it beats
#                   cv_random, occurrence and abundance share information.
#   results         model-fit table, variance partitioning, Beta, Gamma, Rho
#                   (joint: also post hoc Rho for PA and for ABUC),
#                   sample-level residual associations (Omega) and predictions
#                   for every HostType x Season combination
#
# NOTE: each CV refits the model once per fold with the MCMC settings of STEP 1,
# so a CV task takes ~4 times as long as the STEP 1 fit.
#
# USAGE (hmsc_env active):
#   Rscript HMSC_STEP2_cv_results.R --model joint --outdir results_joint --task cv_site
#   Rscript HMSC_STEP2_cv_results.R --model joint --outdir results_joint --task results
#   Test run (after HMSC_STEP1_fit.R --test): use --outdir test_joint
#   Model fitted with another script:
#   Rscript HMSC_STEP2_cv_results.R --model pa --outdir results_phy_PA_SD_ranlev \
#           --model-file results_phy_PA_SD_ranlev/HMSC_PA_with_phylogeny_logDepth_siteRandom.rds --task cv_site
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
model_type <- get_arg_value("--model", "joint")
stopifnot(model_type %in% c("pa", "abund", "joint"))
outdir     <- get_arg_value("--outdir", paste0("results_", model_type))
task       <- get_arg_value("--task", "results")
n_parallel <- as.integer(get_arg_value("--n-parallel", "8"))
cv_folds   <- as.integer(get_arg_value("--cv-folds", "4"))
pred_draws <- as.integer(get_arg_value("--pred-draws", "1000"))   # posterior draws for predictions
# --parallel-folds yes: fit all CV folds at the same time (folds x chains cores,
# e.g. 4 x 8 = 32) with Hmsc::pcomputePredictedValues; "no" fits the folds one
# after the other (8 cores) with computePredictedValues. Same result, ~4x faster.
par_folds  <- tolower(get_arg_value("--parallel-folds", "yes")) == "yes"
# --cv-samples N: posterior samples PER CHAIN kept in each CV fold model. The
# CV models keep the same transient and thin as the full model (same burn-in
# and same autocorrelation between draws) but stop after N samples per chain.
# 0 = same number of samples as the full model.
cv_samples <- as.integer(get_arg_value("--cv-samples", "0"))

if (task == "cv_conditional" && model_type != "joint") stop("cv_conditional is only for the joint model.")

suppressPackageStartupMessages({ library(Hmsc); library(ape); library(coda) })
cat("=== HMSC pipeline - STEP 2 | model:", model_type, "| task:", task, "===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
save_rds <- function(obj, filename) saveRDS(obj, file.path(outdir, filename))

# -----------------------------
# Load the fitted model from STEP 1
# -----------------------------
# By default the model saved by HMSC_STEP1_fit.R; --model-file allows using a
# model fitted earlier with another script (e.g. the separate PA / abundance
# models already fitted on the cluster). Outputs are written to --outdir.
model_file <- get_arg_value("--model-file",
                            file.path(outdir, paste0("HMSC_", model_type, "_phylogeny_siteRandom.rds")))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(model_file)) stop("Model not found: ", model_file, "\nRun STEP 1 first.")
m <- readRDS(model_file)
n_parallel <- min(n_parallel, length(m$postList))   # never more cores than chains (e.g. test runs)

# Column names may differ in models fitted with older scripts
# (e.g. "Species" instead of "HostType", "Site" instead of "site")
if (is.null(m$XData$HostType) && !is.null(m$XData$Species)) m$XData$HostType <- m$XData$Species
site_col <- names(m$studyDesign)[tolower(names(m$studyDesign)) == "site"][1]
if (is.na(site_col)) stop("No 'site' column in studyDesign: ", paste(names(m$studyDesign), collapse = ", "))

if (model_type == "joint") {
  type  <- factor(sub(".*_(PA|ABUC)$", "\\1", m$spNames), levels = c("PA", "ABUC"))
  genus <- sub("_(PA|ABUC)$", "", m$spNames)
} else {
  type  <- factor(rep(if (model_type == "pa") "PA" else "ABUC", m$ns), levels = c("PA", "ABUC"))
  genus <- m$spNames
}
pa <- type == "PA"; ab <- type == "ABUC"
cat("Columns:", m$ns, "| PA:", sum(pa), "| ABUC:", sum(ab), "| samples:", m$ny,
    "| chains:", length(m$postList), "| samples/chain:", m$samples,
    "| thin:", m$thin, "| transient:", m$transient, "\n")

# Posterior thinning used ONLY for predictions (~pred_draws draws in total);
# model estimates (Beta, VP, ...) use the full posterior.
pthin <- max(1, floor(m$samples * length(m$postList) / pred_draws))

fit_summary <- function(MF, label) {
  g <- function(v, sel) if (is.null(v) || !any(sel)) NA else round(mean(v[sel], na.rm = TRUE), 3)
  data.frame(Evaluation = label, PA_AUC = g(MF$AUC, pa), PA_TjurR2 = g(MF$TjurR2, pa),
             ABUC_R2 = g(MF$R2, ab), ABUC_RMSE = g(MF$RMSE, ab))
}

# The random partition is saved and reused, so cv_random and cv_conditional
# use exactly the same folds and can be compared directly.
#
# The folds are STRATIFIED by HostType x Season: the samples of every
# HostType x Season combination are spread evenly across folds, so that no
# fold lacks a level (e.g. a fold without seawater samples).
get_random_partition <- function() {
  f <- file.path(outdir, paste0("cv_partition_random_", model_type, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  set.seed(123)
  strata <- interaction(m$XData$HostType, m$XData$Season, drop = TRUE)
  p <- integer(m$ny); offset <- sample(0:(cv_folds - 1), 1)
  for (lev in levels(strata)) {
    idx <- which(strata == lev)
    labels <- ((seq_along(idx) + offset - 1) %% cv_folds) + 1   # balanced fold labels
    p[idx] <- sample(labels)                                      # shuffled within stratum
    offset <- offset + length(idx)
  }
  saveRDS(p, f); p
}

# Check that every fold (test set) and every training set contains all
# HostType x Season combinations
check_folds <- function(part, label) {
  strata <- interaction(m$XData$HostType, m$XData$Season, drop = TRUE)
  tab <- table(fold = part, stratum = strata)
  cat("\n", label, "- samples per fold and HostType x Season:\n", sep = ""); print(tab)
  train_missing <- sapply(sort(unique(part)), function(k) sum(table(strata[part != k]) == 0))
  if (any(train_missing > 0)) stop("Some training set lacks a HostType x Season combination.")
  if (any(tab == 0)) cat("NOTE: some test folds lack a HostType x Season combination.\n")
  else cat("OK: every fold contains all HostType x Season combinations.\n")
}

run_cv <- function(part, label, file_tag, ...) {
  extra <- list(...)
  mcv <- m; thin_cv <- pthin
  if (cv_samples > 0 && cv_samples < m$samples) {
    # computePredictedValues refits each fold with mcv$samples, mcv$thin and
    # mcv$transient; postList is shortened only so that the array sizes match.
    mcv$samples  <- cv_samples
    mcv$postList <- lapply(mcv$postList, function(ch) ch[seq_len(cv_samples)])
    thin_cv <- max(1, floor(cv_samples * length(m$postList) / pred_draws))
  }
  cat("CV fold models: samples/chain =", mcv$samples, "| thin =", mcv$thin,
      "| transient =", mcv$transient, "| iterations/chain =",
      mcv$transient + mcv$samples * mcv$thin, "\n")
  if (par_folds && is.null(extra$partition.sp) && exists("pcomputePredictedValues")) {
    n_threads <- length(unique(part)) * length(m$postList)       # folds x chains
    cat("Fitting all folds in parallel:", n_threads, "threads (folds x chains)\n")
    preds <- pcomputePredictedValues(mcv, partition = part, thin = thin_cv, nParallel = n_threads, ...)
  } else {
    cat("Fitting folds one after the other,", n_parallel, "cores\n")
    preds <- computePredictedValues(mcv, partition = part, thin = thin_cv, nParallel = n_parallel, ...)
  }
  cat("CV finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  saveRDS(preds, file.path(outdir, paste0("preds_", model_type, "_", file_tag, ".rds")), compress = FALSE)
  MF <- evaluateModelFit(m, preds)
  save_rds(MF, paste0("MF_", model_type, "_", file_tag, ".rds"))
  print(fit_summary(MF, label), row.names = FALSE)
}

# =============================================================================
# Cross-validation tasks
# =============================================================================
if (task == "insample") {
  # In-sample predictions (already done by HMSC_STEP1_fit.R; this task is for
  # models fitted with other scripts)
  preds <- computePredictedValues(m, thin = pthin, nParallel = n_parallel)
  saveRDS(preds, file.path(outdir, paste0("preds_", model_type, "_insample.rds")), compress = FALSE)
  MF <- evaluateModelFit(m, preds)
  save_rds(MF, paste0("MF_", model_type, "_insample.rds"))
  print(fit_summary(MF, "In-sample"), row.names = FALSE)
}

if (task == "cv_random") {
  part <- get_random_partition()
  check_folds(part, "Random CV")
  run_cv(part, paste0("Random ", cv_folds, "-fold CV"), "cv_random")
}

if (task == "cv_site") {
  part <- createPartition(m, nfolds = nlevels(factor(m$studyDesign[[site_col]])), column = site_col)
  print(table(fold = part, site = m$studyDesign[[site_col]]))
  check_folds(part, "Leave-one-site-out CV")
  run_cv(part, "Leave-one-site-out CV", "cv_site")
}

if (task == "cv_conditional") {
  # partition.sp: group 1 = PA columns, group 2 = ABUC columns; each group is
  # predicted conditional on the observed values of the other group.
  part <- get_random_partition()
  run_cv(part, "Conditional CV", "cv_conditional",
         partition.sp = as.integer(type), mcmcStep = 100)
}

# =============================================================================
# Results
# =============================================================================
if (task == "results") {

  # ---- Model fit table (whatever has been computed so far) ------------------
  mf <- c("In-sample" = "insample", "Random CV" = "cv_random",
          "Leave-one-site-out CV" = "cv_site", "Conditional CV" = "cv_conditional")
  fit_tab <- do.call(rbind, lapply(names(mf), function(n) {
    f <- file.path(outdir, paste0("MF_", model_type, "_", mf[n], ".rds"))
    if (file.exists(f)) fit_summary(readRDS(f), n)
  }))
  cat("\n--- Model fit ---\n"); print(fit_tab, row.names = FALSE)
  write.csv(fit_tab, file.path(outdir, paste0("model_fit_", model_type, ".csv")), row.names = FALSE)

  # ---- Variance partitioning ------------------------------------------------
  # Fixed-effect groups only. Site and sample are random levels: they are not
  # columns of X and Hmsc adds them automatically ("Random: sample/site").
  cn <- colnames(m$X)
  grp <- rep(NA_integer_, length(cn))
  grp[cn == "(Intercept)"] <- 1
  grp[grepl("^HostType", cn) & !grepl(":", cn)] <- 2
  grp[grepl("^Season", cn) & !grepl(":", cn)] <- 3
  grp[grepl("HostType.*:Season|Season.*:HostType", cn)] <- 4
  grp[cn == "logDepth_z"] <- 5
  if (any(is.na(grp))) { print(cn[is.na(grp)]); stop("Unassigned columns in m$X.") }
  VP <- computeVariancePartitioning(m, group = grp,
          groupnames = c("(Intercept)", "HostType", "Season", "HostType:Season", "Depth"))
  rn <- sub("^Random: ", "", rownames(VP$vals))
  rn[rn == "sample"] <- "Sample"; rn[rn == "site"] <- "Site"
  rownames(VP$vals) <- rn
  VP_mean <- data.frame(Component = rn)
  if (any(pa)) VP_mean$PA   <- round(rowMeans(VP$vals[, pa, drop = FALSE]), 3)
  if (any(ab)) VP_mean$ABUC <- round(rowMeans(VP$vals[, ab, drop = FALSE]), 3)
  cat("\n--- Mean variance explained ---\n"); print(VP_mean, row.names = FALSE)
  save_rds(VP, paste0("VP_", model_type, ".rds"))
  write.csv(VP_mean, file.path(outdir, paste0("VP_", model_type, "_mean.csv")), row.names = FALSE)
  write.csv(data.frame(Column = m$spNames, Genus = genus, Type = type, t(VP$vals)),
            file.path(outdir, paste0("VP_", model_type, "_by_column.csv")), row.names = FALSE)

  # ---- Beta: response of each column, with support -------------------------
  # support > 0.95 = positive effect; < 0.05 = negative effect
  pb <- getPostEstimate(m, parName = "Beta")
  write.csv(data.frame(Column = rep(m$spNames, each = m$nc), Genus = rep(genus, each = m$nc),
                       Type = rep(type, each = m$nc), Covariate = rep(m$covNames, times = m$ns),
                       Mean = as.vector(pb$mean), Support = as.vector(pb$support)),
            file.path(outdir, paste0("Beta_", model_type, ".csv")), row.names = FALSE)

  # ---- Gamma: community-level responses ------------------------------------
  # In the joint model, the "typeABUC" rows are the difference ABUC minus PA.
  pg <- getPostEstimate(m, parName = "Gamma")
  write.csv(data.frame(Covariate = rep(m$covNames, times = m$nt), Trait = rep(m$trNames, each = m$nc),
                       Mean = as.vector(pg$mean), Support = as.vector(pg$support)),
            file.path(outdir, paste0("Gamma_", model_type, ".csv")), row.names = FALSE)

  # ---- Rho: phylogenetic signal ----------------------------------------------
  rho_val <- unlist(lapply(m$postList, function(ch) sapply(ch, function(s) s$rho)))
  rho_tab <- data.frame(mean = mean(rho_val), q2.5 = unname(quantile(rho_val, 0.025)),
                        q97.5 = unname(quantile(rho_val, 0.975)))
  cat(sprintf("\nRho: mean %.3f, 95%% CI [%.3f, %.3f]\n", rho_tab$mean, rho_tab$q2.5, rho_tab$q97.5))
  write.csv(rho_tab, file.path(outdir, paste0("Rho_", model_type, ".csv")), row.names = FALSE)

  # ---- Joint model: post hoc phylogenetic signal of each component ---------
  # The joint model estimates a single Rho. As suggested by Dave, a separate
  # signal for PA and for ABUC can be obtained post hoc from Beta, Gamma, V and
  # the C matrix: for each posterior draw, the residual responses
  # E = Beta - Gamma %*% t(Tr) of each block (PA or ABUC columns) are evaluated
  # over the Hmsc grid of Rho values with the same likelihood Hmsc uses
  # (covariance among columns = Rho * C + (1 - Rho) * I), and a Rho is drawn for
  # each block. Because C is block-diagonal, this is the conditional posterior
  # of a block-specific Rho given the other parameters of the joint model.
  if (model_type == "joint") {
    blocks <- list(PA = which(pa), ABUC = which(ab))
    eig <- lapply(blocks, function(i) {
      e <- eigen(m$C[i, i], symmetric = TRUE); e$values <- pmax(e$values, 1e-10); e })
    grid <- m$rhopw                                  # column 1 = Rho values, column 2 = prior weights
    post_r <- poolMcmcChains(m$postList, thin = pthin)
    set.seed(1)
    rho_blocks <- t(sapply(post_r, function(sm) {
      E  <- sm$Beta - sm$Gamma %*% t(m$Tr)           # covariates x columns
      iV <- chol2inv(chol(sm$V))
      sapply(names(blocks), function(b) {
        ev <- eig[[b]]
        Fm <- E[, blocks[[b]], drop = FALSE] %*% ev$vectors
        q  <- colSums(Fm * (iV %*% Fm))
        ll <- sapply(grid[, 1], function(r) {
          w <- r * ev$values + 1 - r
          -0.5 * nrow(E) * sum(log(w)) - 0.5 * sum(q / w)
        }) + log(grid[, 2])
        sample(grid[, 1], 1, prob = exp(ll - max(ll)))
      })
    }))
    rho_comp <- data.frame(Component = c("Joint (single Rho)", "PA (post hoc)", "ABUC (post hoc)"),
                           mean  = c(mean(rho_val), colMeans(rho_blocks)),
                           q2.5  = c(quantile(rho_val, 0.025), apply(rho_blocks, 2, quantile, 0.025)),
                           q97.5 = c(quantile(rho_val, 0.975), apply(rho_blocks, 2, quantile, 0.975)))
    cat("\nPhylogenetic signal by component (post hoc):\n"); print(rho_comp, row.names = FALSE)
    write.csv(rho_comp, file.path(outdir, "Rho_by_component_joint.csv"), row.names = FALSE)
    cat("P(Rho_PA > Rho_ABUC):", round(mean(rho_blocks[, "PA"] > rho_blocks[, "ABUC"]), 3), "\n")
  }

  # ---- Residual associations at the sample level -------------------------------
  OmegaCor <- computeAssociations(m, thin = pthin)
  Om <- OmegaCor[[which(m$rLNames == "sample")]]
  sig <- (Om$support > 0.95) - (Om$support < 0.05); diag(sig) <- 0
  blk <- function(i, j) { s <- sig[i, j, drop = FALSE]; if (identical(i, j)) s <- s[upper.tri(s)]
                          c(pairs = length(s), positive = sum(s == 1), negative = sum(s == -1)) }
  om_tab <- rbind(if (any(pa)) c(Block = "PA-PA", blk(pa, pa)),
                  if (any(ab)) c(Block = "ABUC-ABUC", blk(ab, ab)),
                  if (any(pa) && any(ab)) c(Block = "PA-ABUC", blk(pa, ab)))
  cat("\nSample-level residual associations with 95% support:\n"); print(om_tab)
  save_rds(OmegaCor, paste0("Omega_", model_type, ".rds"))
  write.csv(om_tab, file.path(outdir, paste0("Omega_sample_summary_", model_type, ".csv")), row.names = FALSE)

  # ---- Predictions for every HostType x Season -------------------------------
  # Depth at its mean (logDepth_z = 0); new sample and new site, i.e. marginal
  # over sample- and site-level variation.
  #   pa    : P(present)
  #   abund : E[log1p reads | present]
  #   joint : both + hurdle expectation P(present) x E[log1p | present],
  #           computed within each posterior draw (genera present in all
  #           samples have no PA column: P = 1)
  post <- poolMcmcChains(m$postList, thin = pthin)
  rows <- list(); draws_list <- list()
  for (s in levels(m$XData$Season)) {
    G <- constructGradient(m, focalVariable = "HostType",
                           non.focalVariables = list(Season = list(3, s), logDepth_z = list(3, 0)))
    p <- simplify2array(predict(m, Gradient = G, post = post, expected = TRUE))  # rows x columns x draws
    for (h in seq_len(nrow(G$XDataNew))) {
      host <- as.character(G$XDataNew$HostType[h])
      draws <- p[h, , , drop = TRUE]; rownames(draws) <- m$spNames
      if (model_type == "pa") {
        P <- draws
        out <- data.frame(Genus = genus, P_present = rowMeans(P),
                          q2.5 = apply(P, 1, quantile, 0.025), q97.5 = apply(P, 1, quantile, 0.975))
        rich <- mean(colSums(P))
        draws_list[[paste(host, s, sep = "|")]] <- list(P_present = P)
      } else if (model_type == "abund") {
        MU <- draws
        out <- data.frame(Genus = genus, E_log1p_if_present = rowMeans(MU),
                          q2.5 = apply(MU, 1, quantile, 0.025), q97.5 = apply(MU, 1, quantile, 0.975))
        rich <- NA
        draws_list[[paste(host, s, sep = "|")]] <- list(E_log1p_if_present = MU)
      } else {
        gen_ab <- genus[ab]
        has_pa <- paste0(gen_ab, "_PA") %in% m$spNames
        P  <- matrix(1, length(gen_ab), ncol(draws))
        P[has_pa, ] <- draws[paste0(gen_ab[has_pa], "_PA"), ]
        MU <- draws[paste0(gen_ab, "_ABUC"), , drop = FALSE]
        H  <- P * MU
        out <- data.frame(Genus = gen_ab, P_present = rowMeans(P), E_log1p_if_present = rowMeans(MU),
                          E_hurdle = rowMeans(H), q2.5 = apply(H, 1, quantile, 0.025),
                          q97.5 = apply(H, 1, quantile, 0.975))
        rich <- mean(colSums(draws[pa, , drop = FALSE]) + sum(!has_pa))
        draws_list[[paste(host, s, sep = "|")]] <- list(P_present = P, E_log1p_if_present = MU,
                                                        E_hurdle = H, genus = gen_ab)
      }
      rows[[length(rows) + 1]] <- data.frame(HostType = host, Season = s, out, Expected_richness = rich)
    }
  }
  pred_tab <- do.call(rbind, rows)
  write.csv(pred_tab, file.path(outdir, paste0("predictions_", model_type, ".csv")), row.names = FALSE)
  saveRDS(draws_list, file.path(outdir, paste0("predictions_", model_type, "_draws.rds")))
  if (model_type != "abund") {
    cat("\nExpected genus richness by HostType x Season:\n")
    print(unique(pred_tab[, c("HostType", "Season", "Expected_richness")]), row.names = FALSE)
  }
}

writeLines(capture.output(sessionInfo()), file.path(outdir, paste0("sessionInfo_step2_", model_type, ".txt")))
cat("\n=== STEP 2 (", model_type, "|", task, ") finished ===\nEnd time:",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
