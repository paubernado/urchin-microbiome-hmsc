# =========================================================
# GENUS-LEVEL AGGREGATION — build Y (counts) and X (covariates)
# (corresponds to Block 4 of the original analysis notebook)
# =========================================================

library(phyloseq)
library(microbiome)
library(ape)

results_dir <- "results"

psJ_genus <- readRDS(file.path(results_dir, "psJ_genus.rds"))
tree_asv_genus <- readRDS(file.path(results_dir, "tree_asv_genus.rds"))

# ---------------------------------------------------------
# Aggregate to genus
# ---------------------------------------------------------
psJg <- microbiome::aggregate_rare(
  psJ_genus, level = "Genus", detection = 5, prevalence = 0.15
)
psJg <- prune_taxa(taxa_sums(psJg) > 0, psJg)

cat("Genus-level samples:", nsamples(psJg), "\n")
cat("Genus-level taxa:", ntaxa(psJg), "\n")

# ---------------------------------------------------------
# Y: genus x sample count matrix
# ---------------------------------------------------------
Y_counts <- t(as(otu_table(psJg), "matrix"))
if (nrow(Y_counts) != nsamples(psJg)) Y_counts <- t(Y_counts)
Y_counts <- as.matrix(Y_counts)
rownames(Y_counts) <- sample_names(psJg)

keep_taxa <- (colSums(Y_counts) > 0) & (apply(Y_counts, 2, var) > 0)
Y_counts <- Y_counts[, keep_taxa, drop = FALSE]

# ---------------------------------------------------------
# X: covariate table (HostType, Season, Location x Contamination -> Site)
# ---------------------------------------------------------
meta <- data.frame(sample_data(psJg), check.names = FALSE, stringsAsFactors = FALSE)

XData <- data.frame(
  HostType = factor(
    meta$Species,
    levels = c("Water sample", "Arbacia lixula", "Paracentrotus lividus")
  ),
  Season = factor(
    meta$Season,
    levels = c("Spring", "Summer", "Autumn", "Winter")
  ),
  Location = factor(
    meta$Location,
    levels = c("Escala", "Blanes")
  ),
  Contamination = factor(
    meta$Contamination,
    levels = c("Non contaminated", "Contaminated")
  ),
  row.names = rownames(meta)
)

XData$Site <- factor(interaction(XData$Location, XData$Contamination, sep = "_"))

if ("Escala_Non contaminated" %in% levels(XData$Site)) {
  XData$Site <- relevel(XData$Site, ref = "Escala_Non contaminated")
}

# ---------------------------------------------------------
# Sequencing depth (log10, z-scored) — used as a fixed-effect
# covariate in every HMSC model (XFormula includes + logDepth_z)
# to account for sample-to-sample variation in read depth.
# Computed from psJg (post genus-aggregation, pre-genus-variance
# filtering), matching the confirmed production scripts.
# ---------------------------------------------------------
logDepth <- log10(sample_sums(psJg) + 1)

if (any(!is.finite(logDepth))) {
  stop("logDepth contains invalid values.")
}

XData$logDepth_z <- scale(logDepth)[, 1]

vars_needed <- c("HostType", "Season", "Site", "logDepth_z")
keep_samples <- rownames(XData)[complete.cases(XData[, vars_needed, drop = FALSE])]
common_samples <- intersect(rownames(Y_counts), keep_samples)

Y_counts <- Y_counts[common_samples, , drop = FALSE]
XData <- XData[common_samples, , drop = FALSE]

stopifnot(identical(rownames(Y_counts), rownames(XData)))

cat("Final Y_counts dimensions:", paste(dim(Y_counts), collapse = " x "), "\n")

# ---------------------------------------------------------
# Genus-level collapsed phylogenetic tree
# (corresponds to Block 10 of the original notebook)
# ---------------------------------------------------------
tax_asv_for_tree <- as.data.frame(tax_table(psJ_genus), stringsAsFactors = FALSE)
tax_asv_for_tree$Genus <- as.character(tax_asv_for_tree$Genus)

common_tree_taxa <- intersect(tree_asv_genus$tip.label, rownames(tax_asv_for_tree))
tree_asv_genus2 <- keep.tip(tree_asv_genus, common_tree_taxa)
tax_asv_for_tree <- tax_asv_for_tree[common_tree_taxa, , drop = FALSE]

genus_vec2 <- setNames(tax_asv_for_tree$Genus, rownames(tax_asv_for_tree))
genera_in_model <- colnames(Y_counts)
genus_vec2 <- genus_vec2[genus_vec2 %in% genera_in_model]
tree_asv_genus2 <- keep.tip(tree_asv_genus2, names(genus_vec2))

d_asv <- cophenetic(tree_asv_genus2)
genera <- sort(unique(genus_vec2))

Dg <- matrix(
  NA_real_, nrow = length(genera), ncol = length(genera),
  dimnames = list(genera, genera)
)

for (i in seq_along(genera)) {
  for (j in seq_along(genera)) {
    gi <- genera[i]; gj <- genera[j]
    tips_i <- names(genus_vec2)[genus_vec2 == gi]
    tips_j <- names(genus_vec2)[genus_vec2 == gj]
    if (i == j) {
      Dg[i, j] <- 0
    } else {
      Dg[i, j] <- mean(d_asv[tips_i, tips_j, drop = FALSE], na.rm = TRUE)
    }
  }
}

Dg <- as.dist(Dg)
tree_genus <- nj(Dg)

common_genus <- intersect(tree_genus$tip.label, colnames(Y_counts))
tree_genus <- keep.tip(tree_genus, common_genus)

cat("Genus tree tips:", Ntip(tree_genus), "\n")

# ---------------------------------------------------------
# Save outputs for downstream HMSC scripts
# ---------------------------------------------------------
saveRDS(Y_counts, file.path(results_dir, "Y_counts.rds"))
saveRDS(XData, file.path(results_dir, "XData.rds"))
saveRDS(tree_genus, file.path(results_dir, "tree_genus.rds"))
write.tree(tree_genus, file = file.path(results_dir, "Genus_collapsed_tree.nwk"))

cat("\nBlock 4 (genus aggregation) and Block 10 (collapsed tree) complete.\n")
