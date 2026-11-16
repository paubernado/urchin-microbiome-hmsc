#!/usr/bin/env Rscript
# =============================================================================
# HMSC PIPELINE - STEP 0: DATA PREPARATION (shared by ALL models)
# =============================================================================
# The presence/absence, conditional abundance and joint models all read the
# files produced here, so the three models use exactly the same samples,
# genera, covariates and tree.
#
# BLOCKS
#   1  Import (BIOM + tree + metadata + sequences); "Species" -> "HostType"
#   2  Decontam (prevalence method, 3 negative controls)
#   3  Taxonomy cleanup; keep field samples (urchins + water)
#   4  Sequencing depth = log10 of all reads after decontam
#   5  Genus labels:
#        valid genus                                  -> kept
#        uncultured/unassigned genus, valid family    -> "<Family>_uncultured"
#        neither genus nor family valid               -> excluded from Y
#                                                        (reads still counted
#                                                        in sequencing depth)
#        genus name shared by >1 family               -> "<Family>_<Genus>"
#      Archaeal ASVs (if present) are kept apart and used ONLY to root the tree.
#   6  Genus phylogenetic cohesion (diagnostic)
#   7  Aggregate ASVs to genus (base R) and build XData
#   8  Prevalence filter: >= 5 reads in >= 10% of samples
#   9  Genus-level tree: mean patristic distance between ASVs + neighbour
#      joining, then ROOTED (Archaea outgroup if available, else midpoint) and
#      made ULTRAMETRIC (ape::chronos); the outgroup is dropped afterwards
#   10 Final checks and save
#
# OUTPUT (in --outdir):
#   Y_counts.rds         samples x genera read counts (columns = tree tips)
#   XData.rds            HostType, Season, Site (site names), logDepth, logDepth_z
#   tree_genus_rooted.rds
#   data_settings.rds    filters used (read by STEP 1)
#   genus_phylogenetic_cohesion.csv, genus_prevalence_filter.csv
#
# USAGE (from the folder with the QIIME2 files, hmsc_env active):
#   Rscript HMSC_STEP0_prepare_data.R --outdir results_data
#   Optional: --exclude "Cutibacterium,Lawsonella"  (e.g. suspected contaminants)
# =============================================================================

options(stringsAsFactors = FALSE, width = 200)
invisible(Sys.setlocale("LC_COLLATE", "C"))   # same sort order on laptop and cluster

cat("=== HMSC pipeline - STEP 0: data preparation ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

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

biom_file  <- get_arg_value("--biom",  "table-with-taxonomy.biom")
tree_file  <- get_arg_value("--tree",  "tree.nwk")
meta_file  <- get_arg_value("--meta",  "metadata_pbr.tsv")
fasta_file <- get_arg_value("--fasta", "dna-sequences.fasta")
outdir     <- get_arg_value("--outdir", "results_data")
min_reads  <- as.integer(get_arg_value("--min-reads", "5"))      # presence = >= 5 reads
min_prev   <- as.numeric(get_arg_value("--min-prev", "0.10"))    # >= 10% of samples
exclude_genera <- get_arg_value("--exclude", "")
exclude_genera <- if (exclude_genera == "") character(0) else trimws(strsplit(exclude_genera, ",")[[1]])
make_ultrametric <- get_arg_value("--ultrametric", "yes") == "yes"   # chronos (see Bloc 9)

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Output dir:", outdir, "| presence = >=", min_reads, "reads | prevalence >=", min_prev,
    "| excluded genera:", if (length(exclude_genera)) paste(exclude_genera, collapse = ", ") else "none", "\n")

# -----------------------------
# Packages
# -----------------------------
# phangorn is only needed for midpoint rooting (when there are no Archaea)
required_pkgs <- c("phyloseq", "biomformat", "qiime2R", "decontam", "Biostrings", "ape")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them in hmsc_env.")
}
suppressPackageStartupMessages({
  library(phyloseq); library(qiime2R); library(biomformat)
  library(Biostrings); library(decontam); library(ape)
})

tax_as_matrix <- function(ps) {          # taxonomy table keeping ASV IDs as rownames
  m <- as(tax_table(ps), "matrix"); rownames(m) <- taxa_names(ps); m
}
bad_label <- function(x) {               # labels that do not identify a real taxon
  is.na(x) | x %in% c("", "NA", "Unassigned", "Unknown") |
    grepl("uncultured|metagenome|unidentified|Incertae", x, ignore.case = TRUE)
}

# =============================================================================
# Bloc 1: Import
# =============================================================================
cat("\n--- Bloc 1: Import ---\n")
phylo_raw <- merge_phyloseq(import_biom(BIOMfilename = biom_file, treefilename = tree_file),
                            import_qiime_sample_data(mapfilename = meta_file))
colnames(tax_table(phylo_raw)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phylo_raw <- merge_phyloseq(phylo_raw, refseq(readDNAStringSet(fasta_file)))

# Sample type is stored in the metadata column "Species" -> copy to "HostType"
sd_raw <- sample_data(phylo_raw)
if (!"HostType" %in% colnames(sd_raw)) {
  if (!"Species" %in% colnames(sd_raw)) stop("Neither 'HostType' nor 'Species' in the metadata.")
  sd_raw$HostType <- as.character(sd_raw$Species)
  sample_data(phylo_raw) <- sd_raw
}
print(table(sample_data(phylo_raw)$HostType, useNA = "ifany"))

# =============================================================================
# Bloc 2: Decontam (prevalence method)
# =============================================================================
cat("\n--- Bloc 2: Decontam ---\n")
neg_labels <- c("Extraccion control", "PCR control", "Control")
sample_data(phylo_raw)$is.neg <- sample_data(phylo_raw)$HostType %in% neg_labels
if (sum(sample_data(phylo_raw)$is.neg) != 3) {
  stop("Expected 3 negative controls, found ", sum(sample_data(phylo_raw)$is.neg))
}
contamdf.prev <- isContaminant(phylo_raw, method = "prevalence", neg = "is.neg", threshold = 0.2)
cat("Contaminant ASVs removed:", sum(contamdf.prev$contaminant), "\n")
phylo_raw <- prune_taxa(!contamdf.prev$contaminant, phylo_raw)

# =============================================================================
# Bloc 3: Taxonomy cleanup + field samples
# =============================================================================
tax_raw <- as.data.frame(tax_as_matrix(phylo_raw), stringsAsFactors = FALSE)
tax_raw[] <- lapply(tax_raw, function(x) trimws(gsub("[*~]", "", gsub("^[a-z]__", "", as.character(x)))))
colnames(tax_raw)[colnames(tax_raw) == "Kingdom"] <- "Domain"
tax_mat_raw <- as.matrix(tax_raw); rownames(tax_mat_raw) <- taxa_names(phylo_raw)
tax_table(phylo_raw) <- tax_table(tax_mat_raw)

cat("\n--- Bloc 3: Field samples ---\n")
psJ <- subset_samples(phylo_raw,
                      HostType %in% c("Arbacia lixula", "Paracentrotus lividus", "Water sample") & !is.neg)
psJ <- prune_taxa(taxa_sums(psJ) > 0, psJ)
cat("Samples:", nsamples(psJ), "| ASVs:", ntaxa(psJ), "\n")

# =============================================================================
# Bloc 4: Sequencing depth (all bacterial reads after decontam, incl. unclassified ASVs)
# =============================================================================
# Archaeal ASVs (if present) are used only to root the tree, so they are not
# counted in the sequencing depth: the depth is then the same whether or not
# Archaea were re-exported from QIIME2.
is_archaea <- as(tax_table(psJ), "matrix")[, "Domain"] %in% "Archaea"
logDepth_all <- log10(sample_sums(prune_taxa(!is_archaea, psJ)) + 1)

# =============================================================================
# Bloc 5: Genus labels (+ Archaea kept apart for rooting)
# =============================================================================
cat("\n--- Bloc 5: Genus labels ---\n")
tree_asv <- read.tree(tree_file)
tax0 <- as.data.frame(tax_as_matrix(psJ), stringsAsFactors = FALSE)
common_asv <- intersect(tree_asv$tip.label, rownames(tax0))
if (length(common_asv) < 100) stop("Too few matching ASVs between tree and taxonomy.")
tree_asv <- keep.tip(tree_asv, common_asv)
tax0     <- tax0[common_asv, , drop = FALSE]

# Archaea: used ONLY as outgroup to root the genus tree (Bloc 9), never as
# responses. They were filtered out in the bioinformatic pipeline; if archaeal
# ASVs are added back to the QIIME2 table and tree, they are used automatically.
archaea_asv <- rownames(tax0)[tax0$Domain == "Archaea"]
cat("Archaeal ASVs available for rooting:", length(archaea_asv), "\n")
tax0 <- tax0[!rownames(tax0) %in% archaea_asv, , drop = FALSE]

g_raw <- as.character(tax0$Genus)
fam   <- as.character(tax0$Family)
genus_bad <- bad_label(g_raw); fam_ok <- !bad_label(fam)
genus_lab <- g_raw
genus_lab[genus_bad & fam_ok]  <- paste0(fam[genus_bad & fam_ok], "_uncultured")
genus_lab[genus_bad & !fam_ok] <- NA
multi_fam <- tapply(fam, genus_lab, function(f) length(unique(na.omit(f))))
multi_fam <- names(multi_fam)[multi_fam > 1]
idx <- !is.na(genus_lab) & genus_lab %in% multi_fam
genus_lab[idx] <- paste(fam[idx], genus_lab[idx], sep = "_")
cat("ASVs as <Family>_uncultured:", sum(genus_bad & fam_ok),
    "| ASVs excluded (no valid genus/family):", sum(genus_bad & !fam_ok), "\n")

genus_vec <- setNames(genus_lab, rownames(tax0))   # ASV ID -> genus label
genus_vec <- genus_vec[!is.na(genus_vec)]
reads_excl <- 1 - sum(taxa_sums(psJ)[names(genus_vec)]) / sum(taxa_sums(psJ))
cat("Proportion of reads excluded from Y:", round(reads_excl, 3), "\n")
tree_asv_genus <- keep.tip(tree_asv, names(genus_vec))
genus_vec      <- genus_vec[tree_asv_genus$tip.label]

# =============================================================================
# Bloc 6: Genus phylogenetic cohesion (diagnostic)
# =============================================================================
cat("\n--- Bloc 6: Genus phylogenetic cohesion ---\n")
d <- cophenetic(tree_asv_genus)
check_genus <- do.call(rbind, lapply(sort(unique(genus_vec)), function(g) {
  tips_g <- names(genus_vec)[genus_vec == g]
  if (length(tips_g) < 2) return(data.frame(Genus = g, n_ASV = length(tips_g), ratio = NA_real_, monophyletic = NA))
  w <- d[tips_g, tips_g]
  data.frame(Genus = g, n_ASV = length(tips_g),
             ratio = mean(w[upper.tri(w)]) / mean(d[tips_g, setdiff(names(genus_vec), tips_g)]),
             monophyletic = is.monophyletic(tree_asv_genus, tips_g))
}))
rm(d); invisible(gc())
is_unc <- grepl("_uncultured$", check_genus$Genus)
cat("Monophyletic - named genera:", round(mean(check_genus$monophyletic[!is_unc], na.rm = TRUE), 3),
    "| <Family>_uncultured:", round(mean(check_genus$monophyletic[is_unc], na.rm = TRUE), 3),
    "| median within/between distance ratio:", round(median(check_genus$ratio, na.rm = TRUE), 3), "\n")
write.csv(check_genus, file.path(outdir, "genus_phylogenetic_cohesion.csv"), row.names = FALSE)

# =============================================================================
# Bloc 7: Aggregate to genus and build XData
# =============================================================================
cat("\n--- Bloc 7: Genus aggregation and XData ---\n")
otu <- as(otu_table(psJ), "matrix")
if (taxa_are_rows(psJ)) otu <- t(otu)                            # samples x ASVs
otu <- otu[, names(genus_vec), drop = FALSE]
Y_counts <- t(rowsum(t(otu), group = genus_vec[colnames(otu)]))  # samples x genera

meta <- data.frame(sample_data(psJ), check.names = FALSE)[rownames(Y_counts), , drop = FALSE]
XData <- data.frame(
  HostType = factor(meta$HostType, levels = c("Water sample", "Arbacia lixula", "Paracentrotus lividus")),
  Season   = factor(meta$Season,   levels = c("Spring", "Summer", "Autumn", "Winter")),
  Location = factor(meta$Location, levels = c("Escala", "Blanes")),
  Contamination = factor(meta$Contamination, levels = c("Non contaminated", "Contaminated")),
  row.names = rownames(meta))
XData$Site <- droplevels(interaction(XData$Location, XData$Contamination, sep = "_"))
# Site names -- CHECK this mapping against the field records
site_names <- c("Blanes_Contaminated" = "Blanes_Bay", "Blanes_Non contaminated" = "Cala_Sant_Francesc",
                "Escala_Contaminated" = "Escala_Bay", "Escala_Non contaminated" = "Illa_Mateua")
levels(XData$Site) <- unname(site_names[levels(XData$Site)])
XData$logDepth   <- as.numeric(logDepth_all[rownames(XData)])
XData$logDepth_z <- as.numeric(scale(XData$logDepth))

keep_samples <- rownames(XData)[complete.cases(XData[, c("HostType", "Season", "Site", "logDepth_z")])]
Y_counts <- Y_counts[keep_samples, , drop = FALSE]
XData    <- XData[keep_samples, , drop = FALSE]
stopifnot(all(rownames(Y_counts) == rownames(XData)),
          nlevels(XData$HostType) == 3, nlevels(XData$Site) == 4, !anyNA(XData$Site))
print(table(XData$HostType, XData$Season))
print(table(XData$HostType, XData$Site))

# =============================================================================
# Bloc 8: Prevalence filter (>= 5 reads in >= 10% of samples)
# =============================================================================
cat("\n--- Bloc 8: Prevalence filter ---\n")
min_samples <- ceiling(min_prev * nrow(Y_counts))
n_present   <- colSums(Y_counts >= min_reads)
retained    <- n_present >= min_samples & !colnames(Y_counts) %in% exclude_genera
write.csv(data.frame(Genus = colnames(Y_counts), N_present = n_present,
                     Prevalence = round(n_present / nrow(Y_counts), 3),
                     Excluded = colnames(Y_counts) %in% exclude_genera, Retained = retained),
          file.path(outdir, "genus_prevalence_filter.csv"), row.names = FALSE)
Y_counts <- Y_counts[, retained, drop = FALSE]
cat("Samples:", nrow(Y_counts), "| genera retained:", ncol(Y_counts),
    "(", sum(grepl("_uncultured$", colnames(Y_counts))), "<Family>_uncultured )\n")

# =============================================================================
# Bloc 9: Genus-level tree (NJ) + rooting
# =============================================================================
# The NJ tree is unrooted, and the phylogenetic correlation used by Hmsc is
# measured from the root, so the tree must be rooted:
#   Archaea available -> outgroup rooting, then the outgroup is dropped
#   otherwise         -> midpoint rooting
cat("\n--- Bloc 9: Genus-level tree ---\n")
gv_model <- genus_vec[genus_vec %in% colnames(Y_counts)]
stopifnot(setequal(unique(gv_model), colnames(Y_counts)))
if (length(archaea_asv) > 0) {
  gv_model <- c(gv_model, setNames(rep("ARCHAEA_OUTGROUP", length(archaea_asv)), archaea_asv))
}
tree_model_asv <- keep.tip(tree_asv, names(gv_model))
d_asv  <- cophenetic(tree_model_asv)
genera <- sort(unique(gv_model))
tips   <- lapply(genera, function(g) names(gv_model)[gv_model == g])
Dg <- matrix(0, length(genera), length(genera), dimnames = list(genera, genera))
for (i in seq_along(genera)[-1]) for (j in seq_len(i - 1)) {
  Dg[i, j] <- Dg[j, i] <- mean(d_asv[tips[[i]], tips[[j]], drop = FALSE])   # mean patristic distance
}
rm(d_asv); invisible(gc())
tree_genus <- nj(as.dist(Dg))
cat("Negative branch lengths set to 0:", sum(tree_genus$edge.length < 0), "\n")
tree_genus$edge.length[tree_genus$edge.length < 0] <- 0

has_outgroup <- "ARCHAEA_OUTGROUP" %in% tree_genus$tip.label
if (has_outgroup) {
  tree_genus <- root(tree_genus, outgroup = "ARCHAEA_OUTGROUP", resolve.root = TRUE)
  root_method <- "Archaea outgroup"
} else {
  if (!requireNamespace("phangorn", quietly = TRUE))
    stop("No Archaea found: midpoint rooting needs the R package phangorn (not installed).")
  tree_genus <- phangorn::midpoint(tree_genus)
  root_method <- "midpoint"
}
stopifnot(is.rooted(tree_genus))
cat("Tree rooted by:", root_method, "\n")
tree_nonultra <- tree_genus

# Make the tree ULTRAMETRIC (all tips equidistant from the root), as expected
# under the Brownian-motion model used by Hmsc: penalised-likelihood dating
# (ape::chronos, lambda = 1). Zero-length branches get a tiny positive length
# first. The outgroup (if any) is dropped AFTER dating, so the root stays at
# the Archaea-Bacteria split and the tree remains ultrametric.
if (make_ultrametric) {
  eps <- 1e-6 * median(tree_genus$edge.length[tree_genus$edge.length > 0])
  tree_genus$edge.length[tree_genus$edge.length <= 0] <- eps
  tree_genus <- chronos(tree_genus, lambda = 1, quiet = TRUE)
  class(tree_genus) <- "phylo"
  stopifnot(is.ultrametric(tree_genus))
}
if (has_outgroup) {
  tree_genus    <- drop.tip(tree_genus, "ARCHAEA_OUTGROUP")
  tree_nonultra <- drop.tip(tree_nonultra, "ARCHAEA_OUTGROUP")
}
stopifnot(is.rooted(tree_genus))
cat("Ultrametric:", is.ultrametric(tree_genus), "\n")

# Diagnostic: how much does making the tree ultrametric change the
# phylogenetic correlations used by Hmsc?
C_u <- vcv(tree_genus, corr = TRUE)
C_n <- vcv(tree_nonultra, corr = TRUE)[rownames(C_u), colnames(C_u)]
cor_C <- cor(C_u[upper.tri(C_u)], C_n[upper.tri(C_n)])
cat("Correlation between phylogenetic correlations (ultrametric vs non-ultrametric):",
    round(cor_C, 3), "\n")

# =============================================================================
# Bloc 10: Final checks and save
# =============================================================================
stopifnot(setequal(colnames(Y_counts), tree_genus$tip.label))
Y_counts <- Y_counts[, match(tree_genus$tip.label, colnames(Y_counts)), drop = FALSE]
stopifnot(all(colnames(Y_counts) == tree_genus$tip.label))
cat("OK:", nrow(Y_counts), "samples x", ncol(Y_counts), "genera; tree tips:", Ntip(tree_genus), "\n")

saveRDS(Y_counts,   file.path(outdir, "Y_counts.rds"))
saveRDS(XData,      file.path(outdir, "XData.rds"))
saveRDS(tree_genus, file.path(outdir, "tree_genus_rooted.rds"))
saveRDS(tree_nonultra, file.path(outdir, "tree_genus_rooted_nonultrametric.rds"))
write.tree(tree_genus, file.path(outdir, "tree_genus_rooted.nwk"))
saveRDS(list(min_reads = min_reads, min_prev = min_prev, exclude_genera = exclude_genera,
             root_method = root_method, n_archaea_outgroup = length(archaea_asv),
             ultrametric = is.ultrametric(tree_genus), cor_C_ultra_vs_nonultra = cor_C,
             reads_excluded = reads_excl),
        file.path(outdir, "data_settings.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_step0.txt"))

cat("\n=== STEP 0 finished ===\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
