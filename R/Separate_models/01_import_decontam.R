# =========================================================
# IMPORT, DECONTAM, TAXONOMY CLEANUP, SAMPLE FILTERING,
# GENUS PHYLOGENETIC COHESION CHECK
# (corresponds to Blocks 2-4 of the original analysis notebook)
# =========================================================

library(phyloseq)
library(qiime2R)
library(biomformat)
library(Biostrings)
library(decontam)
library(microbiome)
library(ggplot2)
library(ape)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# ---------------------------------------------------------
# INPUT FILE PATHS -- edit these if your filenames differ.
# ---------------------------------------------------------
biom_file     <- "data/Phyloseq/table-with-taxonomy.biom"
tree_file     <- "data/Phyloseq/tree.nwk"
metadata_file <- "data/metadata_pbr2.tsv"
fasta_file    <- "data/Phyloseq/dna-sequences.fasta"

# ---------------------------------------------------------
# Import biom + tree + metadata
# ---------------------------------------------------------
biom_data_raw <- import_biom(
  BIOMfilename = biom_file,
  treefilename = tree_file
)

mapping_file_raw <- import_qiime_sample_data(
  mapfilename = metadata_file
)

phylo_raw <- merge_phyloseq(biom_data_raw, mapping_file_raw)

colnames(tax_table(phylo_raw)) <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)

rep_seqs <- Biostrings::readDNAStringSet(fasta_file)
phylo_raw <- merge_phyloseq(phylo_raw, refseq(rep_seqs))

cat("Initial samples:", nsamples(phylo_raw), "\n")
cat("Initial taxa:", ntaxa(phylo_raw), "\n")
# Expected: 187 samples (184 biological/field samples + 3 negative
# controls: NCEXT, NCPR, Control), now that metadata_pbr2.tsv
# includes all 3 control rows.

# ---------------------------------------------------------
# Decontam (prevalence method)
# ---------------------------------------------------------
sample_data(phylo_raw)$is.neg <- sample_data(phylo_raw)$HostType %in% c(
  "Extraccion control", "PCR control", "Control"
)

cat("\nNegative control flag counts:\n")
print(table(sample_data(phylo_raw)$is.neg, useNA = "ifany"))
# Expected: FALSE 184, TRUE 3

contamdf.prev <- isContaminant(
  phylo_raw, method = "prevalence", neg = "is.neg", threshold = 0.2
)

cat("\nContaminants detected:\n")
print(table(contamdf.prev$contaminant))

# Presence/absence plot for prevalence diagnostics
ps.pa <- transform_sample_counts(phylo_raw, function(abund) 1 * (abund > 0))
neg_labels <- c("Extraccion control", "PCR control", "Control")
ps.pa.neg <- prune_samples(sample_data(phylo_raw)$HostType %in% neg_labels, ps.pa)
ps.pa.pos <- prune_samples(!(sample_data(phylo_raw)$HostType %in% neg_labels), ps.pa)

df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos),
  pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev$contaminant
)

p_decontam <- ggplot(df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence in negative controls") +
  ylab("Prevalence in true samples")
print(p_decontam)
ggsave(file.path(results_dir, "decontam_prevalence_diagnostic.png"), p_decontam,
       width = 6, height = 5, dpi = 300)

phylo_raw <- prune_taxa(!contamdf.prev$contaminant, phylo_raw)

cat("\nSamples after decontam:", nsamples(phylo_raw), "\n")
cat("Taxa after decontam:", ntaxa(phylo_raw), "\n")

# ---------------------------------------------------------
# Report real sequencing depth of the 3 negative controls
# (QC / methods reporting only -- not used downstream, since
# these samples are excluded from psJ below)
# ---------------------------------------------------------
control_ids <- c("NCEXT", "NCPR", "Control")
control_ids_present <- intersect(control_ids, sample_names(phylo_raw))
if (length(control_ids_present) > 0) {
  depths <- sample_sums(phylo_raw)[control_ids_present]
  cat("\nNegative control read depths (post-decontam):\n")
  print(depths)
}

# ---------------------------------------------------------
# Taxonomy cleanup (runs ONCE, after decontam, directly on
# phylo_raw -- no separate "ps_raw" copy)
# ---------------------------------------------------------
tax_raw <- as.data.frame(tax_table(phylo_raw), stringsAsFactors = FALSE)

tax_raw[] <- lapply(tax_raw, function(x) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("[*]", "", x)
  x <- gsub("[~]", "", x)
  x
})

colnames(tax_raw)[colnames(tax_raw) == "Kingdom"] <- "Domain"
tax_table(phylo_raw) <- as.matrix(tax_raw)

cat("\nTaxonomy cleaned.\n")

# ---------------------------------------------------------
# Filter to samples of interest
# ---------------------------------------------------------
ps0 <- phylo_raw
keep_types <- c("Arbacia lixula", "Paracentrotus lividus", "Water sample")

# Negative controls are excluded again here explicitly (on top of
# the decontam taxa removal above) as a belt-and-braces check.
psJ <- subset_samples(
  ps0,
  HostType %in% keep_types & (is.na(is.neg) | is.neg == FALSE)
)
psJ <- prune_taxa(taxa_sums(psJ) > 0, psJ)

cat("\nSamples after keeping target sample types:", nsamples(psJ), "\n")
cat("Taxa after keeping target sample types:", ntaxa(psJ), "\n")
print(table(sample_data(psJ)$HostType, useNA = "ifany"))

# ---------------------------------------------------------
# Match tree tips with taxonomy-assigned ASVs
# ---------------------------------------------------------
tree_asv <- read.tree(tree_file)
tax0 <- as.data.frame(tax_table(psJ), stringsAsFactors = FALSE)

common_asv <- intersect(tree_asv$tip.label, rownames(tax0))
cat("\nASVs shared between tree and taxonomy:", length(common_asv), "\n")

if (length(common_asv) < 100) {
  stop("Too few matching ASVs between tree and taxonomy.")
}

tree_asv <- keep.tip(tree_asv, common_asv)
tax0 <- tax0[common_asv, , drop = FALSE]

tax0$Genus <- as.character(tax0$Genus)
tax0$Genus[tax0$Genus %in% c("", "NA", "Unassigned")] <- NA

keep_asv_genus <- rownames(tax0)[!is.na(tax0$Genus)]

psJ_genus <- prune_taxa(keep_asv_genus, psJ)
tree_asv_genus <- keep.tip(tree_asv, keep_asv_genus)

cat("Samples after keeping assigned genera:", nsamples(psJ_genus), "\n")
cat("Taxa after keeping assigned genera:", ntaxa(psJ_genus), "\n")
cat("Tree tips after keeping assigned genera:", Ntip(tree_asv_genus), "\n")

# ---------------------------------------------------------
# Genus phylogenetic cohesion diagnostic
# ---------------------------------------------------------
genus_phylo_check <- function(tree, genus_vec) {
  d <- cophenetic(tree)
  genera <- sort(unique(genus_vec))

  out <- lapply(genera, function(g) {
    tips_g <- names(genus_vec)[genus_vec == g]
    tips_other <- names(genus_vec)[genus_vec != g]

    if (length(tips_g) < 2 || length(tips_other) < 1) {
      return(data.frame(
        Genus = g,
        n_ASV = length(tips_g),
        mean_within = NA_real_,
        mean_between = NA_real_,
        ratio_within_between = NA_real_
      ))
    }

    within_mat <- d[tips_g, tips_g, drop = FALSE]
    within_vals <- within_mat[upper.tri(within_mat)]
    between_vals <- as.vector(d[tips_g, tips_other, drop = FALSE])

    data.frame(
      Genus = g,
      n_ASV = length(tips_g),
      mean_within = mean(within_vals, na.rm = TRUE),
      mean_between = mean(between_vals, na.rm = TRUE),
      ratio_within_between = mean(within_vals, na.rm = TRUE) /
        mean(between_vals, na.rm = TRUE)
    )
  })

  do.call(rbind, out)
}

tax_genus_check <- as.data.frame(tax_table(psJ_genus), stringsAsFactors = FALSE)
genus_vec <- setNames(as.character(tax_genus_check$Genus), rownames(tax_genus_check))

check_genus <- genus_phylo_check(tree_asv_genus, genus_vec)

cat("\nSummary of genus phylogenetic cohesion:\n")
print(summary(check_genus$ratio_within_between))

cat("\nProportion genera with ratio < 0.5:\n")
print(mean(check_genus$ratio_within_between < 0.5, na.rm = TRUE))

cat("\nProportion genera with ratio < 0.75:\n")
print(mean(check_genus$ratio_within_between < 0.75, na.rm = TRUE))

cat("\nProportion genera with ratio < 1:\n")
print(mean(check_genus$ratio_within_between < 1, na.rm = TRUE))

saveRDS(check_genus, file.path(results_dir, "genus_phylo_cohesion_check.rds"))

# ---------------------------------------------------------
# Save outputs for downstream scripts
# ---------------------------------------------------------
saveRDS(phylo_raw, file.path(results_dir, "phylo_raw.rds"))
saveRDS(psJ, file.path(results_dir, "psJ.rds"))
saveRDS(psJ_genus, file.path(results_dir, "psJ_genus.rds"))
saveRDS(tree_asv_genus, file.path(results_dir, "tree_asv_genus.rds"))

cat("\nBlock 2-4 (import, decontam, taxonomy cleanup, filtering,",
    "phylogenetic cohesion check) complete.\n")

cat("\nBlock 2-4 (import, decontam, taxonomy cleanup, filtering,",
    "phylogenetic cohesion check) complete.\n")
