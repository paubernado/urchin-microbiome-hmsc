# =========================================================
# FIGURE 2 — PHYLUM BARPLOTS + ALPHA DIVERSITY
# Columns: A = Water sample, B = A. lixula, C = P. lividus
# Top row: phylum-level relative abundance (Season x Location)
# Rows 2-4: Shannon diversity, Pielou's evenness, Faith's PD
#           (rarefied at 5000 reads; shared y-axis per metric)
#
# Requires (loaded by R/06_figures.R): phylo_raw
# =========================================================

library(phyloseq)
library(microbiome)
library(picante)
library(ape)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(grid)

stopifnot(exists("phylo_raw"))

season_levels <- c("Spring", "Summer", "Autumn", "Winter")

# ---------------------------------------------------------
# PART 1 — Phylum composition barplots
# ---------------------------------------------------------

tax_clean <- as.data.frame(tax_table(phylo_raw), stringsAsFactors = FALSE)
tax_clean[] <- lapply(tax_clean, function(x) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("[*]", "", x)
  x <- gsub("[~]", "", x)
  x
})
tax_clean$Phylum <- ifelse(
  is.na(tax_clean$Phylum) | tax_clean$Phylum %in% c("", "Unassigned"),
  "Unknown", tax_clean$Phylum
)
tax_table(phylo_raw) <- tax_table(as.matrix(tax_clean))

psFg <- tax_glom(phylo_raw, taxrank = "Phylum", NArm = FALSE)

psFg_rel <- microbiome::transform(psFg, "compositional")
df <- psmelt(psFg_rel)
df$Sample_short <- sub("_.*", "", df$Sample)

# Phyla present in < 15% of samples (prevalence) or detected in < 5
# samples are pooled into "Other" before ranking.
prepare_phylum_df <- function(df, species_name,
                              prevalence_thr = 0.15,
                              detection_thr = 5) {

  df_sub <- df %>%
    dplyr::filter(HostType == species_name) %>%
    dplyr::mutate(Season = factor(Season, levels = season_levels))

  n_samples <- length(unique(df_sub$Sample))

  phylum_summary <- df_sub %>%
    dplyr::group_by(Phylum) %>%
    dplyr::summarise(
      prevalence = sum(Abundance > 0) / n_samples,
      detection = sum(Abundance > 0),
      .groups = "drop"
    )

  rare_phyla <- phylum_summary$Phylum[
    phylum_summary$prevalence < prevalence_thr | phylum_summary$detection < detection_thr
  ]

  df_sub %>%
    dplyr::mutate(
      Phylum_plot = as.character(Phylum),
      Phylum_plot = ifelse(is.na(Phylum_plot), "Unknown", Phylum_plot),
      Phylum_plot = ifelse(Phylum_plot %in% rare_phyla, "Other", Phylum_plot)
    )
}

df_water         <- prepare_phylum_df(df, "Water sample")
df_arbacia       <- prepare_phylum_df(df, "Arbacia lixula")
df_paracentrotus <- prepare_phylum_df(df, "Paracentrotus lividus")

all_df <- dplyr::bind_rows(df_water, df_arbacia, df_paracentrotus)

phylum_order <- all_df %>%
  dplyr::filter(!Phylum_plot %in% c("Unknown", "Other")) %>%
  dplyr::group_by(Phylum_plot) %>%
  dplyr::summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_abundance)) %>%
  dplyr::pull(Phylum_plot)

top15_phyla <- phylum_order[1:min(15, length(phylum_order))]

collapse_to_top15 <- function(df_sub, top_phyla) {
  df_sub %>%
    dplyr::mutate(
      Phylum_plot = as.character(Phylum_plot),
      Phylum_plot = ifelse(!Phylum_plot %in% c(top_phyla, "Unknown"), "Other", Phylum_plot)
    )
}

df_water_top15         <- collapse_to_top15(df_water, top15_phyla)
df_arbacia_top15       <- collapse_to_top15(df_arbacia, top15_phyla)
df_paracentrotus_top15 <- collapse_to_top15(df_paracentrotus, top15_phyla)

global_levels <- c(top15_phyla, "Other", "Unknown")

df_water_top15$Phylum_plot         <- factor(df_water_top15$Phylum_plot, levels = global_levels)
df_arbacia_top15$Phylum_plot       <- factor(df_arbacia_top15$Phylum_plot, levels = global_levels)
df_paracentrotus_top15$Phylum_plot <- factor(df_paracentrotus_top15$Phylum_plot, levels = global_levels)

marine_palette <- c(
  "#7F9AA2", "#8AA892", "#A6A38A", "#B89C84", "#B58E7A",
  "#9A94A8", "#9EB4BE", "#B8AE9A", "#8E8179", "#97A89B",
  "#C4B29B", "#A79895", "#B9C3CC", "#9EA083", "#B7A3AD",
  "#F0F0F0", "#787878"
)
global_palette <- marine_palette[1:length(global_levels)]
names(global_palette) <- global_levels

if ("Bacteroidota" %in% names(global_palette))     global_palette["Bacteroidota"]     <- "#7FAF8C"
if ("Proteobacteria" %in% names(global_palette))   global_palette["Proteobacteria"]   <- "#6F88A8"
if ("Fusobacteriota" %in% names(global_palette))   global_palette["Fusobacteriota"]   <- "#8EB6B1"
if ("Firmicutes" %in% names(global_palette))       global_palette["Firmicutes"]       <- "#C97B7B"
if ("Actinobacteriota" %in% names(global_palette)) global_palette["Actinobacteriota"] <- "#C89A64"
if ("Other" %in% names(global_palette))            global_palette["Other"]            <- "#F0F0F0"
if ("Unknown" %in% names(global_palette))          global_palette["Unknown"]          <- "#787878"

theme_pub <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8.5, colour = "black"),
      axis.text.y = element_text(size = 7.2, colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
      strip.text = element_text(face = "bold", size = 7.8),
      panel.spacing.x = unit(0.14, "lines"),
      panel.spacing.y = unit(0.14, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      legend.key.height = unit(0.18, "cm"),
      legend.key.width = unit(0.25, "cm"),
      plot.margin = margin(2, 2, 2, 2),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

plot_phylum <- function(df_sub, palette_map, bar_width = 0.85) {
  ggplot(df_sub, aes(x = Sample_short, y = Abundance, fill = Phylum_plot)) +
    geom_col(width = bar_width, position = position_stack(reverse = TRUE), colour = NA) +
    facet_nested(. ~ Season + Location, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette_map, drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(x = NULL, y = "Relative abundance", fill = "Phylum") +
    theme_pub() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left",
                                label.position = "right", override.aes = list(size = 2.5)))
}

pA <- plot_phylum(df_water_top15, global_palette, bar_width = 0.85) + labs(tag = "A")
pB <- plot_phylum(df_arbacia_top15, global_palette, bar_width = 0.85) + labs(tag = "B")
pC <- plot_phylum(df_paracentrotus_top15, global_palette, bar_width = 0.85) + labs(tag = "C")

# ---------------------------------------------------------
# PART 2 — Alpha diversity (rarefied at 5000 reads)
# ---------------------------------------------------------

# Rarefaction is applied ONLY for alpha diversity estimation here.
# It does NOT affect the count data used to fit the HMSC models,
# which instead account for sequencing depth as a fixed covariate
# (log-transformed and standardized; see Statistical analysis).
set.seed(123)
rare_depth <- 5000

phylo_rarefied <- rarefy_even_depth(
  phylo_raw, sample.size = rare_depth, rngseed = 123,
  replace = FALSE, trimOTUs = TRUE, verbose = TRUE
)

shannon_vals  <- microbiome::diversity(phylo_rarefied, index = "shannon")
richness_vals <- microbiome::richness(phylo_rarefied, index = "observed")

# Columns extracted with $ lose their row names, so they must be
# explicitly re-named before indexing by sample ID.
shannon_named  <- setNames(shannon_vals$shannon, rownames(shannon_vals))
richness_named <- setNames(richness_vals$observed, rownames(richness_vals))
evenness_named <- shannon_named / log(richness_named)

otu_mat <- as(otu_table(phylo_rarefied), "matrix")
if (taxa_are_rows(phylo_rarefied)) otu_mat <- t(otu_mat)

tree <- phy_tree(phylo_rarefied)
common_taxa <- intersect(colnames(otu_mat), tree$tip.label)
otu_mat_pd <- otu_mat[, common_taxa, drop = FALSE]
tree_pd <- ape::keep.tip(tree, common_taxa)

faith_pd <- picante::pd(otu_mat_pd, tree_pd, include.root = TRUE)

meta_df <- data.frame(sample_data(phylo_rarefied), check.names = FALSE, stringsAsFactors = FALSE)
meta_df$Sample <- rownames(meta_df)

div_df <- data.frame(
  Sample = rownames(otu_mat),
  Shannon = shannon_named[rownames(otu_mat)],
  Evenness = evenness_named[rownames(otu_mat)],
  Faith_PD = faith_pd$PD[match(rownames(otu_mat), rownames(faith_pd))]
) %>%
  left_join(meta_df, by = "Sample") %>%
  mutate(Season = factor(Season, levels = season_levels))

add_margin <- function(r, frac = 0.05) {
  pad <- diff(r) * frac
  c(r[1] - pad, r[2] + pad)
}

shannon_range  <- add_margin(range(div_df$Shannon, na.rm = TRUE))
evenness_range <- add_margin(range(div_df$Evenness, na.rm = TRUE))
faithpd_range  <- add_margin(range(div_df$Faith_PD, na.rm = TRUE))

div_water         <- div_df %>% filter(HostType == "Water sample")
div_arbacia       <- div_df %>% filter(HostType == "Arbacia lixula")
div_paracentrotus <- div_df %>% filter(HostType == "Paracentrotus lividus")

host_cols <- c("Water sample" = "#4A7FB5", "Arbacia lixula" = "#55A868", "Paracentrotus lividus" = "#C44E52")

theme_div_notop <- function() {
  theme_classic(base_size = 9) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(3, 5, 3, 5)
    )
}

theme_div_bottom <- function() {
  theme_classic(base_size = 9) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(3, 5, 3, 5)
    )
}

plot_diversity_box <- function(df_sub, metric_col, y_label, y_range, box_colour, bottom_row = FALSE) {
  p <- ggplot(df_sub, aes(x = Season, y = .data[[metric_col]])) +
    geom_boxplot(fill = box_colour, alpha = 0.5, outlier.shape = NA,
                 width = 0.6, linewidth = 0.35, colour = "black") +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.6, colour = "black") +
    coord_cartesian(ylim = y_range) +
    labs(x = NULL, y = y_label)

  if (bottom_row) p + theme_div_bottom() else p + theme_div_notop()
}

p_shannon_water <- plot_diversity_box(div_water, "Shannon", "Shannon diversity", shannon_range, host_cols["Water sample"])
p_shannon_arb   <- plot_diversity_box(div_arbacia, "Shannon", NULL, shannon_range, host_cols["Arbacia lixula"])
p_shannon_par   <- plot_diversity_box(div_paracentrotus, "Shannon", NULL, shannon_range, host_cols["Paracentrotus lividus"])

p_even_water <- plot_diversity_box(div_water, "Evenness", "Pielou's evenness", evenness_range, host_cols["Water sample"])
p_even_arb   <- plot_diversity_box(div_arbacia, "Evenness", NULL, evenness_range, host_cols["Arbacia lixula"])
p_even_par   <- plot_diversity_box(div_paracentrotus, "Evenness", NULL, evenness_range, host_cols["Paracentrotus lividus"])

p_faith_water <- plot_diversity_box(div_water, "Faith_PD", "Faith's PD", faithpd_range, host_cols["Water sample"], bottom_row = TRUE)
p_faith_arb   <- plot_diversity_box(div_arbacia, "Faith_PD", NULL, faithpd_range, host_cols["Arbacia lixula"], bottom_row = TRUE)
p_faith_par   <- plot_diversity_box(div_paracentrotus, "Faith_PD", NULL, faithpd_range, host_cols["Paracentrotus lividus"], bottom_row = TRUE)

fig2_unified <- (pA | pB | pC) /
  (p_shannon_water | p_shannon_arb | p_shannon_par) /
  (p_even_water | p_even_arb | p_even_par) /
  (p_faith_water | p_faith_arb | p_faith_par) +
  plot_layout(heights = c(1.4, 1, 1, 1.2), guides = "collect") &
  theme(legend.position = "bottom")

print(fig2_unified)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure2_composition_and_diversity.pdf", fig2_unified, width = 18, height = 13, bg = "white")
ggsave("figures/Figure2_composition_and_diversity.png", fig2_unified, width = 18, height = 13, dpi = 600, bg = "white")
