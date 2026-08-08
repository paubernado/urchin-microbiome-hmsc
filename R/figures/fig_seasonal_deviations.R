# =========================================================
# SEASONAL DEVIATIONS
# A-C: Conditional abundance (Water, A. lixula, P. lividus)
# D-F: Presence/absence (Water, A. lixula, P. lividus)
# Criterion: top-15 genera per group by magnitude of seasonal range
#
# Requires (loaded by R/06_figures.R):
#   preds_abundance_with_phylogeny, preds_PA_with_phylogeny,
#   XData_final
#
# WHY THINNING IS USED HERE:
# The raw posterior prediction arrays (preds_abundance_with_phylogeny,
# preds_PA_with_phylogeny) contain every stored MCMC draw
# (samples x nChains = 80,000 and 96,000 draws respectively, per
# the mcmc_config in R/00_mcmc_config.R: abundance = 10000 x 8,
# PA = 12000 x 8). Each full array occupies ~13-16 GB in memory,
# so any operation that copies it (e.g. apply() across all draws)
# can exceed available RAM and crash the R session. We therefore
# compute posterior means from a THINNED subset of draws (every
# 20th draw retained) rather than the full set. This is a purely
# computational thinning step, separate from the thinning already
# applied during MCMC sampling to reduce autocorrelation. Given
# the high effective sample sizes already achieved for both
# models (ESS ~57,000-96,000; see Table SXX), using ~4,000-4,800
# draws for the posterior mean is still far more than sufficient
# for a stable estimate, so this additional thinning does not
# compromise the precision of the reported means.
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

stopifnot(exists("preds_abundance_with_phylogeny"))
stopifnot(exists("preds_PA_with_phylogeny"))
stopifnot(exists("XData_final"))

clean_taxon <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\(.*?\\)", "", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

shorten_taxon <- function(x, max_chars = 30) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

theme_heatmap_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      legend.title = element_text(face = "plain"),
      legend.text = element_text(size = 8),
      plot.title = element_text(face = "plain", hjust = 0.5, size = 10),
      plot.margin = margin(5, 8, 5, 5)
    )
}

# Computes the posterior mean of a large 3D prediction array
# (samples x taxa x MCMC draws) using only a thinned subset of
# draws, processed in small chunks to avoid loading the full
# array into memory at once. See explanation above.
compute_thinned_mean <- function(predY_array, thin_step = 20, chunk_size = 200) {
  n_draws_total <- dim(predY_array)[3]
  draw_idx <- seq(1, n_draws_total, by = thin_step)

  pred_sum <- matrix(0, nrow = dim(predY_array)[1], ncol = dim(predY_array)[2])
  chunks <- split(draw_idx, ceiling(seq_along(draw_idx) / chunk_size))

  for (ch in chunks) {
    block <- predY_array[, , ch, drop = FALSE]
    pred_sum <- pred_sum + apply(block, c(1, 2), sum, na.rm = TRUE)
    rm(block)
    gc()
  }

  pred_mean <- pred_sum / length(draw_idx)
  rownames(pred_mean) <- dimnames(predY_array)[[1]]
  colnames(pred_mean) <- dimnames(predY_array)[[2]]
  pred_mean
}

build_group_season_df <- function(pred_mean, XData_final, value_name = "Value") {

  sample_ids <- rownames(pred_mean)

  meta_df <- XData_final[sample_ids, , drop = FALSE] %>%
    mutate(
      SampleID = rownames(.),
      Group = case_when(
        HostType == "Water sample" ~ "Water sample",
        HostType == "Arbacia lixula" ~ "A. lixula",
        HostType == "Paracentrotus lividus" ~ "P. lividus",
        TRUE ~ NA_character_
      ),
      Group = factor(Group, levels = c("Water sample", "A. lixula", "P. lividus")),
      Season = factor(as.character(Season), levels = c("Spring", "Summer", "Autumn", "Winter"))
    ) %>%
    filter(!is.na(Group))

  pred_long <- as.data.frame(pred_mean) %>%
    mutate(SampleID = rownames(.)) %>%
    pivot_longer(cols = -SampleID, names_to = "Taxon", values_to = "Value") %>%
    left_join(meta_df, by = "SampleID") %>%
    mutate(Taxon_clean = clean_taxon(Taxon))

  pred_long %>%
    group_by(Group, Taxon_clean, Season) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group, Taxon_clean) %>%
    mutate(
      Delta = Value - mean(Value, na.rm = TRUE),
      Seasonal_range = max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)
    ) %>%
    ungroup()
}

make_heatmap <- function(df, group_name, tag, top_n = 15, global_lim, fill_label) {

  group_df <- df %>% filter(Group == group_name)

  top_taxa <- group_df %>%
    group_by(Taxon_clean) %>%
    summarise(r = unique(Seasonal_range)[1], .groups = "drop") %>%
    arrange(desc(r)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon_clean)

  plot_df <- group_df %>% filter(Taxon_clean %in% top_taxa)

  plot_df$Taxon_plot <- factor(
    shorten_taxon(plot_df$Taxon_clean),
    levels = shorten_taxon(top_taxa)
  )

  ggplot(plot_df, aes(x = Season, y = forcats::fct_rev(Taxon_plot), fill = Delta)) +
    geom_tile(width = 0.95, height = 0.95, colour = NA) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-global_lim, global_lim),
      name = fill_label
    ) +
    labs(title = group_name, tag = tag) +
    theme_heatmap_pub(base_size = 9)
}

# --- Conditional abundance (A, B, C) ---

pred_mean_abund <- compute_thinned_mean(preds_abundance_with_phylogeny, thin_step = 20, chunk_size = 200)
pred_group_season_abund <- build_group_season_df(pred_mean_abund, XData_final)
global_lim_abund <- max(abs(pred_group_season_abund$Delta), na.rm = TRUE)

p_A <- make_heatmap(pred_group_season_abund, "Water sample", "A", top_n = 15, global_lim = global_lim_abund, fill_label = expression(Delta))
p_B <- make_heatmap(pred_group_season_abund, "A. lixula", "B", top_n = 15, global_lim = global_lim_abund, fill_label = expression(Delta))
p_C <- make_heatmap(pred_group_season_abund, "P. lividus", "C", top_n = 15, global_lim = global_lim_abund, fill_label = expression(Delta))

rm(pred_mean_abund); gc()

# --- Presence/absence (D, E, F) ---

pred_mean_pa <- compute_thinned_mean(preds_PA_with_phylogeny, thin_step = 20, chunk_size = 200)
pred_group_season_pa <- build_group_season_df(pred_mean_pa, XData_final)
global_lim_pa <- max(abs(pred_group_season_pa$Delta), na.rm = TRUE)

p_D <- make_heatmap(pred_group_season_pa, "Water sample", "D", top_n = 15, global_lim = global_lim_pa, fill_label = expression(Delta~"occ. prob."))
p_E <- make_heatmap(pred_group_season_pa, "A. lixula", "E", top_n = 15, global_lim = global_lim_pa, fill_label = expression(Delta~"occ. prob."))
p_F <- make_heatmap(pred_group_season_pa, "P. lividus", "F", top_n = 15, global_lim = global_lim_pa, fill_label = expression(Delta~"occ. prob."))

rm(pred_mean_pa); gc()

fig_seasonal_deviations <- (p_A + p_B + p_C) / (p_D + p_E + p_F) +
  plot_layout(guides = "collect") &
  theme(plot.tag = element_text(face = "bold", size = 12), legend.position = "right")

print(fig_seasonal_deviations)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/FIGURE_seasonal_deviations_FINAL.pdf", fig_seasonal_deviations, width = 14, height = 9, bg = "white", device = cairo_pdf)
ggsave("figures/FIGURE_seasonal_deviations_FINAL.png", fig_seasonal_deviations, width = 14, height = 9, dpi = 800, bg = "white")
