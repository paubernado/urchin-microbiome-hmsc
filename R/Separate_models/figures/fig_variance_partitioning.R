# =========================================================
# VARIANCE PARTITIONING FIGURE (VP)
# Two panels: conditional abundance model | presence-absence model
# Excludes the Intercept row; bootstrap across taxa for CI
#
# Requires (loaded by R/06_figures.R):
#   VP_abundance_with_phylogeny, VP_PA_with_phylogeny
# =========================================================

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

stopifnot(exists("VP_abundance_with_phylogeny"))
stopifnot(exists("VP_PA_with_phylogeny"))

build_VP_plot <- function(VP_obj, label_map, exclude = "Intercept", plot_title = NULL) {

  vp_taxa <- t(VP_obj$vals)
  vp_taxa <- vp_taxa[, !(colnames(vp_taxa) %in% exclude), drop = FALSE]

  set.seed(1)
  B <- 2000
  boot <- replicate(B, {
    idx <- sample(seq_len(nrow(vp_taxa)), replace = TRUE)
    colMeans(vp_taxa[idx, , drop = FALSE], na.rm = TRUE)
  })

  ci <- t(apply(boot, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  colnames(ci) <- c("lo", "mid", "hi")

  df <- data.frame(
    Predictor = rownames(ci),
    lo = ci[, "lo"], mid = ci[, "mid"], hi = ci[, "hi"],
    row.names = NULL
  )

  df$Label <- label_map[df$Predictor]
  df$Label[is.na(df$Label)] <- df$Predictor[is.na(df$Label)]

  df <- df %>%
    arrange(mid) %>%
    mutate(Label = factor(Label, levels = Label))

  ggplot(df, aes(x = Label, y = mid)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), linewidth = 0.7, colour = "black") +
    geom_label(
      aes(label = paste0(round(mid * 100, 1), "%")),
      nudge_y = 0.02, hjust = 0, size = 3.3, label.size = 0, fill = "white"
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, max(df$hi, na.rm = TRUE) * 1.25)
    ) +
    labs(x = NULL, y = "Explained variance (median ± 95% CI across taxa)", title = plot_title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11))
}

label_map <- c(
  "HostType" = "Host type",
  "Season"   = "Season",
  "Site"     = "Site",
  "HostType:Season" = "Host type × Season",
  "Depth"    = "Sequencing depth",
  "logDepth" = "Sequencing depth"
)

p_vp_abund <- build_VP_plot(
  VP_abundance_with_phylogeny,
  label_map = label_map,
  exclude = c("Intercept", "(Intercept)"),
  plot_title = "Conditional abundance model"
)

p_vp_pa <- build_VP_plot(
  VP_PA_with_phylogeny,
  label_map = label_map,
  exclude = c("Intercept", "(Intercept)"),
  plot_title = "Presence–absence model"
)

fig_VP <- p_vp_abund + p_vp_pa +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

print(fig_VP)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/FIGURE_VP_both_models.pdf", fig_VP, width = 12, height = 4.5, bg = "white")
ggsave("figures/FIGURE_VP_both_models.png", fig_VP, width = 12, height = 4.5, dpi = 800, bg = "white")
