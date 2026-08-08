# =========================================================
# HOST-SPECIFIC CONTRASTS
# A/D: seasonal heatmap within each host (posterior support >= 0.90)
# B/E: direct host contrast forest plot, directional axis label
# C/F: predicted seasonal trajectories for the top contrasted genera
#      (Water sample shown for visual reference; not tested statistically)
# Top row (A-C): conditional abundance model
# Bottom row (D-F): presence-absence model
#
# Requires (loaded by R/06_figures.R):
#   XData_final,
#   mpost_abundance_with_phylogeny, mpost_PA_with_phylogeny,
#   preds_abundance_with_phylogeny, preds_PA_with_phylogeny
#
# Panel C/F reuse the same thinning strategy explained in
# fig_seasonal_deviations.R (compute_thinned_mean-style logic
# inside build_panelC_from_array), for the same memory-related
# reasons: the full posterior arrays are too large to process in
# full, and the already-high ESS makes a thinned subset of draws
# statistically adequate for posterior summaries.
# =========================================================

library(coda)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)
library(tidyr)
library(purrr)
library(tibble)

stopifnot(exists("XData_final"))
stopifnot(exists("mpost_abundance_with_phylogeny"))
stopifnot(exists("mpost_PA_with_phylogeny"))
stopifnot(exists("preds_abundance_with_phylogeny"))
stopifnot(exists("preds_PA_with_phylogeny"))

clean_taxon_light <- function(x) {
  x <- gsub("_", " ", x)
  trimws(x)
}

shorten_taxon <- function(x, max_chars = 30) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.position = "right",
      legend.title = element_text(face = "plain"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      strip.text = element_text(colour = "black"),
      panel.border = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(4, 8, 4, 4)
    )
}

host_cols_water <- c("Water sample" = "#4A7FB5", "A. lixula" = "#55A868", "P. lividus" = "#C44E52")

posterior_support <- function(x, support_level = 0.90) {
  p_pos <- mean(x > 0, na.rm = TRUE)
  p_neg <- mean(x < 0, na.rm = TRUE)
  if (p_pos >= support_level) return(1)
  if (p_neg >= support_level) return(-1)
  return(0)
}

summarise_samples <- function(x) {
  c(median = median(x), lo = quantile(x, 0.025, names = FALSE), hi = quantile(x, 0.975, names = FALSE))
}

build_panelsAB <- function(mpost_obj, top_n = 12, support_level = 0.90) {

  beta_mat <- as.matrix(mpost_obj$Beta)

  param_info <- data.frame(Parameter = colnames(beta_mat), stringsAsFactors = FALSE)
  parsed <- stringr::str_match(
    param_info$Parameter,
    "^B\\[(.*) \\(C\\d+\\), (.*) \\(S\\d+\\)\\]$"
  )
  param_info$Predictor_raw <- parsed[, 2]
  param_info$Taxon <- gsub("_", " ", parsed[, 3])

  get_beta_samples <- function(predictor_raw, taxon_name) {
    idx <- which(param_info$Predictor_raw == predictor_raw & param_info$Taxon == taxon_name)
    stopifnot(length(idx) == 1)
    as.numeric(beta_mat[, idx])
  }

  taxa_all <- sort(unique(param_info$Taxon))
  season_levels <- c("Summer", "Autumn", "Winter")
  host_map <- data.frame(
    Host_label = c("A. lixula", "P. lividus"),
    Host_raw = c("Arbacia lixula", "Paracentrotus lividus"),
    stringsAsFactors = FALSE
  )

  panelA_df <- purrr::map_dfr(seq_len(nrow(host_map)), function(i) {
    host_label <- host_map$Host_label[i]
    host_raw <- host_map$Host_raw[i]
    purrr::map_dfr(season_levels, function(season) {
      main_pred <- paste0("Season", season)
      int_pred  <- paste0("HostType", host_raw, ":Season", season)
      purrr::map_dfr(taxa_all, function(tx) {
        samp_tot <- get_beta_samples(main_pred, tx) + get_beta_samples(int_pred, tx)
        ss <- summarise_samples(samp_tot)
        supp <- posterior_support(samp_tot, support_level)
        tibble::tibble(
          Taxon = tx, Host = host_label, Season = season,
          median = ss["median"], lo = ss["lo"], hi = ss["hi"],
          supported = supp != 0
        )
      })
    })
  })

  panelB_df <- purrr::map_dfr(taxa_all, function(tx) {
    samp_contrast <- get_beta_samples("HostTypeParacentrotus lividus", tx) -
      get_beta_samples("HostTypeArbacia lixula", tx)
    ss <- summarise_samples(samp_contrast)
    supp <- posterior_support(samp_contrast, support_level)
    tibble::tibble(Taxon = tx, contrast = ss["median"], lo = ss["lo"], hi = ss["hi"], supported = supp != 0)
  })

  rank_df <- panelB_df %>%
    mutate(abs_B = abs(contrast)) %>%
    arrange(desc(supported), desc(abs_B))

  selected_taxa <- head(rank_df$Taxon, top_n)

  panelA_plot <- panelA_df %>%
    filter(Taxon %in% selected_taxa) %>%
    mutate(
      Taxon_short = factor(shorten_taxon(Taxon), levels = shorten_taxon(selected_taxa)),
      Host = factor(Host, levels = c("A. lixula", "P. lividus")),
      Season = factor(Season, levels = season_levels),
      Beta_plot = ifelse(supported, median, NA_real_)
    )

  panelB_plot <- panelB_df %>%
    filter(Taxon %in% selected_taxa) %>%
    arrange(contrast) %>%
    mutate(Taxon_short = factor(shorten_taxon(Taxon), levels = shorten_taxon(Taxon)))

  list(panelA = panelA_plot, panelB = panelB_plot, selected_taxa = selected_taxa)
}

plot_panelA <- function(panelA_plot, y_label = expression(beta)) {
  ggplot(panelA_plot, aes(x = Season, y = forcats::fct_rev(Taxon_short), fill = Beta_plot)) +
    geom_tile(width = 0.95, height = 0.95, colour = NA) +
    facet_grid(. ~ Host, labeller = as_labeller(function(x) paste0(x, "\n(ref: Spring)"))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "white", name = y_label) +
    theme_pub(9) +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 7))
}

plot_panelB <- function(panelB_plot) {
  ggplot(panelB_plot, aes(x = Taxon_short, y = contrast)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "black") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.55, colour = "black") +
    geom_point(size = 2.2, shape = 21, stroke = 0.6, colour = "black", fill = "white") +
    coord_flip() +
    labs(x = NULL, y = "Higher in A. lixula  <\u2014  0  \u2014>  Higher in P. lividus") +
    theme_pub(9) +
    theme(legend.position = "none", axis.text.y = element_text(size = 7))
}

# Same thinning logic as compute_thinned_mean() in
# fig_seasonal_deviations.R, applied here to build seasonal x
# host predicted trajectories directly from the large posterior
# prediction arrays.
build_panelC_from_array <- function(predY_array, XData_final, taxa_keep,
                                     thin_step = 20, chunk_size = 200,
                                     y_axis_label = "Predicted value") {

  stopifnot(is.array(predY_array))
  stopifnot(length(dim(predY_array)) == 3)

  sample_ids <- dimnames(predY_array)[[1]]
  taxa_names <- dimnames(predY_array)[[2]]
  n_draws_total <- dim(predY_array)[3]

  draw_idx <- seq(1, n_draws_total, by = thin_step)
  n_draws <- length(draw_idx)
  cat("Panel C: using", n_draws, "thinned draws out of", n_draws_total, "\n")

  meta_df <- XData_final[sample_ids, , drop = FALSE] %>%
    mutate(
      SampleID = rownames(.),
      Host = case_when(
        HostType == "Water sample" ~ "Water sample",
        HostType == "Arbacia lixula" ~ "A. lixula",
        HostType == "Paracentrotus lividus" ~ "P. lividus",
        TRUE ~ NA_character_
      ),
      Season = factor(as.character(Season), levels = c("Spring", "Summer", "Autumn", "Winter"))
    ) %>%
    filter(!is.na(Host))

  keep <- sample_ids %in% meta_df$SampleID
  sample_ids <- sample_ids[keep]
  meta_df <- meta_df[match(sample_ids, meta_df$SampleID), ]

  taxa_names_clean <- clean_taxon_light(taxa_names)
  taxa_idx <- match(taxa_keep, taxa_names_clean)
  if (any(is.na(taxa_idx))) {
    stop("Some taxa_keep not found in array taxa names. Missing: ",
         paste(taxa_keep[is.na(taxa_idx)], collapse = ", "))
  }

  seasons <- levels(meta_df$Season)
  hosts <- c("Water sample", "A. lixula", "P. lividus")
  n_taxa_keep <- length(taxa_keep)

  value_store <- array(NA_real_, dim = c(n_draws, length(seasons), length(hosts), n_taxa_keep),
                        dimnames = list(NULL, seasons, hosts, taxa_keep))

  chunks <- split(draw_idx, ceiling(seq_along(draw_idx) / chunk_size))
  draw_pos <- 0

  for (ch in chunks) {
    block <- predY_array[keep, taxa_idx, ch, drop = FALSE]
    n_block <- length(ch)

    for (h in hosts) {
      h_rows <- meta_df$Host == h
      if (sum(h_rows) == 0) next
      sub_block <- block[h_rows, , , drop = FALSE]
      sub_season <- meta_df$Season[h_rows]

      for (s in seasons) {
        s_rows <- sub_season == s
        if (sum(s_rows) == 0) next
        vals <- apply(sub_block[s_rows, , , drop = FALSE], c(2, 3), mean, na.rm = TRUE)
        value_store[(draw_pos + 1):(draw_pos + n_block), s, h, ] <- t(vals)
      }
    }
    draw_pos <- draw_pos + n_block
    rm(block)
    gc()
  }

  traj_summary <- purrr::map_dfr(taxa_keep, function(tx) {
    purrr::map_dfr(hosts, function(h) {
      purrr::map_dfr(seasons, function(s) {
        v <- value_store[, s, h, tx]
        tibble::tibble(
          Taxon_clean = tx, Host = h, Season = s,
          median = median(v, na.rm = TRUE),
          lo = quantile(v, 0.025, na.rm = TRUE),
          hi = quantile(v, 0.975, na.rm = TRUE)
        )
      })
    })
  })

  plot_df <- traj_summary %>%
    mutate(
      Taxon_plot = factor(shorten_taxon(Taxon_clean), levels = shorten_taxon(taxa_keep)),
      Host = factor(Host, levels = hosts),
      Season = factor(Season, levels = seasons)
    )

  ggplot(plot_df, aes(x = Season, y = median, group = Host, colour = Host, fill = Host)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.14, linewidth = 0, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ Taxon_plot, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = host_cols_water) +
    scale_fill_manual(values = host_cols_water) +
    labs(x = NULL, y = y_axis_label, colour = "Group", fill = "Group") +
    theme_pub(9) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

resAB_abund <- build_panelsAB(mpost_abundance_with_phylogeny, top_n = 12, support_level = 0.90)
p_A_abund <- plot_panelA(resAB_abund$panelA, y_label = expression(beta))
p_B_abund <- plot_panelB(resAB_abund$panelB)

resAB_pa <- build_panelsAB(mpost_PA_with_phylogeny, top_n = 12, support_level = 0.90)
p_A_pa <- plot_panelA(resAB_pa$panelA, y_label = expression(beta))
p_B_pa <- plot_panelB(resAB_pa$panelB)

p_C_abund <- build_panelC_from_array(
  predY_array = preds_abundance_with_phylogeny,
  XData_final = XData_final,
  taxa_keep = resAB_abund$selected_taxa[1:8],
  thin_step = 20, chunk_size = 200,
  y_axis_label = "Predicted conditional abundance"
)

p_C_pa <- build_panelC_from_array(
  predY_array = preds_PA_with_phylogeny,
  XData_final = XData_final,
  taxa_keep = resAB_pa$selected_taxa[1:8],
  thin_step = 20, chunk_size = 200,
  y_axis_label = "Predicted probability of occurrence"
)

fig_abundance <- (p_A_abund + theme(legend.position = "right")) +
  p_B_abund +
  (p_C_abund + theme(legend.position = "right")) +
  plot_layout(ncol = 3, widths = c(1, 0.7, 2)) +
  plot_annotation(tag_levels = "A")

print(fig_abundance)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/FIGURE_block4_ABUNDANCE.pdf", fig_abundance, width = 20, height = 6, bg = "white")
ggsave("figures/FIGURE_block4_ABUNDANCE.png", fig_abundance, width = 20, height = 6, dpi = 800, bg = "white")

fig_pa <- (p_A_pa + theme(legend.position = "right")) +
  p_B_pa +
  (p_C_pa + theme(legend.position = "right")) +
  plot_layout(ncol = 3, widths = c(1, 0.7, 2)) +
  plot_annotation(tag_levels = "A")

print(fig_pa)

ggsave("figures/FIGURE_block4_PA.pdf", fig_pa, width = 20, height = 6, bg = "white")
ggsave("figures/FIGURE_block4_PA.png", fig_pa, width = 20, height = 6, dpi = 800, bg = "white")
