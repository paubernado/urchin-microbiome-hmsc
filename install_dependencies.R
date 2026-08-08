# =========================================================
# INSTALL ALL R DEPENDENCIES FOR THIS PIPELINE
# Run once before sourcing any script in R/.
#
# For reproducing the EXACT package versions used in the paper,
# prefer renv instead (see README "Setup" section):
#   install.packages("renv")
#   renv::init()
#   renv::snapshot()   # run this yourself once everything works,
#                       # to create renv.lock with your exact versions
# =========================================================

# ---------------------------------------------------------
# CRAN packages
# ---------------------------------------------------------
cran_pkgs <- c(
  "ape",
  "coda",
  "dplyr",
  "forcats",
  "ggh4x",
  "ggplot2",
  "patchwork",
  "picante",
  "purrr",
  "scales",
  "stringr",
  "tibble",
  "tidyr"
)

new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran) > 0) install.packages(new_cran)

# ---------------------------------------------------------
# Bioconductor packages
# ---------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "Biostrings",
  "biomformat",
  "decontam",
  "microbiome",
  "phyloseq"
)

new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc) > 0) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# ---------------------------------------------------------
# GitHub-only packages
# ---------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Hmsc: joint species distribution models
if (!requireNamespace("Hmsc", quietly = TRUE)) {
  remotes::install_github("hmsc-r/HMSC")
}

# qiime2R: import QIIME2 artifacts into R / phyloseq
if (!requireNamespace("qiime2R", quietly = TRUE)) {
  remotes::install_github("jbisanz/qiime2R")
}

cat("\nAll dependencies installed (or already present).\n")
