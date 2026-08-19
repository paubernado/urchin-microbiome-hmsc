# Sea urchin coelomic fluid microbiome: QIIME2 + HMSC pipeline

Bioinformatic and statistical pipeline for the **coelomic fluid microbiome** of the sea urchins *Arbacia lixula* and *Paracentrotus lividus* (with
surrounding seawater samples as environmental controls), from raw amplicon
reads to joint species distribution models (HMSC) and manuscript figures.

## Pipeline overview

[![Pipeline overview](https://github.com/paubernado/urchin-microbiome-hmsc/raw/main/docs/pipeline_overview.svg)](https://github.com/paubernado/urchin-microbiome-hmsc/blob/main/docs/pipeline_overview.svg)

The diagram above shows the pipeline as reported in the manuscript (QIIME2
→ R/phyloseq → HMSC models with phylogeny + matching null models →
figures). The repository also includes a no-phylogeny variant of the HMSC
models (`R/03_hmsc_models_no_phylogeny.R`), kept for reference only — it
is **not** part of the manuscript. See [Repository structure](#repository-structure) and [`METHODS.md`](https://github.com/paubernado/urchin-microbiome-hmsc/blob/main/METHODS.md) for details.

## Repository structure

```
.
├── README.md
├── METHODS.md                   # manuscript Methods ↔ code mapping (single source of truth)
├── install_dependencies.R       # one-shot R package installer
├── environment.yml              # QIIME2 conda environment
├── renv.lock                    # R package versions (see Setup) — generate with renv::snapshot()
├── docs/
│   └── pipeline_overview.svg    # diagram shown above
├── pipeline/
│   └── 01_qiime2_asv_pipeline.sh    # raw reads -> ASVs, taxonomy, tree, phyloseq export
├── R/
│   ├── 00_mcmc_config.R             # single source of truth for all MCMC settings
│   ├── 01_import_decontam.R         # import phyloseq, decontam, taxonomy cleanup, filtering
│   ├── 02_genus_aggregation.R       # aggregate to genus, build Y (counts) and X (covariates incl. logDepth_z)
│   ├── 03_hmsc_models_no_phylogeny.R
│   ├── 04_hmsc_models_with_phylogeny.R
│   ├── 05_null_models.R             # intercept-only null models + model-vs-null table
│   ├── 06_figures.R                 # master script: loads results, sources figures/*.R
│   ├── 07_alpha_diversity_models.R  # Shannon/evenness/Faith's PD: rarefaction (5000 reads,
│   │                                 #   seed=123), LMMs (Shannon, log-Faith's PD) + beta
│   │                                 #   regression (evenness) with Site as random intercept,
│   │                                 #   Type II Wald tests, emmeans pairwise contrasts;
│   │                                 #   writes Table SXX (alpha diversity summary)
│   └── figures/
│       ├── fig2_composition_diversity.R
│       ├── fig3_alpha_diversity.R   # boxplot + significance brackets (host contrast) and
│       │                             #   seasonal reaction-norm panels, reusing the models
│       │                             #   fitted in R/07_alpha_diversity_models.R
│       ├── fig_variance_partitioning.R
│       ├── fig_seasonal_deviations.R
│       └── fig_host_contrasts.R
├── results/                     # .rds objects saved by R/01-05, 07 (not tracked in git; see below)
├── figures/                     # .pdf / .png output of R/06_figures.R
└── data/
    └── README.md                 # how to obtain metadata_pbr2.tsv, raw reads, etc.
```

Large intermediate files (`results/*.rds`, raw reads, reference databases) are
not tracked in this repository. See `data/README.md` for how to regenerate
or download them.

## Setup

### QIIME2 (bash pipeline)

```
conda env create -f environment.yml
conda activate qiime2-2022.11
```

### R (statistical pipeline)

**Quick start** (installs latest versions, no version pinning):

```
source("install_dependencies.R")
```

**Exact reproduction of the paper's package versions** (recommended once
you have a working setup — this captures what you actually have installed):

```
install.packages("renv")
renv::init()
renv::snapshot()   # writes renv.lock; commit it so others can renv::restore()
```

If a `renv.lock` file already exists in the repo (generated this way by the
authors), others can instead run:

```
install.packages("renv")
renv::restore()   # installs the exact versions listed in renv.lock
```

Key R packages: `phyloseq`, `decontam`, `microbiome`, `Hmsc`, `coda`, `picante`, `ape`, `ggplot2`, `patchwork`, `ggh4x`, `qiime2R`, `lme4`, `lmerTest`, `glmmTMB`, `car`, `emmeans`, `multcomp`, `multcompView` (see `install_dependencies.R` for the full list, including Bioconductor and
GitHub-only packages).

## Running the pipeline

1. **QIIME2** (on a machine/cluster with the raw fastq files):

```
bash pipeline/01_qiime2_asv_pipeline.sh
```

    Produces the `Phyloseq/` export folder (biom table, rooted/unrooted tree,
taxonomy, representative sequences) described in `data/README.md`.

2. **R / phyloseq + decontam + genus aggregation**:

```
source("R/01_import_decontam.R")
source("R/02_genus_aggregation.R")
```

3. **HMSC models** (MCMC settings are centralized in `R/00_mcmc_config.R` — edit that single file if you need to change `samples` / `thin` / `transient` / `nChains` / `nParallel`):

```
source("R/03_hmsc_models_no_phylogeny.R")
source("R/04_hmsc_models_with_phylogeny.R")
source("R/05_null_models.R")
```

4. **Alpha diversity models** (Shannon diversity, Pielou's evenness, Faith's PD; independent of the HMSC models above — only needs `results/phylo_raw.rds` from step 2):

```
source("R/07_alpha_diversity_models.R")
```

    Writes `results/TableSXX_alpha_diversity_summary.csv` (descriptive
statistics, fixed-effect significance, and pairwise host/season
contrasts for all three metrics).

5. **Figures**:

```
source("R/06_figures.R")
```

    This loads every object the figure scripts need directly from the `.rds` files saved in steps 3–4 (so figures can be regenerated without
re-running MCMC), then produces, in order:

  - `figures/Figure2_composition_and_diversity.{pdf,png}`
  - `figures/Figure3_alpha_diversity.{pdf,png}`
  - `figures/FIGURE_VP_both_models.{pdf,png}`
  - `figures/FIGURE_seasonal_deviations_FINAL.{pdf,png}`
  - `figures/FIGURE_block4_ABUNDANCE.{pdf,png}` and `FIGURE_block4_PA.{pdf,png}`

## Methods ↔ code mapping

[`METHODS.md`](https://github.com/paubernado/urchin-microbiome-hmsc/blob/main/METHODS.md) links every sentence that will appear in the
manuscript's Methods section to the exact parameter and line of code that
produced it (trimming, DADA2 truncation lengths, classifier confidence,
rarefaction depths, HMSC MCMC settings, alpha diversity model
specification, etc.). Check it before editing either the manuscript text
or the pipeline, so the two cannot silently drift apart.

## MCMC settings (confirmed from fitted model objects)

| Model                                  | samples | thin | transient | nChains |
| --------------------------------------- | ------- | ---- | --------- | ------- |
| Presence/absence (with phylogeny)      | 12000   | 50   | 4000      | 8       |
| Conditional abundance (with phylogeny) | 10000   | 20   | 6000      | 8       |

These values live in one place, `R/00_mcmc_config.R`, and are sourced by
every fitting script (`03`, `04`, `05`) so the manuscript's Methods /
Table SXX and the code cannot drift apart.

## Alpha diversity model specification

| Metric            | Model                                              | Response scale |
| ------------------ | --------------------------------------------------- | -------------- |
| Shannon diversity  | Linear mixed model (`lme4`/`lmerTest`)              | Raw            |
| Faith's PD         | Linear mixed model (`lme4`/`lmerTest`)              | log            |
| Pielou's evenness  | Beta-regression mixed model (`glmmTMB`, logit link) | Raw (0–1)      |

All three: `metric ~ HostType * Season + (1 | Site)`, rarefied to 5,000
reads/sample (seed = 123) before calculating diversity indices. Fixed-effect
significance via Type II Wald chi-square tests (`car::Anova`); pairwise
host/season contrasts via `emmeans`, Benjamini–Hochberg (FDR)-corrected.
Defined in `R/07_alpha_diversity_models.R`.

## Status

This repository accompanies a manuscript currently in preparation. A
citation (and DOI, once available as a preprint or published article) will
be added here.

## License

This project is licensed under the MIT License — see LICENSE for details.
