# Methods ↔ code mapping

This file is the single source of truth linking sentences that will appear
in the manuscript's Methods section to the exact parameter and line of code
that produced them. When the Methods text changes, update the code first
(or vice versa) and then update this table in the same commit, so the two
can never silently drift apart.

Every row below has been checked against the actual pipeline script or R
script; update the "Confirmed?" column if a value here was taken from a
draft rather than a verified run.

## 1. Read trimming and denoising

| Methods sentence | Parameter | Value | Location | Confirmed? |
|---|---|---|---|---|
| "trimmed with Cutadapt, allowing no mismatches" | `--p-error-rate` | `0` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 2 | Yes |
| "denoised using DADA2" | `qiime dada2 denoise-paired` | — | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 3 | Yes |
| "truncated at 250 bp (forward reads)" | `--p-trunc-len-f` | `250` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 3 | Yes |
| "and 240 bp (reverse reads)" | `--p-trunc-len-r` | `240` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 3 | Yes |
| (no additional 5' trimming) | `--p-trim-left-f` / `--p-trim-left-r` | `0` / `0` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 3 | Yes |

Primers removed by cutadapt (16S V3-V4, 341F/805R):
- Forward: `CCTACGGGNGGCWGCAG`
- Reverse: `GACTACHVGGGTATCTAATCC`

## 2. Phylogeny

| Step | Tool | Location |
|---|---|---|
| Multiple sequence alignment | MAFFT (`qiime alignment mafft`) | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 4 |
| Alignment masking | `qiime alignment mask` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 4 |
| Tree building | FastTree (`qiime phylogeny fasttree`) | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 4 |
| Rooting | Midpoint rooting (`qiime phylogeny midpoint-root`) | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 4 |
| Genus-level phylogeny used by HMSC | Neighbor-joining on mean cophenetic distance between ASVs sharing a genus | `R/02_genus_aggregation.R` |

## 3. Taxonomic classification

| Methods sentence | Parameter | Value | Location |
|---|---|---|---|
| Classifier | SILVA 138.1, 99% OTUs, pre-trained naive Bayes | `silva-138-99-nb-classifier.qza` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 5 |
| Confidence threshold | `--p-confidence` | `0.7` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 5 |
| Read orientation | `--p-read-orientation` | `same` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 5 |

> **Note:** QIIME2's data-resources page names this classifier file
> `silva-138-99-nb-classifier.qza` (dropping the point release), but the
> underlying reference database is **SILVA v138.1** — use "SILVA v138.1"
> in the manuscript text, matching the actual release, not the QIIME2
> filename.

## 4. Filtering

| Methods sentence | Parameter | Value | Location |
|---|---|---|---|
| Archaea excluded | `--p-exclude Archaea` | — | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 6 |
| Eukaryota excluded | `--p-exclude Eukaryota` | — | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 6 |
| Chloroplast and mitochondria excluded | `--p-exclude chloroplast,mitochondria` | — | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 6 |
| Singleton / rare features removed | `--p-min-frequency` | `2` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 6 |
| Decontam (prevalence method, negative controls) | `isContaminant(method = "prevalence", threshold = 0.2)` | `0.2` | `R/01_import_decontam.R` |
| Sample-type filter + explicit re-exclusion of negative controls | `subset_samples(Species %in% keep_types & (is.na(is.neg) \| is.neg == FALSE))` | — | `R/01_import_decontam.R` |

## 5. Diversity / rarefaction

| Methods sentence | Parameter | Value | Location |
|---|---|---|---|
| Alpha-rarefaction curve max depth | `--p-max-depth` | `7000` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 7 |
| Rarefaction depth for core diversity metrics (QIIME2) | `--p-sampling-depth` | `1103` | `pipeline/01_qiime2_asv_pipeline.sh`, STEP 7 |
| Rarefaction depth for Figure 2 alpha diversity (R) | `rare_depth` | `5000` | `R/figures/fig2_composition_diversity.R` |

Note the two rarefaction depths above serve different purposes and are
**not** meant to match: the QIIME2 core-metrics depth (1103) reflects the
minimum sequencing depth across all retained samples at that stage of the
pipeline, while the Figure 2 alpha-diversity depth (5000) is applied later,
after further filtering, and is documented inline in the figure script.

## 6. Genus-level aggregation

| Methods sentence | Parameter | Value | Location |
|---|---|---|---|
| Genus aggregation | `microbiome::aggregate_rare(level = "Genus")` | — | `R/02_genus_aggregation.R` |
| Minimum detection | `detection` | `5` | `R/02_genus_aggregation.R` |
| Minimum prevalence | `prevalence` | `0.15` | `R/02_genus_aggregation.R` |
| Genus phylogenetic cohesion diagnostic | ratio of mean within-genus to mean between-genus cophenetic distance | — | `R/01_import_decontam.R` |
| Sequencing depth covariate | `logDepth_z = scale(log10(sample_sums + 1))` | — | `R/02_genus_aggregation.R` |

## 7. HMSC models

See `R/00_mcmc_config.R` for the single source of truth on MCMC settings.

| Model | samples | thin | transient | nChains | Confirmed? |
|---|---|---|---|---|---|
| Presence/absence (with phylogeny) | 12000 | 50 | 4000 | 8 | Yes, from fitted object |
| Conditional abundance (with phylogeny) | 10000 | 20 | 6000 | 8 | Yes, from fitted object |
| Presence/absence (no phylogeny) | 12000 | 50 | 4000 | 8 | No — assumed identical to the phylogeny model for consistency; update `R/00_mcmc_config.R` if your no-phylogeny fitted object differs |
| Conditional abundance (no phylogeny) | 10000 | 20 | 6000 | 8 | No — same caveat as above |
| Null models (intercept-only) | same as corresponding full model | — | — | — | Inherited from `R/00_mcmc_config.R` |

Fixed effects: `~ HostType * Season + Site + logDepth_z`
(`logDepth_z`: log10 sequencing depth per sample, z-scored, computed in
`R/02_genus_aggregation.R`)
Random effect (abundance models only): sample-level random intercept
Distributions: `probit` (presence/absence), `normal` (log1p conditional abundance)
Cross-validation: **4-fold**, `set.seed(123)` (confirmed: matches
`Step2_..._pred_cv.R`, the production script that computes the CV results
reported in the manuscript; `R/03`, `R/04`, and `R/05` all use `nfolds = 4`).

## 8. Figures

| Figure | Script | Output |
|---|---|---|
| Figure 2 (phylum composition + alpha diversity) | `R/figures/fig2_composition_diversity.R` | `figures/Figure2_composition_and_diversity.{pdf,png}` |
| Variance partitioning | `R/figures/fig_variance_partitioning.R` | `figures/FIGURE_VP_both_models.{pdf,png}` |
| Seasonal deviations | `R/figures/fig_seasonal_deviations.R` | `figures/FIGURE_seasonal_deviations_FINAL.{pdf,png}` |
| Host-specific contrasts | `R/figures/fig_host_contrasts.R` | `figures/FIGURE_block4_ABUNDANCE.{pdf,png}`, `figures/FIGURE_block4_PA.{pdf,png}` |

## 9. Software versions and package citations

| Software | Version | Where used |
|---|---|---|
| QIIME2 | 2022.11 | `pipeline/01_qiime2_asv_pipeline.sh`, `environment.yml` |
| R | 4.5 | all `R/*.R` scripts |

| Package | Citation | Where used |
|---|---|---|
| phyloseq | McMurdie & Holmes (2013) | `R/01_import_decontam.R` |
| decontam | Davis et al. (2018) | `R/01_import_decontam.R` |
| microbiome | Lahti & Shetty (2017) | `R/02_genus_aggregation.R` |
| Hmsc | Ovaskainen et al. (2017) | `R/03`, `R/04`, `R/05` |
| coda | Plummer et al. (2006) | convergence diagnostics in `R/03`, `R/04`, `R/05` |

## 10. Models reported in the manuscript

- **Presence/absence and conditional abundance, both with phylogeny**
  (`R/04_hmsc_models_with_phylogeny.R`): reported in the manuscript.
- **Null models** (`R/05_null_models.R`): reported in the manuscript
  (model-vs-null comparison). Methods text for this part is still being
  written.
- **No-phylogeny models** (`R/03_hmsc_models_no_phylogeny.R`): **not**
  part of the manuscript. Kept in the repo for reference only.

---


