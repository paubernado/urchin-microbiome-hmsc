#!/usr/bin/env bash
# =========================================================
# QIIME2 PIPELINE — Arbacia lixula & Paracentrotus lividus microbiota
# Raw paired-end reads -> ASVs, taxonomy, phylogeny, phyloseq export
# Requires: conda env qiime2-2022.11 (see environment.yml)
# =========================================================
set -euo pipefail

conda activate qiime2-2022.11

FASTQ_DIR="${FASTQ_DIR:-./data/sequences}"
cd "$FASTQ_DIR" || { echo "Folder not found: $FASTQ_DIR"; exit 1; }

# ---------------------------------------------------------
# Build manifest
# ---------------------------------------------------------
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

for fwd in *_R1.fastq.gz; do
  sample=${fwd%_R1.fastq.gz}
  rev="${sample}_R2.fastq.gz"
  if [[ -f "$rev" ]]; then
    printf "%s\t%s\t%s\n" "$sample" "$(realpath "$fwd")" "$(realpath "$rev")" >> manifest.tsv
  else
    echo "Warning: reverse file missing for sample $sample"
  fi
done
echo "Manifest created as manifest.tsv in $FASTQ_DIR"

# ---------------------------------------------------------
# STEP 1: Import raw reads into QIIME2
# ---------------------------------------------------------
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux-seqs.qza

# ---------------------------------------------------------
# STEP 2: Trim primers (16S V3-V4, 341F/805R)
# Methods text: "trimmed with Cutadapt, allowing no mismatches"
# -> --p-error-rate 0
# ---------------------------------------------------------
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-seqs.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqs.qza

# ---------------------------------------------------------
# STEP 3: Denoise and cluster into ASVs (DADA2)
# Methods text: "denoised using DADA2, and truncated at 250 bp
# (forward reads) and 240 bp (reverse reads)"
# -> --p-trunc-len-f 250, --p-trunc-len-r 240
# These two lines are the single source of truth for the
# truncation lengths reported in the manuscript's Methods
# section; if they ever change, update here first.
# ---------------------------------------------------------
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 240 \
  --p-n-threads 4 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2-summary.qzv \
  --m-sample-metadata-file ../metadata_pbr2.tsv

# ---------------------------------------------------------
# STEP 4: Phylogeny (MAFFT alignment -> mask -> FastTree -> midpoint root)
# ---------------------------------------------------------
qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment alignment-rep-seqs.qza

qiime alignment mask \
  --i-alignment alignment-rep-seqs.qza \
  --o-masked-alignment masked-alignment-rep-seqs.qza

qiime phylogeny fasttree \
  --i-alignment masked-alignment-rep-seqs.qza \
  --o-tree not-rooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree not-rooted-tree.qza \
  --o-rooted-tree rooted-tree-rep-seqs.qza

# ---------------------------------------------------------
# STEP 5: Taxonomy (SILVA 138 classifier)
# ---------------------------------------------------------
wget -O "silva-138-99-nb-classifier.qza" \
  "https://data.qiime2.org/2022.11/common/silva-138-99-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 1 \
  --p-confidence 0.7 \
  --p-read-orientation same

# ---------------------------------------------------------
# STEP 6: Filter Archaea / Eukaryota / singletons
# ---------------------------------------------------------
qiime taxa filter-table \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Archaea \
  --p-mode contains \
  --o-filtered-table table-dada2-noAr.qza

qiime taxa filter-table \
  --i-table table-dada2-noAr.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Eukaryota \
  --p-mode contains \
  --o-filtered-table table-dada2-noAr-noEu.qza

qiime feature-table filter-features \
  --i-table table-dada2-noAr-noEu.qza \
  --p-min-frequency 2 \
  --o-filtered-table table-clean.qza

# ---------------------------------------------------------
# STEP 7: Rarefaction curve + taxa barplot + diversity
# ---------------------------------------------------------
qiime diversity alpha-rarefaction \
  --i-table table-dada2.qza \
  --i-phylogeny rooted-tree-rep-seqs.qza \
  --p-max-depth 7000 \
  --m-metadata-file ../metadata_pbr2.tsv \
  --o-visualization alpha-rarefaction.qzv

qiime taxa barplot \
  --i-table table-clean.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ../metadata_pbr2.tsv \
  --o-visualization taxa-bar-plots.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-rep-seqs.qza \
  --i-table table-clean.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file ../metadata_pbr2.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../metadata_pbr2.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file ../metadata_pbr2.tsv \
  --o-visualization core-metrics-results/bray_curtis_pcoa_results.qzv

# ---------------------------------------------------------
# STEP 8: Export to phyloseq (unrooted tree + taxonomy + table + fasta)
# ---------------------------------------------------------
mkdir -p Phyloseq

qiime tools export --input-path not-rooted-tree.qza --output-path Phyloseq
qiime tools export --input-path taxonomy.qza --output-path Phyloseq
qiime tools export --input-path table-clean.qza --output-path Phyloseq

sed 's/Feature ID/#OTUID/' Phyloseq/taxonomy.tsv \
  | sed 's/Taxon/taxonomy/' \
  | sed 's/Confidence/confidence/' \
  > Phyloseq/biom-taxonomy.tsv

biom add-metadata \
  -i Phyloseq/feature-table.biom \
  -o Phyloseq/table-with-taxonomy.biom \
  --observation-metadata-fp Phyloseq/biom-taxonomy.tsv \
  --sc-separated taxonomy

qiime tools export \
  --input-path rep-seqs-dada2.qza \
  --output-path Phyloseq

echo "Phyloseq export complete: $FASTQ_DIR/Phyloseq"
echo "Files needed by R/01_import_decontam.R:"
echo "  Phyloseq/table-with-taxonomy.biom"
echo "  Phyloseq/tree.nwk (rename the exported .nwk from not-rooted-tree.qza)"
echo "  Phyloseq/dna-sequences.fasta"
echo "  ../metadata_pbr2.tsv (as metadata_pbr2.tsv, see R script)"

# ---------------------------------------------------------
# STEP 9 (optional): FAPROTAX functional annotation
# ---------------------------------------------------------
# See data/README.md for how to obtain FAPROTAX.txt and collapse_table.py
#
# biom convert -i Phyloseq/table-with-taxonomy.biom \
#   -o Phyloseq/table-with-taxonomy.txt --to-tsv --header-key taxonomy
#
# python3 collapse_table.py \
#   -i Phyloseq/table-with-taxonomy.txt -g FAPROTAX.txt -f \
#   -o functional_otu_table.tsv -r report.txt \
#   --column_names_are_in last_comment_line --keep_header_comments \
#   --non_numeric consolidate -v --row_names_are_in_column "taxonomy" \
#   --omit_columns 0 --normalize_collapsed columns_before_collapsing \
#   --group_leftovers_as 'other'
