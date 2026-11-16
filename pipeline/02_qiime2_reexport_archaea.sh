#!/bin/bash
# =============================================================================
# RE-EXPORT THE QIIME2 DATA KEEPING ARCHAEA (to root the phylogeny)
# =============================================================================
#
# WHY
#   The genus-level tree used by Hmsc must be rooted, because the phylogenetic
#   correlation between two genera is the branch length they share from the
#   ROOT to their common ancestor. The most biologically meaningful root is the
#   split between Archaea and Bacteria, so Archaea are used as an OUTGROUP.
#   Archaea were removed in the original pipeline, so the table, sequences and
#   ASV tree are regenerated here WITHOUT removing them.
#
# WHAT THIS SCRIPT DOES
#   1. Filters the table as in the original pipeline (mitochondria,
#      chloroplasts...) but WITHOUT excluding Archaea
#   2. Keeps only the representative sequences that remain in the table
#   3. Rebuilds the ASV tree (MAFFT + FastTree), now including the Archaea
#   4. Exports the three files read by HMSC_STEP0_prepare_data.R:
#        table-with-taxonomy.biom, tree.nwk, dna-sequences.fasta
#   5. Reports how many archaeal ASVs are available as outgroup
#
# WHAT HAPPENS NEXT (in R, HMSC_STEP0_prepare_data.R, nothing to change)
#   - archaeal ASVs are detected from the taxonomy (Domain == "Archaea") and set
#     apart: they never enter the response matrix Y, and their reads are not
#     counted in the sequencing depth
#   - the genus tree is built as before (mean patristic distance + NJ) with the
#     Archaea as an extra tip "ARCHAEA_OUTGROUP"
#   - the tree is rooted on that tip, made ultrametric (ape::chronos), and the
#     outgroup tip is then dropped (the root stays at the Archaea-Bacteria split)
#   - if no Archaea are found, the tree is midpoint-rooted instead
#
# USAGE
#   cd ~/data/data/microbiota_pbr
#   conda activate qiime2-2022.11      # the QIIME2 version used originally
#   bash QIIME2_reexport_with_archaea.sh
#
#   Original pipeline (from the .qza provenance):
#     table-dada2.qza
#       -> taxa filter-table --p-exclude Archaea,Eukaryota,mitochondria,chloroplast
#       -> table-dada2-clean.qza
#       -> feature-table filter-features --p-min-frequency 2
#       -> table-clean.qza
#   Here the same steps are repeated WITHOUT "Archaea" in the exclude list.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# INPUT FILES (edit these)
# -----------------------------------------------------------------------------
TABLE="table-dada2.qza"       # DADA2 output, before any filter (contains Archaea)
REPSEQS="rep-seqs-dada2.qza"  # DADA2 representative sequences
TAXONOMY="taxonomy.qza"       # taxonomy assigned to the ASVs
METADATA="metadata_pbr.tsv"  # sample metadata (copied to the output folder)

# Taxa removed as in the original pipeline, but WITHOUT "Archaea"
EXCLUDE="Eukaryota,mitochondria,chloroplast"
MIN_FREQ=2                    # filter-features --p-min-frequency, as originally

OUTDIR="qiime_with_archaea"  # new folder: the original files are not overwritten
THREADS=4                    # threads for MAFFT / FastTree

mkdir -p "$OUTDIR"
cd "$OUTDIR"
for f in "$TABLE" "$REPSEQS" "$TAXONOMY" "$METADATA"; do
  [ -f "../$f" ] || { echo "ERROR: ../$f not found"; exit 1; }
done

# -----------------------------------------------------------------------------
# 1. Taxonomic filter WITHOUT removing Archaea
# -----------------------------------------------------------------------------
echo "--- 1. Filtering table (excluding: $EXCLUDE; min frequency: $MIN_FREQ) ---"
qiime taxa filter-table \
  --i-table "../$TABLE" \
  --i-taxonomy "../$TAXONOMY" \
  --p-exclude "$EXCLUDE" \
  --o-filtered-table table-taxfilt.qza

# Same frequency filter as the original table-clean.qza
qiime feature-table filter-features \
  --i-table table-taxfilt.qza \
  --p-min-frequency "$MIN_FREQ" \
  --o-filtered-table table-filt.qza

# -----------------------------------------------------------------------------
# 2. Keep only the sequences of the ASVs that remain in the table
# -----------------------------------------------------------------------------
echo "--- 2. Filtering representative sequences ---"
qiime feature-table filter-seqs \
  --i-data "../$REPSEQS" \
  --i-table table-filt.qza \
  --o-filtered-data rep-seqs-filt.qza

# -----------------------------------------------------------------------------
# 3. Rebuild the ASV tree, now including the archaeal ASVs
#    (the R script roots the GENUS tree itself, so the unrooted tree is used)
# -----------------------------------------------------------------------------
echo "--- 3. Building ASV tree (MAFFT + FastTree) ---"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filt.qza \
  --p-n-threads "$THREADS" \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# -----------------------------------------------------------------------------
# 4. Export the files read by HMSC_STEP0_prepare_data.R
# -----------------------------------------------------------------------------
echo "--- 4. Exporting ---"
qiime tools export --input-path table-filt.qza     --output-path exported
qiime tools export --input-path "../$TAXONOMY"     --output-path exported
qiime tools export --input-path rep-seqs-filt.qza  --output-path exported
qiime tools export --input-path unrooted-tree.qza  --output-path exported

# Taxonomy header in the format expected by biom add-metadata
sed -e '1s/Feature ID/#OTUID/' -e '1s/Taxon/taxonomy/' -e '1s/Confidence/confidence/' \
  exported/taxonomy.tsv > exported/taxonomy_biom.tsv

# BIOM table with the taxonomy attached (same format as the original file)
biom add-metadata \
  -i exported/feature-table.biom \
  -o table-with-taxonomy.biom \
  --observation-metadata-fp exported/taxonomy_biom.tsv \
  --sc-separated taxonomy

cp exported/tree.nwk            tree.nwk
cp exported/dna-sequences.fasta dna-sequences.fasta
cp "../$METADATA"               "$METADATA"

# -----------------------------------------------------------------------------
# 5. Checks
# -----------------------------------------------------------------------------
echo "--- 5. Checks ---"
N_ASV=$(grep -c ">" dna-sequences.fasta)
IDS_TABLE=$(biom summarize-table -i table-with-taxonomy.biom | grep "Num observations" | awk '{print $3}')
# archaeal ASVs that are still in the filtered table
N_ARCH=$(awk -F'\t' 'NR==FNR {keep[$1]; next} ($1 in keep) && $2 ~ /Archaea/' \
  <(grep ">" dna-sequences.fasta | sed 's/>//') exported/taxonomy.tsv | wc -l)

echo "ASVs in sequences: $N_ASV | ASVs in table: $IDS_TABLE"
echo "Archaeal ASVs available as outgroup: $N_ARCH"
if [ "$N_ARCH" -eq 0 ]; then
  echo "WARNING: no Archaea found. STEP 0 will use midpoint rooting."
  echo "         Check that 'Archaea' is not in EXCLUDE and that the input table is the unfiltered one."
fi

echo
echo "Done. Files for HMSC (in $OUTDIR/):"
echo "  table-with-taxonomy.biom  tree.nwk  dna-sequences.fasta  $METADATA"
echo
echo "Next step (from $OUTDIR/, with hmsc_env active):"
echo "  Rscript ../HMSC_STEP0_prepare_data.R --outdir ../results_data"
echo "  (the log should say: 'Tree rooted by: Archaea outgroup')"
