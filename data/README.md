# Data

Raw sequencing data and large intermediate files are not tracked in this
repository. This folder documents what is needed to reproduce the pipeline
end to end.

## Raw reads

- Paired-end 16S amplicon fastq files (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
- Deposit / accession: [add SRA/ENA accession once available].
- Expected layout for `pipeline/01_qiime2_asv_pipeline.sh`:
  ```
  sequences/
  ├── sample1_R1.fastq.gz
  ├── sample1_R2.fastq.gz
  ├── sample2_R1.fastq.gz
  └── ...
  ```

## Reference database

- SILVA 138 pre-trained naive Bayes classifier (QIIME2 2022.11 format):
  ```
  wget -O silva-138-99-nb-classifier.qza \
    "https://data.qiime2.org/2022.11/common/silva-138-99-nb-classifier.qza"
  ```

## Metadata

- `metadata_pbr2.tsv`: sample metadata (HostType/Species, Season, Location,
  Contamination, sequencing depth, negative-control flags, etc.), one row
  per sample, tab-separated, first column `sample-id` matching the fastq
  manifest. Expected at `data/metadata_pbr2.tsv`.
- Add the actual file (or a link to where it is archived, e.g. a data
  repository / supplementary material) here.

## Phyloseq export (QIIME2 output consumed by `R/01_import_decontam.R`)

The QIIME2 pipeline (`pipeline/01_qiime2_asv_pipeline.sh`) exports these
files into `data/Phyloseq/`:

```
data/Phyloseq/
├── table-with-taxonomy.biom   # ASV x sample count table + taxonomy metadata
├── tree.nwk                   # unrooted ASV phylogenetic tree
└── dna-sequences.fasta        # representative sequences per ASV
```

`R/01_import_decontam.R` reads these three files by name at the top of the
script (`biom_file`, `tree_file`, `fasta_file`) along with
`data/metadata_pbr2.tsv`. If your own export uses different filenames or a
different folder layout, only those four lines need to change — the rest
of the script does not reference filenames directly.

## Functional annotation (optional)

- FAPROTAX database and `collapse_table.py` script:
  <http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Download>

## Intermediate `.rds` objects

Scripts `R/01`–`R/05` save intermediate objects (`phylo_raw.rds`,
`XData_final.rds`, fitted `Hmsc` objects, posterior prediction arrays,
variance-partitioning objects, etc.) into `results/`. These are not
version-controlled because the posterior prediction arrays alone can be
several GB. Regenerate them by running the R scripts in order, or contact
the authors for the archived copies used in the manuscript.
