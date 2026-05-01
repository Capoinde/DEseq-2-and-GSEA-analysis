#!/usr/bin/env bash
# =============================================================================
# config/config.template.sh
#
# Copy this file to config/config.sh and fill in your values.
# config/config.sh is gitignored — your paths stay local.
#
# Usage:
#   source config/config.sh
#   bash pipeline/rnaseq_pipeline.sh
# =============================================================================

# ---------------------------------------------------------------------------
# Data paths  (required)
# ---------------------------------------------------------------------------

# Directory containing your FASTQ files.
# Files must be named: <Sample_ID>_R1.fastq.gz  and  <Sample_ID>_R2.fastq.gz
# Sample_ID must match the Sample_ID column in your samplesheet exactly.
export FASTQ_DIR="/path/to/your/fastq"

# Where all pipeline outputs will be written.
export OUTDIR="./results"

# Your completed samplesheet (see samplesheets/samplesheet_template.csv).
export SAMPLESHEET="./samplesheets/your_project.csv"

# ---------------------------------------------------------------------------
# Reference files  (populated by scripts/download_references.sh)
# ---------------------------------------------------------------------------

# Ensembl cDNA FASTA — used to build the Salmon index.
export TRANSCRIPTOME_FA="/ref/cdna.all.fa.gz"

# Ensembl GTF — used by tximeta to aggregate transcripts to genes.
export GTF="/ref/genome.gtf"

# Salmon transcriptome index (built from TRANSCRIPTOME_FA — one-time).
export SALMON_INDEX="/ref/salmon_index"

# ---------------------------------------------------------------------------
# GSEA
# ---------------------------------------------------------------------------

# Any .gmt file from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/).
# Common choices:
#   Hallmark:  h.all.v2023.2.Hs.symbols.gmt
#   KEGG:      c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt
#   GO BP:     c5.go.bp.v2023.2.Hs.symbols.gmt
export HALLMARK_GMT="/ref/h.all.v2023.2.Hs.symbols.gmt"

# GSEA CLI jar (downloaded by setup.sh or manually from gsea-msigdb.org).
export GSEA_JAR="/tools/gsea/gsea-cli.jar"

# ---------------------------------------------------------------------------
# Comparisons
# Format: "TreatmentLabel,ControlLabel,OutputLabel"
#
# The pipeline runs every comparison for every unique Donor_ID in the
# samplesheet. With 3 comparisons and 4 donors you get 12 output sets.
#
# TreatmentLabel and ControlLabel must match values in the Treatment column
# of your samplesheet exactly.
# ---------------------------------------------------------------------------
export COMPARISONS=(
  "DrugA,Vehicle,DrugA_vs_Vehicle"
  "DrugB,Vehicle,DrugB_vs_Vehicle"
  # "DrugA,DrugB,DrugA_vs_DrugB"    # uncomment to add direct comparisons
)

# ---------------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------------

# "true" for paired-end (R1 + R2); "false" for single-end (R1 only).
export PAIRED_END="true"

# DESeq2 significance thresholds.
export PADJ_THRESHOLD="0.05"
export LFC_THRESHOLD="1"

# GSEA minimum / maximum gene set size.
export GSEA_SET_MIN="15"
export GSEA_SET_MAX="500"
export GSEA_NPERM="1000"

# CPU threads for Salmon and other parallel steps.
export THREADS="16"
