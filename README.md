# Bulk RNA-Seq Pipeline
### FASTQ → Salmon → DESeq2 → GSEA
**Transcript-level quantification | Per-group differential expression | Hallmark pathway enrichment**

---

## Overview

A portable, reproducible bulk RNA-seq analysis pipeline for paired-end sequencing data. Designed to run on a local machine, HPC, cloud instance, or GitHub Codespaces with minimal configuration.

```
FASTQ (R1 / R2)
      │
      ▼
Salmon quant          ← quasi-mapping, GC/seq bias correction, bootstrapped
      │
      ▼
tximeta → DESeq2      ← transcript → gene aggregation, differential expression
      │
      ▼
GSEA PreRanked        ← Wald stat ranking, MSigDB gene sets
```

**Key design principles:**
- **Samplesheet-driven** — all sample metadata, groupings, and comparisons are defined in one CSV
- **Idempotent** — safely re-runnable; completed steps are detected and skipped automatically
- **Reference-agnostic** — works with any Ensembl genome release and any organism
- **Comparison-flexible** — supports any number of treatment vs. control pairings
- **Per-donor / per-replicate aware** — donors are never silently pooled; each is run independently

---

## Quickstart

```bash
# 1. Clone
git clone https://github.com/YOUR_USERNAME/rnaseq-pipeline.git
cd rnaseq-pipeline

# 2. Configure (edit paths in one place)
cp config/config.template.sh config/config.sh
nano config/config.sh

# 3. Download reference files (one-time, ~10 min)
bash scripts/download_references.sh --species human --release 115

# 4. Run
bash pipeline/rnaseq_pipeline.sh 2>&1 | tee results/pipeline.log
```

---

## Repository Structure

```
rnaseq-pipeline/
├── README.md
├── .devcontainer/
│   ├── devcontainer.json          ← GitHub Codespaces environment definition
│   └── setup.sh                   ← auto-installs all tools + reference files
├── config/
│   ├── config.template.sh         ← copy to config.sh and fill in your paths
│   └── config.sh                  ← gitignored — your local/project settings
├── pipeline/
│   └── rnaseq_pipeline.sh         ← main pipeline entry point
├── r_scripts/
│   └── run_deseq2.R               ← DESeq2 analysis (called by pipeline)
├── python_scripts/
│   └── make_rnk.py                ← DESeq2 output → GSEA ranked list converter
├── scripts/
│   ├── download_references.sh     ← one-time reference file downloader
│   └── sb_submit.py               ← Seven Bridges batch job submission template
├── samplesheets/
│   ├── samplesheet_template.csv   ← blank template with column descriptions
│   └── example_project.csv        ← working example with dummy data
├── docker/
│   └── Dockerfile                 ← for Seven Bridges / containerized runs
├── tests/
│   └── data/                      ← small example FASTQs for CI testing
├── results/                       ← gitignored; all outputs written here at runtime
└── .gitignore
```

---

## Samplesheet Format

All sample metadata and comparison groupings are defined in a single CSV.
Copy `samplesheets/samplesheet_template.csv` and fill in your values.

### Required columns

| Column | Description | Example |
|--------|-------------|---------|
| `Sample_ID` | Unique identifier — must match your FASTQ filename prefix exactly | `ctrl_donor1` |
| `Treatment` | What was applied to this sample | `DrugA`, `Vehicle`, `siKD` |
| `Condition` | Role in the experiment — `Control` or a group label | `Control`, `Group1` |
| `Donor_ID` | Biological replicate or individual identifier | `D1`, `Mouse3` |

### Optional columns (used if present)

| Column | Description |
|--------|-------------|
| `Timepoint` | For time-course experiments |
| `Batch` | For batch correction in downstream analysis |
| `Notes` | Free text — ignored by the pipeline |

### Example

```csv
Sample_ID,Treatment,Condition,Donor_ID
ctrl_D1,Vehicle,Control,D1
drugA_D1,DrugA,Group1,D1
drugB_D1,DrugB,Group2,D1
ctrl_D2,Vehicle,Control,D2
drugA_D2,DrugA,Group1,D2
drugB_D2,DrugB,Group2,D2
ctrl_D3,Vehicle,Control,D3
drugA_D3,DrugA,Group1,D3
drugB_D3,DrugB,Group2,D3
```

### Defining comparisons

Comparisons are specified in `config/config.sh` as a plain list.
The pipeline runs every comparison for every unique `Donor_ID` in the samplesheet.

```bash
# Format: "TreatmentLabel,ControlLabel,OutputLabel"
export COMPARISONS=(
  "DrugA,Vehicle,DrugA_vs_Vehicle"
  "DrugB,Vehicle,DrugB_vs_Vehicle"
  "DrugA,DrugB,DrugA_vs_DrugB"
)
```

With the example samplesheet above and 3 comparisons across 3 donors,
the pipeline produces **9 DESeq2 outputs and 9 GSEA runs**.

---

## Configuration

```bash
# config/config.sh — the only file you need to edit

# Directory containing your FASTQ files
# Files must be named:  <Sample_ID>_R1.fastq.gz  and  <Sample_ID>_R2.fastq.gz
export FASTQ_DIR="/path/to/your/fastq"

# Where all pipeline outputs will be written
export OUTDIR="./results"

# Your completed samplesheet
export SAMPLESHEET="./samplesheets/your_project.csv"

# Reference files — populated by scripts/download_references.sh
export TRANSCRIPTOME_FA="/ref/cdna.all.fa.gz"
export GTF="/ref/genome.gtf"
export SALMON_INDEX="/ref/salmon_index"

# Gene set database for GSEA — any .gmt from MSigDB works
export HALLMARK_GMT="/ref/h.all.v2023.2.Hs.symbols.gmt"

# GSEA CLI jar (downloaded by setup.sh or manually)
export GSEA_JAR="/tools/gsea/gsea-cli.jar"

# Comparisons: "Treatment,Control,OutputLabel"
export COMPARISONS=(
  "DrugA,Vehicle,DrugA_vs_Vehicle"
  "DrugB,Vehicle,DrugB_vs_Vehicle"
)

# Set to "false" for single-end data
export PAIRED_END="true"

# CPU threads
export THREADS="16"
```

---

## FASTQ Naming Convention

```
<Sample_ID>_R1.fastq.gz
<Sample_ID>_R2.fastq.gz
```

`Sample_ID` must match the `Sample_ID` column in your samplesheet exactly.
Avoid spaces and special characters — use underscores or hyphens.

**Single-end data:** set `PAIRED_END="false"` in `config.sh`.
The pipeline adjusts Salmon flags automatically.

---

## Downloading Reference Files

```bash
# Human GRCh38 / Ensembl 115
bash scripts/download_references.sh --species human --release 115

# Mouse GRCm39 / Ensembl 112
bash scripts/download_references.sh --species mouse --release 112

# Any organism — provide Ensembl FTP URLs directly
bash scripts/download_references.sh \
  --cdna-url "https://ftp.ensembl.org/.../your_organism.cdna.all.fa.gz" \
  --gtf-url  "https://ftp.ensembl.org/.../your_organism.gtf.gz"
```

MSigDB gene sets by organism (download from https://www.gsea-msigdb.org/gsea/msigdb/):

| Organism | Recommended GMT |
|----------|----------------|
| Human | `h.all.v2023.2.Hs.symbols.gmt` |
| Mouse | `h.all.v2023.2.Mm.symbols.gmt` |
| Rat | Use human orthologs or C2/C5 collections |

---

## Output Structure

```
results/
├── pipeline.log
├── salmon/
│   └── <Sample_ID>/
│       ├── quant.sf               ← transcript-level quantification (TPM + counts)
│       ├── quant.genes.sf         ← gene-level summary
│       └── aux_info/              ← bias models, mapping rate stats
├── deseq2/
│   ├── DESeq2_summary.csv         ← sig-gene counts across all comparisons/donors
│   └── <Donor_ID>/
│       └── <CompLabel>_<DonorID>_deseq2_out.csv
└── gsea/
    └── <Donor_ID>/
        └── <CompLabel>/
            ├── <CompLabel>_<DonorID>.rnk
            └── <CompLabel>_<DonorID>.GseaPreranked.<timestamp>/
                ├── gsea_report_for_na_pos_*.html
                ├── gsea_report_for_na_neg_*.html
                └── index.html
```

---

## DESeq2 Output Columns

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID |
| `gene_name` | Gene symbol from Ensembl annotation via tximeta |
| `baseMean` | Mean normalized count across both samples |
| `log2FoldChange` | Treatment vs control (positive = upregulated in treatment) |
| `lfcSE` | Standard error of log2FoldChange |
| `stat` | Wald statistic — also used as GSEA ranking metric |
| `pvalue` | Nominal Wald test p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `significant` | `TRUE` if padj < 0.05 and \|log2FoldChange\| ≥ 1 |

---

## Method Notes

**Why Salmon instead of featureCounts or HTSeq?**
Salmon uses quasi-mapping and expectation-maximization to probabilistically assign reads across all compatible transcripts, preserving isoform-level information that read-counting approaches discard. GC and sequence-specific bias correction are enabled by default. Bootstrapped quantification uncertainty estimates are also generated and can be used for downstream sleuth or swish analyses.

**Why per-donor/replicate instead of pooled?**
When donors, animals, or cell lines are biologically distinct individuals, pooling them as replicates treats inter-individual variation as within-group noise and inflates false positives. Running each independently lets you assess the consistency of findings across individuals — a much stronger basis for biological conclusions. For experiments with true technical replicates of the same biological source, set `PER_DONOR="false"` in `config.sh` to pool normally.

**DESeq2 with n=2 (one treatment + one control per donor):**
The pipeline uses `fitType="local"` and `useT=TRUE` for robust dispersion estimation at small sample sizes, and disables `independentFiltering` which requires replicates to be meaningful. For experiments with ≥3 samples per group, these settings are relaxed automatically.

**GSEA ranking metric:** DESeq2 Wald statistic — signed (direction-aware), continuous, and handles ties better than fold change alone.

---

## Dependencies

| Tool | Minimum version | Purpose |
|------|----------------|---------|
| [Salmon](https://github.com/COMBINE-lab/salmon) | 1.10 | Transcript quantification |
| R | 4.2 | Statistical analysis |
| DESeq2 | 1.38 | Differential expression |
| tximeta | 1.16 | Transcript → gene aggregation |
| Java | 11 | GSEA CLI runtime |
| Python | 3.9 | Pipeline orchestration |
| pandas | latest | Data manipulation |

### Install R packages

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximeta", "tximport", "GenomicFeatures"))
install.packages(c("jsonlite", "ggplot2", "ggrepel", "pheatmap"))
```

All dependencies are installed automatically in GitHub Codespaces — see `.devcontainer/`.

---

## Running on Seven Bridges

See [`docker/Dockerfile`](docker/Dockerfile) for a container image with all dependencies pre-installed, suitable for use as a custom tool in the Seven Bridges Workflow Editor.

For batch job submission via the SB Python API, see [`scripts/sb_submit.py`](scripts/sb_submit.py).
The script reads your `config.sh` comparisons and samplesheet, then submits one task per `(Donor_ID, CompLabel)` pair programmatically — no manual clicking required.

---

## Reproducibility

When publishing or sharing results, record the exact pipeline version used:

```bash
git rev-parse HEAD
```

Include this commit hash alongside the tool versions (logged automatically to `results/pipeline.log`) in your methods section or electronic lab notebook.

---

## Contributing

Pull requests are welcome. Please open an issue first to discuss proposed changes to core pipeline logic. For bug reports, include the relevant section of `results/pipeline.log`.

---

## License

MIT — see [LICENSE](LICENSE)
