# RNA-Seq Analysis Pipeline
### Salmon → DESeq2 → GSEA | Per-Donor Paired Analysis | GRCh38 / Ensembl 115

[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-blue)]()
[![Salmon](https://img.shields.io/badge/Salmon-1.10%2B-green)]()
[![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-red)]()
[![GSEA](https://img.shields.io/badge/GSEA-PreRanked-orange)]()

---

## Overview

Bulk RNA-seq pipeline for paired-end sequencing data. Produces per-donor differential expression and pathway enrichment results across any number of treatment vs. control comparisons.

```
FASTQ (R1 / R2)
      │
      ▼
Salmon quant        ← transcript-level quantification, GC/seq bias correction
      │
      ▼
tximeta → DESeq2    ← transcript → gene aggregation, per-donor DE analysis
      │
      ▼
GSEA PreRanked      ← Wald stat ranking, MSigDB Hallmark gene sets
```

**Key design:** donors are never pooled as replicates. Each donor runs independently — 4 donors × 3 comparisons = 12 DESeq2 outputs and 12 GSEA runs.

---

## Quickstart

```bash
# 1. Clone
git clone https://github.com/Capoinde/DEseq-2-and-GSEA-analysis.git
cd DEseq-2-and-GSEA-analysis

# 2. Install tools (see Platform Setup below)

# 3. Configure
cp config/config.template.sh config/config.sh
# Edit config/config.sh with your paths and comparisons

# 4. Run
source config/config.sh
bash pipeline/rnaseq_pipeline_salmon.sh 2>&1 | tee results/pipeline.log
```

---

## Platform Setup

> **Full step-by-step instructions with error fixes and troubleshooting are in:**
> `docs/RNAseq_MultiPlatform_Setup_Guide.docx`

### Windows

Windows requires WSL (Windows Subsystem for Linux) since bioinformatics tools are built for Linux.

```powershell
# Run in PowerShell as Administrator
wsl --install
```

After restart, open Ubuntu from the Start menu. All subsequent commands run inside the WSL terminal in VS Code (connect via the green >< button, bottom-left).

**Critical Windows note:** FASTQ files must be copied into the Linux filesystem (`~/fastq_data/`) before running. Files accessed through `/mnt/c/` cause read errors.

```bash
# Copy FASTQs from Windows to Linux filesystem
mkdir -p ~/fastq_data
find "/mnt/c/Users/YourName/path/to/fastq" \
     -mindepth 2 -name "*.fastq" -not -name "*Zone*" | while read f; do
    base=$(basename "$f")
    sid=$(echo "$base" | sed 's/_S[0-9]*_L[0-9]*_R\([12]\)_001\.fastq//')
    read_num=$(echo "$base" | grep -o '_R[12]_' | tr -d '_')
    sid_hyphen=$(echo "$sid" | sed 's/_/-/g')
    sudo cp "$f" ~/fastq_data/"${sid_hyphen}_${read_num}.fastq"
    sudo chown $USER:$USER ~/fastq_data/"${sid_hyphen}_${read_num}.fastq"
done
```

### macOS

Terminal is natively Unix-based — no WSL needed.

```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Java
brew install openjdk@17
echo 'export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

Point `FASTQ_DIR` directly at your FASTQ location — no copying required on macOS.

### Linux (Ubuntu/Debian)

Native environment — most straightforward setup.

```bash
sudo apt-get update
sudo apt-get install -y git default-jdk python3 python3-pip wget unzip
```

---

## Installing Tools (All Platforms)

### Miniconda + Salmon

```bash
# Install Miniconda (use MacOSX-arm64 version on Apple Silicon)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda3
$HOME/miniconda3/bin/conda init bash   # use 'conda init zsh' on macOS
source ~/.bashrc

# Accept Terms of Service
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install Salmon
conda install -y -c bioconda -c conda-forge salmon
salmon --version
```

### R + DESeq2

```bash
# Write to file to avoid bash special character issues with the ! operator
cat > /tmp/install_r_packages.R << 'EOF'
options(repos = c(CRAN = 'https://cloud.r-project.org'))
lib_path <- file.path(Sys.getenv('HOME'), 'R', 'library')
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', lib = lib_path)
BiocManager::install(
    c('DESeq2','tximeta','tximport','GenomicFeatures'),
    lib = lib_path, ask = FALSE, update = FALSE)
install.packages(c('jsonlite','ggplot2','ggrepel','pheatmap'), lib = lib_path)
message('R packages done')
EOF
Rscript /tmp/install_r_packages.R
echo 'export R_LIBS_USER="$HOME/R/library"' >> ~/.bashrc && source ~/.bashrc
```

### GSEA + Reference Files

```bash
# GSEA jar
mkdir -p ~/tools/gsea
wget 'https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.3/GSEA_Linux_4.3.3.zip' \
     -O /tmp/gsea.zip
unzip -j /tmp/gsea.zip '*/gsea-cli.jar' -d ~/tools/gsea/

# Reference files
mkdir -p ~/ref
wget 'https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' \
     -O ~/ref/cdna.all.fa.gz
wget 'https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz' \
     -O ~/ref/genome.gtf.gz && gunzip ~/ref/genome.gtf.gz
wget 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt' \
     -O ~/ref/h.all.v2023.2.Hs.symbols.gmt

# Build Salmon index (one-time, ~10 minutes)
salmon index --transcripts ~/ref/cdna.all.fa.gz \
             --index ~/ref/salmon_index --threads 4
```

---

## FASTQ Naming Convention

The pipeline matches filenames to `Sample_ID` values in your samplesheet.

```
# Required format
<Sample_ID>_R1.fastq         uncompressed paired-end R1
<Sample_ID>_R1.fastq.gz      gzip compressed
<Sample_ID>_R2.fastq         paired-end R2

# Absci internal format — auto-detected by organize_fastqs.py
<Sample_ID>_S<N>_L<N>_R1_001.fastq
Example: seqRNA49_HF_CUT_26_010_6hr_IgG_S1_L001_R1_001.fastq
```

The `organize_fastqs.py` script handles hyphen-to-underscore conversion automatically between samplesheet Sample_IDs and actual filenames.

---

## Samplesheet Format

| Column | Description |
|--------|-------------|
| `Sample_ID` | Must match FASTQ filename prefix exactly. No spaces. |
| `Treatment` | Must match values used in `COMPARISONS_STR` exactly. |
| `Condition` | Use `Control` for baseline; `Group1`, `Group2` etc for treatments. |
| `Donor_ID` | Individual identifier — each donor runs independently. |

---

## Configuration

```bash
cp config/config.template.sh config/config.sh
```

Key settings in `config/config.sh`:

```bash
export FASTQ_DIR="/home/username/fastq_data"      # Windows/Linux
# export FASTQ_DIR="/Users/username/fastq_data"  # macOS

export OUTDIR="$(pwd)/results"
export SAMPLESHEET="$(pwd)/samplesheets/my_project.csv"
export TRANSCRIPTOME_FA="$HOME/ref/cdna.all.fa.gz"
export GTF="$HOME/ref/genome.gtf"
export SALMON_INDEX="$HOME/ref/salmon_index"
export HALLMARK_GMT="$HOME/ref/h.all.v2023.2.Hs.symbols.gmt"
export GSEA_JAR="$HOME/tools/gsea/gsea-cli.jar"
export THREADS="8"

# Use COMPARISONS_STR (plain string) — NOT a bash array
export COMPARISONS_STR="ABS-201,IgG,ABS201_vs_IgG|||IgG+PRL,IgG,IgGPRL_vs_IgG|||ABS201+PRL,IgG+PRL,ABS201PRL_vs_IgGPRL"
```

> **Always use `COMPARISONS_STR`** with `|||` separators, not a bash array. Bash arrays cannot be exported between scripts and cause `COMPARISONS is empty` errors.

---

## Output Structure

```
results/
├── pipeline.log
├── salmon/
│   └── <Sample_ID>/quant.sf
├── deseq2/
│   ├── DESeq2_summary.csv
│   └── donor_<N>/
│       ├── ABS201_vs_IgG_donor<N>_deseq2_out.csv
│       ├── IgGPRL_vs_IgG_donor<N>_deseq2_out.csv
│       └── ABS201PRL_vs_IgGPRL_donor<N>_deseq2_out.csv
└── gsea/
    └── donor_<N>/
        └── <CompLabel>/
            └── index.html    ← open in browser
```

---

## Common Errors

| Error | Fix |
|-------|-----|
| `chmod: Operation not permitted` | Use `bash script.sh` not `./script.sh` |
| `COMPARISONS is empty` | Use `COMPARISONS_STR` string, not a bash array |
| `Error reading FASTA/Q stream` | Copy FASTQs to `~/fastq_data/` not `/mnt/c/` |
| `Is a directory` on FASTQ | Use `-mindepth 2` in find — files nested in same-named folders |
| `locale::facet error` | Run `sudo locale-gen en_US.UTF-8` then rebuild Salmon index |
| Salmon stalled for hours | Normal — 100 bootstraps = 1-2 hrs/sample. Reduce to 10 for testing |
| `conda: command not found` | Run `source ~/.bashrc` after Miniconda install |
| `CondaToSNonInteractiveError` | Run `conda tos accept` commands before installing packages |

---

## Repository Structure

```
DEseq-2-and-GSEA-analysis/
├── README.md
├── docs/
│   └── RNAseq_MultiPlatform_Setup_Guide.docx
├── config/
│   ├── config.template.sh
│   └── config.sh                  ← gitignored
├── pipeline/
│   └── rnaseq_pipeline_salmon.sh
├── r_scripts/
│   └── run_deseq2_salmon.R
├── python_scripts/
│   └── make_rnk.py
├── scripts/
│   └── organize_fastqs.py
├── samplesheets/
│   ├── samplesheet_template.csv
│   └── diversity_samplesheet.csv
└── results/                       ← gitignored
```

---

## Reproducibility

Record the pipeline version used for published analyses:

```bash
git rev-parse HEAD
```

Include this commit hash alongside tool versions (logged to `results/pipeline.log`) in your methods section.

---

## License

MIT
