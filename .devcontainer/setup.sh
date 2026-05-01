#!/usr/bin/env bash
# =============================================================================
# .devcontainer/setup.sh
#
# Runs once when a GitHub Codespace is created (postCreateCommand).
# Safe to re-run — all steps check if already complete before executing.
#
# Installs:
#   1. System packages (wget, curl, build tools, R)
#   2. Salmon (transcript quantification)
#   3. Python packages (pandas, numpy, matplotlib, etc.)
#   4. R packages (DESeq2, tximeta, tximport, ggplot2, etc.)
#   5. GSEA CLI jar
#   6. Reference files (Ensembl cDNA FASTA + GTF + MSigDB GMT)
#   7. Salmon index (built from downloaded FASTA)
#   8. Project config (writes ~/.rnaseq_env for all terminal sessions)
#
# To use a different genome, edit the SPECIES and ENSEMBL_RELEASE variables
# below, or override them as environment variables before running:
#   SPECIES=mouse ENSEMBL_RELEASE=112 bash .devcontainer/setup.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration — edit these for your organism / release
# ---------------------------------------------------------------------------
SPECIES="${SPECIES:-human}"                  # human | mouse | zebrafish | etc.
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-115}"
SALMON_VERSION="${SALMON_VERSION:-1.10.3}"
GSEA_VERSION="${GSEA_VERSION:-4.3.3}"

REF_DIR="${REF_DIR:-${HOME}/ref}"
TOOLS_DIR="${TOOLS_DIR:-${HOME}/tools}"
R_LIB="${HOME}/.R/library"

# ---------------------------------------------------------------------------
# Resolve Ensembl URLs for the selected species
# ---------------------------------------------------------------------------
case "${SPECIES}" in
  human)
    LATIN="homo_sapiens"
    ASSEMBLY="GRCh38"
    CDNA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${LATIN}/cdna/${LATIN^}.${ASSEMBLY}.cdna.all.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${LATIN}/${LATIN^}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf.gz"
    GMT_FILE="h.all.v2023.2.Hs.symbols.gmt"
    ;;
  mouse)
    LATIN="mus_musculus"
    ASSEMBLY="GRCm39"
    CDNA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${LATIN}/cdna/${LATIN^}.${ASSEMBLY}.cdna.all.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${LATIN}/${LATIN^}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf.gz"
    GMT_FILE="h.all.v2023.2.Mm.symbols.gmt"
    ;;
  *)
    echo "  Species '${SPECIES}' not pre-configured."
    echo "  Set CDNA_URL and GTF_URL manually and re-run."
    CDNA_URL="${CDNA_URL:-}"
    GTF_URL="${GTF_URL:-}"
    GMT_FILE="h.all.v2023.2.Hs.symbols.gmt"
    ;;
esac

GMT_URL="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/${GMT_FILE}"
CDNA_FA="${REF_DIR}/cdna.all.fa.gz"
GTF_FILE="${REF_DIR}/genome.gtf"
SALMON_INDEX="${REF_DIR}/salmon_index"
GSEA_JAR="${TOOLS_DIR}/gsea/gsea-cli.jar"
HALLMARK_GMT="${REF_DIR}/${GMT_FILE}"

echo "============================================================"
echo " RNA-Seq Codespace Setup"
echo " Species:  ${SPECIES}  (${LATIN:-custom})"
echo " Ensembl:  release ${ENSEMBL_RELEASE}"
echo " Salmon:   ${SALMON_VERSION}"
echo "============================================================"

mkdir -p "${REF_DIR}" "${TOOLS_DIR}/gsea" "${R_LIB}"

# ---------------------------------------------------------------------------
# 1. System packages
# ---------------------------------------------------------------------------
echo "[1/7] System packages..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
    wget curl unzip git \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    r-base \
    r-base-dev
echo "  Done."

# ---------------------------------------------------------------------------
# 2. Salmon
# ---------------------------------------------------------------------------
echo "[2/7] Salmon ${SALMON_VERSION}..."
if command -v salmon &>/dev/null; then
    echo "  Already installed: $(salmon --version 2>&1)"
else
    wget -q \
      "https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz" \
      -O /tmp/salmon.tar.gz
    sudo mkdir -p /usr/local/salmon
    sudo tar -xzf /tmp/salmon.tar.gz -C /usr/local/salmon --strip-components=1
    sudo ln -sf /usr/local/salmon/bin/salmon /usr/local/bin/salmon
    echo "  Installed: $(salmon --version 2>&1)"
fi

# ---------------------------------------------------------------------------
# 3. Python packages
# ---------------------------------------------------------------------------
echo "[3/7] Python packages..."
pip install --quiet --upgrade pip
pip install --quiet pandas numpy scipy matplotlib seaborn
echo "  Done."

# ---------------------------------------------------------------------------
# 4. R / Bioconductor packages
# ---------------------------------------------------------------------------
echo "[4/7] R packages (first run takes ~10-15 min)..."

Rscript - <<REOF
lib <- "${R_LIB}"
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager",
                     repos="https://cloud.r-project.org",
                     lib=lib, quiet=TRUE)
}

bioc_pkgs <- c("DESeq2", "tximeta", "tximport",
               "GenomicFeatures", "AnnotationHub", "ensembldb")

cran_pkgs <- c("jsonlite", "ggplot2", "ggrepel",
               "pheatmap", "RColorBrewer", "dplyr")

installed <- rownames(installed.packages(lib.loc=lib))

to_bioc <- bioc_pkgs[!bioc_pkgs %in% installed]
to_cran <- cran_pkgs[!cran_pkgs %in% installed]

if (length(to_bioc) > 0) {
    message("Installing Bioconductor: ", paste(to_bioc, collapse=", "))
    BiocManager::install(to_bioc, lib=lib, ask=FALSE, update=FALSE)
}
if (length(to_cran) > 0) {
    message("Installing CRAN: ", paste(to_cran, collapse=", "))
    install.packages(to_cran, repos="https://cloud.r-project.org", lib=lib, quiet=TRUE)
}
message("R packages ready.")
REOF

# ---------------------------------------------------------------------------
# 5. GSEA CLI jar
# ---------------------------------------------------------------------------
echo "[5/7] GSEA CLI..."
if [ -f "${GSEA_JAR}" ]; then
    echo "  Already present: ${GSEA_JAR}"
else
    wget -q \
      "https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/${GSEA_VERSION%.*}/GSEA_Linux_${GSEA_VERSION}.zip" \
      -O /tmp/gsea.zip
    unzip -q -j /tmp/gsea.zip '*/gsea-cli.jar' -d "${TOOLS_DIR}/gsea"
    echo "  Installed: ${GSEA_JAR}"
fi

# ---------------------------------------------------------------------------
# 6. Reference files
# ---------------------------------------------------------------------------
echo "[6/7] Reference files..."

if [ ! -f "${CDNA_FA}" ] && [ -n "${CDNA_URL:-}" ]; then
    echo "  Downloading cDNA FASTA (~300-600 MB)..."
    wget -q "${CDNA_URL}" -O "${CDNA_FA}"
else
    echo "  cDNA FASTA: $([ -f "${CDNA_FA}" ] && echo 'present' || echo 'SKIPPED — CDNA_URL not set')"
fi

if [ ! -f "${GTF_FILE}" ] && [ -n "${GTF_URL:-}" ]; then
    echo "  Downloading GTF (~50 MB)..."
    wget -q "${GTF_URL}" -O "${GTF_FILE}.gz"
    gunzip "${GTF_FILE}.gz"
else
    echo "  GTF: $([ -f "${GTF_FILE}" ] && echo 'present' || echo 'SKIPPED — GTF_URL not set')"
fi

if [ ! -f "${HALLMARK_GMT}" ]; then
    echo "  Downloading Hallmark GMT..."
    wget -q "${GMT_URL}" -O "${HALLMARK_GMT}"
else
    echo "  Hallmark GMT: present"
fi

# ---------------------------------------------------------------------------
# 7. Build Salmon index
# ---------------------------------------------------------------------------
echo "[7/7] Salmon index..."
if [ -d "${SALMON_INDEX}" ] && [ -f "${SALMON_INDEX}/info.json" ]; then
    echo "  Already built: ${SALMON_INDEX}"
elif [ -f "${CDNA_FA}" ]; then
    echo "  Building from ${CDNA_FA} (~10 min)..."
    salmon index \
        --transcripts "${CDNA_FA}" \
        --index       "${SALMON_INDEX}" \
        --threads     4 \
        2>&1 | tail -5
    echo "  Index built: ${SALMON_INDEX}"
else
    echo "  SKIPPED — cDNA FASTA not present. Run after downloading:"
    echo "    salmon index --transcripts ${CDNA_FA} --index ${SALMON_INDEX} --threads 4"
fi

# ---------------------------------------------------------------------------
# Write environment file — sourced by every terminal session
# ---------------------------------------------------------------------------
ENV_FILE="${HOME}/.rnaseq_env"
cat > "${ENV_FILE}" <<EOF
# Auto-generated by .devcontainer/setup.sh — do not edit manually
# Re-run setup.sh to regenerate.

export THREADS="4"
export REF_DIR="${REF_DIR}"
export OUTDIR="\${OUTDIR:-/workspaces/rnaseq-pipeline/results}"
export SAMPLESHEET="\${SAMPLESHEET:-/workspaces/rnaseq-pipeline/samplesheets/your_project.csv}"
export TRANSCRIPTOME_FA="${CDNA_FA}"
export GTF="${GTF_FILE}"
export SALMON_INDEX="${SALMON_INDEX}"
export GSEA_JAR="${GSEA_JAR}"
export HALLMARK_GMT="${HALLMARK_GMT}"
export R_LIBS_USER="${R_LIB}"
EOF

# Source automatically in every new terminal
if ! grep -q "rnaseq_env" "${HOME}/.bashrc" 2>/dev/null; then
    echo "" >> "${HOME}/.bashrc"
    echo "# RNA-Seq pipeline environment" >> "${HOME}/.bashrc"
    echo "[ -f ${ENV_FILE} ] && source ${ENV_FILE}" >> "${HOME}/.bashrc"
fi

echo ""
echo "============================================================"
echo " Setup complete!"
echo ""
echo " Environment variables written to: ${ENV_FILE}"
echo " Sourced automatically in every new terminal."
echo ""
echo " Next steps:"
echo "   1. Copy config/config.template.sh → config/config.sh"
echo "   2. Set FASTQ_DIR and SAMPLESHEET in config/config.sh"
echo "   3. Add your comparisons to COMPARISONS= in config/config.sh"
echo "   4. Run: bash pipeline/rnaseq_pipeline.sh"
echo ""
echo " NOTE: FASTQs must not be stored in this repository."
echo "       Point FASTQ_DIR at an external location or mount."
echo "============================================================"
