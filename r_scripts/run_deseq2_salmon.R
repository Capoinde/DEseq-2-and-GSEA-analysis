#!/usr/bin/env Rscript
# =============================================================================
# r_scripts/run_deseq2_salmon.R
#
# Runs DESeq2 differential expression for ONE donor × ONE comparison.
# Called by the pipeline script once per (donor, comparison) pair.
#
# Usage:
#   Rscript r_scripts/run_deseq2_salmon.R <pair_meta.json> <output_dir>
#
# pair_meta.json must contain:
#   {
#     "treat_sample": "seqRNA50",          <- Sample_ID of the treatment sample
#     "ctrl_sample":  "seqRNA49",          <- Sample_ID of the control sample
#     "donor":        "1",                 <- Donor_ID label
#     "group_label":  "ABS201_vs_IgG",     <- comparison name used for output file
#     "salmon_dir":   "/home/user/results/salmon",   <- folder containing per-sample quant.sf dirs
#     "gtf_path":     "/home/user/ref/genome.gtf"    <- Ensembl GTF for gene annotation
#   }
#
# Output:
#   <output_dir>/<group_label>_donor<donor>_deseq2_out.csv
# =============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(tximeta)
    library(jsonlite)
    library(GenomicFeatures)
})

# ---------------------------------------------------------------------------
# 1. Read arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript run_deseq2_salmon.R <pair_meta.json> <output_dir>")
}

meta_file  <- args[1]
out_dir    <- args[2]

if (!file.exists(meta_file)) stop(paste("pair_meta.json not found:", meta_file))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- fromJSON(meta_file)

treat_s <- meta$treat_sample
ctrl_s  <- meta$ctrl_sample
donor   <- meta$donor
label   <- meta$group_label
salmon  <- meta$salmon_dir
gtf     <- meta$gtf_path

message(sprintf("=== DESeq2: %s | Donor %s ===", label, donor))
message(sprintf("  Treatment: %s", treat_s))
message(sprintf("  Control:   %s", ctrl_s))

# ---------------------------------------------------------------------------
# 2. Locate quant.sf files
# ---------------------------------------------------------------------------
treat_sf <- file.path(salmon, treat_s, "quant.sf")
ctrl_sf  <- file.path(salmon, ctrl_s,  "quant.sf")

for (f in c(treat_sf, ctrl_sf)) {
    if (!file.exists(f)) {
        stop(sprintf("quant.sf not found: %s\nCheck that Salmon ran successfully for this sample.", f))
    }
}

# ---------------------------------------------------------------------------
# 3. Build coldata for tximeta
# ---------------------------------------------------------------------------
coldata <- data.frame(
    names     = c(treat_s,     ctrl_s),
    files     = c(treat_sf,    ctrl_sf),
    condition = factor(
        c("treatment", "control"),
        levels = c("control", "treatment")   # control is the reference level
    ),
    donor     = donor,
    stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 4. Load quantification data with tximeta
#    tximeta auto-detects the Ensembl release from the Salmon index provenance.
#    If it cannot reach the internet, it falls back to a local linkedTxome.
# ---------------------------------------------------------------------------
message("  Loading quantification data via tximeta...")

se <- tryCatch({
    tximeta(coldata, type = "salmon", dropInfReps = FALSE)
}, error = function(e) {
    message("  tximeta auto-link failed; trying local linkedTxome registration...")

    salmon_index <- Sys.getenv("SALMON_INDEX")
    cdna_fa      <- Sys.getenv("TRANSCRIPTOME_FA")

    if (nchar(salmon_index) == 0 || nchar(cdna_fa) == 0) {
        stop(paste(
            "tximeta failed and SALMON_INDEX / TRANSCRIPTOME_FA env vars are not set.",
            "Set these before running, or ensure internet access for auto-linking.",
            sep = "\n"
        ))
    }

    makeLinkedTxome(
        indexDir  = salmon_index,
        source    = "Ensembl",
        organism  = "Homo sapiens",
        release   = "115",
        genome    = "GRCh38",
        fasta     = cdna_fa,
        gtf       = gtf,
        write     = FALSE
    )
    tximeta(coldata, type = "salmon", dropInfReps = FALSE)
})

# ---------------------------------------------------------------------------
# 5. Summarize transcript-level estimates to gene level
# ---------------------------------------------------------------------------
message("  Summarizing to gene level...")
gse <- summarizeToGene(se)

# ---------------------------------------------------------------------------
# 6. DESeq2
# ---------------------------------------------------------------------------
message("  Running DESeq2...")

dds <- DESeqDataSet(gse, design = ~ condition)

# Pre-filter: remove genes with fewer than 10 total counts across both samples
# This reduces multiple testing burden and speeds up the run
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message(sprintf("  Genes after pre-filter (>= 10 total counts): %d", nrow(dds)))

dds <- DESeq(
    dds,
    fitType              = "local",   # more robust than 'parametric' for n=2
    useT                 = TRUE,      # t-distribution tails for small sample sizes
    minmu                = 0.5        # minimum mean for numerical stability
)

res <- results(
    dds,
    contrast             = c("condition", "treatment", "control"),
    independentFiltering = FALSE,     # requires replicates; disabled for n=2
    alpha                = 0.05
)

# ---------------------------------------------------------------------------
# 7. Add gene symbols from tximeta row annotation
# ---------------------------------------------------------------------------
res_df         <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

rd <- as.data.frame(rowData(gse))
if ("gene_name" %in% colnames(rd)) {
    res_df$gene_name <- rd[res_df$gene_id, "gene_name"]
} else if ("symbol" %in% colnames(rd)) {
    res_df$gene_name <- rd[res_df$gene_id, "symbol"]
} else {
    res_df$gene_name <- res_df$gene_id
    message("  Note: gene_name not found in row annotation — using Ensembl IDs as names")
}

# ---------------------------------------------------------------------------
# 8. Flag significant genes and sort
# ---------------------------------------------------------------------------
padj_thresh <- as.numeric(Sys.getenv("PADJ_THRESHOLD", unset = "0.05"))
lfc_thresh  <- as.numeric(Sys.getenv("LFC_THRESHOLD",  unset = "1"))

res_df$significant <- (
    !is.na(res_df$padj) &
    res_df$padj < padj_thresh &
    abs(res_df$log2FoldChange) >= lfc_thresh
)

res_df <- res_df[order(res_df$pvalue, na.last = TRUE), ]

# ---------------------------------------------------------------------------
# 9. Write output
# ---------------------------------------------------------------------------
out_cols <- c(
    "gene_id", "gene_name",
    "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj",
    "significant"
)

out_file <- file.path(
    out_dir,
    sprintf("%s_donor%s_deseq2_out.csv", label, donor)
)

write.csv(res_df[, out_cols], out_file, row.names = FALSE)

# ---------------------------------------------------------------------------
# 10. Print summary
# ---------------------------------------------------------------------------
n_sig <- sum(res_df$significant, na.rm = TRUE)
n_up  <- sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE)
n_dn  <- sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)

message(sprintf("  Significant genes (padj < %.2f, |LFC| >= %.1f): %d",
                padj_thresh, lfc_thresh, n_sig))
message(sprintf("    Upregulated:   %d", n_up))
message(sprintf("    Downregulated: %d", n_dn))
message(sprintf("  Output: %s", out_file))
message("  Done.")
