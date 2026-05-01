#!/usr/bin/env python3
"""
python_scripts/make_rnk.py

Converts a DESeq2 output CSV into a GSEA PreRanked .rnk file.

The ranking metric used is the DESeq2 Wald statistic (column: stat).
This is preferred over fold change because it is:
  - Signed  (positive = upregulated in treatment)
  - Continuous  (no ties from rounding)
  - Already accounts for variance / uncertainty in the estimate

Usage:
    python3 python_scripts/make_rnk.py <deseq2_output.csv> <output.rnk>

Example:
    python3 python_scripts/make_rnk.py \\
        results/deseq2/donor_1/ABS201_vs_IgG_donor1_deseq2_out.csv \\
        results/gsea/donor_1/ABS201_vs_IgG/ABS201_vs_IgG_donor1.rnk

Output format (tab-separated, no header):
    GENE_SYMBOL    Wald_stat
    MYC            4.823
    TP53           3.201
    ...
    BRCA1          -2.441
"""

import sys
import os
import pandas as pd


def make_rnk(deseq_csv: str, rnk_out: str) -> None:
    """
    Read a DESeq2 CSV and write a sorted .rnk file for GSEA PreRanked.

    Parameters
    ----------
    deseq_csv : str
        Path to DESeq2 output CSV. Must contain columns: stat, gene_id.
        Optionally contains: gene_name (preferred — GSEA Hallmark uses symbols).
    rnk_out : str
        Path to write the .rnk file.
    """

    # -- Load DESeq2 results --------------------------------------------------
    if not os.path.exists(deseq_csv):
        sys.exit(f"ERROR: DESeq2 output not found: {deseq_csv}")

    df = pd.read_csv(deseq_csv)

    required_cols = {"stat", "gene_id"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"ERROR: Required columns missing from DESeq2 output: {missing}")

    # -- Drop rows with no Wald statistic (filtered genes) --------------------
    n_before = len(df)
    df = df.dropna(subset=["stat"])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(f"  Dropped {n_dropped} genes with NA stat (low-count filtered genes)")

    # -- Choose gene name column ----------------------------------------------
    # GSEA Hallmark gene sets use HGNC symbols (e.g. MYC, TP53).
    # If gene_name is available, use it. Fall back to gene_id (Ensembl IDs)
    # only if necessary — Ensembl IDs will NOT match Hallmark gene sets.
    if "gene_name" in df.columns:
        name_series = df["gene_name"].fillna(df["gene_id"])
        n_ensg = name_series.str.startswith("ENSG").sum()
        if n_ensg > 100:
            print(f"  WARNING: {n_ensg} entries still have Ensembl IDs as gene_name.")
            print("  GSEA Hallmark gene sets require HGNC symbols.")
            print("  Consider re-running tximeta with a full Ensembl annotation.")
    else:
        name_series = df["gene_id"]
        print("  WARNING: 'gene_name' column not found — using Ensembl IDs.")
        print("  GSEA Hallmark gene sets require HGNC symbols.")
        print("  Results may show very few or no enriched gene sets.")

    # -- Build ranked dataframe -----------------------------------------------
    df_rnk = pd.DataFrame({
        "gene": name_series,
        "stat": df["stat"]
    })

    # Remove duplicate gene names — keep the one with the highest absolute stat
    df_rnk["abs_stat"] = df_rnk["stat"].abs()
    df_rnk = (
        df_rnk
        .sort_values("abs_stat", ascending=False)
        .drop_duplicates(subset="gene", keep="first")
        .drop(columns="abs_stat")
        .sort_values("stat", ascending=False)
        .reset_index(drop=True)
    )

    # -- Write .rnk file ------------------------------------------------------
    os.makedirs(os.path.dirname(rnk_out) if os.path.dirname(rnk_out) else ".", exist_ok=True)
    df_rnk.to_csv(rnk_out, sep="\t", header=False, index=False)

    # -- Summary --------------------------------------------------------------
    n_pos = (df_rnk["stat"] > 0).sum()
    n_neg = (df_rnk["stat"] < 0).sum()
    print(f"  RNK file written: {rnk_out}")
    print(f"  Total genes ranked: {len(df_rnk)}")
    print(f"    Positive stat (upregulated in treatment): {n_pos}")
    print(f"    Negative stat (downregulated in treatment): {n_neg}")
    print(f"  Top 5 upregulated:   {list(df_rnk['gene'].head(5))}")
    print(f"  Top 5 downregulated: {list(df_rnk['gene'].tail(5))}")


# ---------------------------------------------------------------------------
# Entry point — accepts command-line arguments
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    deseq_csv = sys.argv[1]
    rnk_out   = sys.argv[2]

    print(f"\nConverting DESeq2 output → GSEA ranked list")
    print(f"  Input:  {deseq_csv}")
    print(f"  Output: {rnk_out}\n")

    make_rnk(deseq_csv, rnk_out)
