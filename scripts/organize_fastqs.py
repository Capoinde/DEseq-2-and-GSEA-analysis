#!/usr/bin/env python3
"""
scripts/organize_fastqs.py

Prepares FASTQ files for the RNA-seq pipeline by creating symbolic links
(or copies) in the expected naming format:
    <Sample_ID>_R1.fastq.gz
    <Sample_ID>_R2.fastq.gz

Handles multiple naming conventions automatically, including:
  - Absci internal format:  seqRNA49_HF_CUT_26_010_6hr_IgG_S1_L001_R1_001.fastq
  - Compressed format:      seqRNA49-HF-CUT-26-010-6hr-IgG_R1.fastq.gz
  - Simple format:          seqRNA49_R1.fastq.gz
  - Azenta format:          seqRNA49_HF_CUT_26_010_6hr_IgG_R1_001.fastq.gz

Usage:
    # Dry run first — shows what would happen without creating anything
    python3 scripts/organize_fastqs.py \\
        --samplesheet samplesheets/diversity_samplesheet.csv \\
        --source-dir  "/home/cpoindexter/rnaseq-pipeline/RNA-seq FastQ and Data" \\
        --output-dir  ~/fastq \\
        --dry-run

    # Run for real
    python3 scripts/organize_fastqs.py \\
        --samplesheet samplesheets/diversity_samplesheet.csv \\
        --source-dir  "/home/cpoindexter/rnaseq-pipeline/RNA-seq FastQ and Data" \\
        --output-dir  ~/fastq \\
        --mode symlink
"""

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Naming patterns — ordered from most specific to most general.
# {sid} is replaced with the Sample_ID from the samplesheet.
# The sid is also tried with hyphens converted to underscores automatically.
# ---------------------------------------------------------------------------
R1_PATTERNS = [
    # Absci internal format: seqRNA49_HF_CUT_26_010_6hr_IgG_S1_L001_R1_001.fastq
    "{sid}_S*_L*_R1_001.fastq",
    "{sid}_S*_L*_R1_001.fastq.gz",
    # Absci internal without lane: seqRNA49_HF_CUT_26_010_6hr_IgG_R1_001.fastq
    "{sid}_R1_001.fastq",
    "{sid}_R1_001.fastq.gz",
    # Standard format
    "{sid}_R1.fastq.gz",
    "{sid}_R1.fastq",
    "{sid}_R1.fq.gz",
    "{sid}_R1.fq",
    # Numbered format
    "{sid}_1.fastq.gz",
    "{sid}_1.fastq",
    "{sid}_1.fq.gz",
    # Dot separator
    "{sid}.R1.fastq.gz",
    "{sid}.R1.fastq",
    "{sid}.1.fastq.gz",
    # Read spelled out
    "{sid}_read1.fastq.gz",
    "{sid}_read1.fastq",
]

R2_PATTERNS = [
    p.replace("_R1_", "_R2_")
     .replace("_R1.", "_R2.")
     .replace(".R1.", ".R2.")
     .replace("_1.", "_2.")
     .replace("_read1", "_read2")
    for p in R1_PATTERNS
]


def normalise_sid(sample_id: str) -> list[str]:
    """
    Return a list of sid variants to try when searching for files.
    Your samplesheet uses hyphens  (seqRNA49-HF-CUT-26-010-6hr-IgG)
    but Absci FASTQ filenames use underscores (seqRNA49_HF_CUT_26_010_6hr_IgG).
    We try both so neither format needs manual editing.
    """
    variants = [sample_id]
    # hyphens → underscores
    underscored = sample_id.replace("-", "_")
    if underscored != sample_id:
        variants.append(underscored)
    # underscores → hyphens
    hyphenated = sample_id.replace("_", "-")
    if hyphenated != sample_id:
        variants.append(hyphenated)
    return variants


def find_fastq(source_dir: Path, sample_id: str, patterns: list[str]) -> Path | None:
    """
    Search source_dir for a FASTQ matching sample_id.
    Tries every combination of sid variant × pattern.
    Falls back to a partial filename scan if nothing matches exactly.
    """
    sid_variants = normalise_sid(sample_id)

    for sid_variant in sid_variants:
        for pat in patterns:
            glob_pat = pat.replace("{sid}", sid_variant)
            matches = sorted(source_dir.glob(glob_pat))
            if matches:
                return matches[0]

    # -----------------------------------------------------------------
    # Last resort — scan every file in the directory and look for a
    # file whose name contains the sample number (e.g. "seqRNA49") and
    # the R1/R2 tag.  This catches unusual naming conventions.
    # -----------------------------------------------------------------
    # Extract just the sample number prefix e.g. "seqRNA49"
    sample_prefix = re.split(r'[-_]', sample_id)[0].lower()
    r1_tags = ["_r1_", "_r1.", ".r1.", "_read1", "_1."]

    for f in sorted(source_dir.iterdir()):
        name_lower = f.name.lower()
        if sample_prefix in name_lower:
            # Make sure the full sample ID core matches, not just the number
            # e.g. seqRNA49 should not match seqRNA490
            sid_core = sample_id.lower().replace("-", "_")
            # Try the underscore version of the ID
            if sid_core in name_lower or sid_core.replace("_","-") in name_lower:
                for tag in r1_tags:
                    if tag in name_lower:
                        return f

    return None


def output_extension(src: Path) -> str:
    """
    Decide what extension the output symlink/copy should have.
    The pipeline expects .fastq.gz — if the source is uncompressed
    we still name the link .fastq.gz (Salmon handles both transparently
    when --readFilesCommand is not set, or use --no-gz flag at quant time).
    """
    return ".fastq.gz"


def load_samplesheet(path: Path) -> list[dict]:
    """Load samplesheet CSV, skip comment lines starting with #."""
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(
            line for line in fh if not line.strip().startswith("#")
        )
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Organise FASTQs into pipeline-expected naming format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--samplesheet", required=True,
                        help="Path to your samplesheet CSV")
    parser.add_argument("--source-dir",  required=True,
                        help="Directory containing your original FASTQ files")
    parser.add_argument("--output-dir",  required=True,
                        help="Directory where renamed/linked files will be placed")
    parser.add_argument("--mode", default="symlink", choices=["symlink", "copy"],
                        help="symlink = fast, no extra disk space (default); "
                             "copy = duplicate files (use if symlinks cause issues)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would happen without creating any files")
    parser.add_argument("--paired-end", default="true", choices=["true", "false"],
                        help="true for R1+R2 paired-end (default); false for single-end R1 only")
    args = parser.parse_args()

    samplesheet = Path(args.samplesheet)
    source_dir  = Path(args.source_dir)
    output_dir  = Path(args.output_dir)
    paired      = args.paired_end.lower() == "true"
    dry_run     = args.dry_run

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not samplesheet.exists():
        sys.exit(f"ERROR: Samplesheet not found: {samplesheet}")
    if not source_dir.exists():
        sys.exit(f"ERROR: Source directory not found: {source_dir}")

    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    samples = load_samplesheet(samplesheet)
    if not samples:
        sys.exit("ERROR: No samples found in samplesheet.")

    sample_ids = [s["Sample_ID"] for s in samples if s.get("Sample_ID")]

    print(f"\nSamplesheet : {samplesheet}")
    print(f"Source dir  : {source_dir}")
    print(f"Output dir  : {output_dir}")
    print(f"Mode        : {'DRY RUN — no files will be created' if dry_run else args.mode}")
    print(f"Paired-end  : {paired}")
    print(f"Samples     : {len(sample_ids)}")
    print()
    print("-" * 70)

    found     = 0
    not_found = []

    for sid in sample_ids:

        r1_src = find_fastq(source_dir, sid, R1_PATTERNS)
        r2_src = find_fastq(source_dir, sid, R2_PATTERNS) if paired else None

        # Output files always end in .fastq.gz regardless of source extension
        r1_dst = output_dir / f"{sid}_R1.fastq.gz"
        r2_dst = output_dir / f"{sid}_R2.fastq.gz"

        # Report
        r1_status = f"FOUND  → {r1_src.name}" if r1_src else "NOT FOUND"
        r2_status = (f"FOUND  → {r2_src.name}" if r2_src else "NOT FOUND") if paired else "N/A"

        print(f"\n  {sid}")
        print(f"    R1: {r1_status}")
        if paired:
            print(f"    R2: {r2_status}")

        missing = []
        if not r1_src:
            missing.append("R1")
        if paired and not r2_src:
            missing.append("R2")

        if missing:
            print(f"    ⚠  Missing: {', '.join(missing)}")
            not_found.append(sid)
            continue

        # Create symlinks or copies
        if not dry_run:
            for src, dst in [(r1_src, r1_dst)] + ([(r2_src, r2_dst)] if paired else []):
                if dst.exists():
                    print(f"    [SKIP] Already exists: {dst.name}")
                    continue
                if args.mode == "symlink":
                    os.symlink(src.resolve(), dst)
                    print(f"    ✓  Linked : {dst.name}")
                else:
                    shutil.copy2(src, dst)
                    print(f"    ✓  Copied : {dst.name}")
        else:
            print(f"    [DRY RUN] Would create: {r1_dst.name}  ←  {r1_src.name}")
            if paired and r2_src:
                print(f"    [DRY RUN] Would create: {r2_dst.name}  ←  {r2_src.name}")

        found += 1

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    print("-" * 70)
    print(f"\nSummary:")
    print(f"  ✓  Processed : {found} / {len(sample_ids)} samples")

    if not_found:
        print(f"\n  ✗  Could not find FASTQs for {len(not_found)} sample(s):")
        for s in not_found:
            print(f"       {s}")
        print(f"\n  Tip: Check that the sample number (e.g. seqRNA49) appears in")
        print(f"       your FASTQ filenames inside: {source_dir}")
    else:
        print(f"\n  All FASTQs found and organised.")
        if not dry_run:
            print(f"\n  Next step — set this in config/config.sh:")
            print(f"    export FASTQ_DIR=\"{output_dir}\"")

    if dry_run:
        print(f"\n  This was a DRY RUN — no files were created.")
        print(f"  Remove --dry-run to apply changes.")


if __name__ == "__main__":
    main()
