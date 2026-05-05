#!/usr/bin/env python3
"""
scripts/organize_fastqs.py

Prepares FASTQ files for the RNA-seq pipeline by creating symbolic links
(or copies) in the expected naming format:
    <Sample_ID>_R1.fastq.gz
    <Sample_ID>_R2.fastq.gz

Reads your samplesheet to know which Sample_IDs to expect, then searches
your source FASTQ directory for matching files using common naming patterns.

Usage:
    python3 scripts/organize_fastqs.py \\
        --samplesheet samplesheets/your_project.csv \\
        --source-dir  /mnt/c/Users/YourName/Downloads/fastq_files \\
        --output-dir  /home/yourname/fastq \\
        --mode        symlink        # or 'copy' if you want actual copies

Dry run (shows what would happen without doing anything):
    python3 scripts/organize_fastqs.py \\
        --samplesheet samplesheets/your_project.csv \\
        --source-dir  /mnt/c/Users/YourName/Downloads/fastq_files \\
        --output-dir  /home/yourname/fastq \\
        --dry-run
"""

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Naming patterns to try when searching for a FASTQ matching a Sample_ID.
# Add more patterns here if your sequencing core uses a different convention.
# ---------------------------------------------------------------------------
R1_PATTERNS = [
    "{sid}_R1.fastq.gz",
    "{sid}_R1.fq.gz",
    "{sid}_1.fastq.gz",
    "{sid}_1.fq.gz",
    "{sid}_read1.fastq.gz",
    "{sid}.R1.fastq.gz",
    "{sid}.1.fastq.gz",
    # Seven Bridges often appends run info — match loosely
    "{sid}*_R1*.fastq.gz",
    "{sid}*_R1*.fq.gz",
]
R2_PATTERNS = [p.replace("R1", "R2").replace("_1.", "_2.").replace("read1", "read2").replace(".1.", ".2.")
               for p in R1_PATTERNS]


def find_fastq(source_dir: Path, sample_id: str, patterns: list[str]) -> Path | None:
    """Search source_dir for a FASTQ matching sample_id using the pattern list."""
    for pat in patterns:
        # Exact match first
        candidate = source_dir / pat.format(sid=sample_id)
        if candidate.exists():
            return candidate
        # Glob match (for wildcard patterns)
        glob_pat = pat.format(sid=sample_id)
        if "*" in glob_pat:
            matches = sorted(source_dir.glob(glob_pat))
            if matches:
                return matches[0]
    # Last resort: case-insensitive partial scan
    sid_lower = sample_id.lower().replace("-", "_").replace(" ", "_")
    for f in sorted(source_dir.iterdir()):
        name_lower = f.name.lower().replace("-", "_").replace(" ", "_")
        if sid_lower in name_lower and f.suffix in (".gz",):
            for tag in ["_r1", "_1.", ".r1", "read1"]:
                if tag in name_lower:
                    return f
    return None


def load_samplesheet(path: Path) -> list[dict]:
    """Load samplesheet, skip comment lines."""
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(line for line in fh if not line.strip().startswith("#"))
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Organize FASTQs into pipeline-expected naming format."
    )
    parser.add_argument("--samplesheet",  required=True,  help="Path to your samplesheet CSV")
    parser.add_argument("--source-dir",   required=True,  help="Directory containing your original FASTQ files")
    parser.add_argument("--output-dir",   required=True,  help="Directory where renamed/linked files will be placed")
    parser.add_argument("--mode",         default="symlink", choices=["symlink", "copy"],
                        help="'symlink' creates links (fast, no extra disk space); 'copy' duplicates files")
    parser.add_argument("--dry-run",      action="store_true",
                        help="Show what would happen without actually doing anything")
    parser.add_argument("--paired-end",   default="true",  choices=["true", "false"],
                        help="'true' for R1+R2; 'false' for R1 only")
    args = parser.parse_args()

    samplesheet = Path(args.samplesheet)
    source_dir  = Path(args.source_dir)
    output_dir  = Path(args.output_dir)
    paired      = args.paired_end.lower() == "true"
    dry_run     = args.dry_run

    # Validate inputs
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
    print(f"\nSamplesheet:  {samplesheet}")
    print(f"Source dir:   {source_dir}")
    print(f"Output dir:   {output_dir}")
    print(f"Mode:         {'DRY RUN — no files will be created' if dry_run else args.mode}")
    print(f"Paired-end:   {paired}")
    print(f"Samples:      {len(sample_ids)}\n")
    print("-" * 70)

    found     = 0
    not_found = []

    for sid in sample_ids:
        r1_src = find_fastq(source_dir, sid, R1_PATTERNS)
        r2_src = find_fastq(source_dir, sid, R2_PATTERNS) if paired else None

        r1_dst = output_dir / f"{sid}_R1.fastq.gz"
        r2_dst = output_dir / f"{sid}_R2.fastq.gz"

        # Report what was found
        r1_status = f"FOUND  → {r1_src.name}" if r1_src else "NOT FOUND"
        r2_status = (f"FOUND  → {r2_src.name}" if r2_src else "NOT FOUND") if paired else "N/A (SE)"
        print(f"  {sid}")
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
                    print(f"    ✓  Linked:  {dst.name}")
                else:
                    shutil.copy2(src, dst)
                    print(f"    ✓  Copied:  {dst.name}")
        else:
            print(f"    [DRY RUN] Would create: {r1_dst.name}")
            if paired:
                print(f"    [DRY RUN] Would create: {r2_dst.name}")

        found += 1
        print()

    # Summary
    print("-" * 70)
    print(f"\nSummary:")
    print(f"  ✓  Processed: {found} / {len(sample_ids)} samples")

    if not_found:
        print(f"\n  ✗  Could not find FASTQs for {len(not_found)} sample(s):")
        for sid in not_found:
            print(f"       {sid}")
        print(f"\n  Check that these Sample_IDs appear (at least partially) in your")
        print(f"  FASTQ filenames, or rename the files manually to:")
        print(f"    <Sample_ID>_R1.fastq.gz")
        print(f"    <Sample_ID>_R2.fastq.gz")
    else:
        print(f"\n  All FASTQs found and organized.")
        print(f"  Set FASTQ_DIR={output_dir} in config/config.sh")

    if dry_run:
        print(f"\n  This was a DRY RUN — no files were created.")
        print(f"  Remove --dry-run to apply changes.")


if __name__ == "__main__":
    main()
