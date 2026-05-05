#!/usr/bin/env bash
# =============================================================================
# pipeline/rnaseq_pipeline_salmon.sh
#
# Main pipeline entry point.
# Calls r_scripts/run_deseq2_salmon.R and python_scripts/make_rnk.py
# as standalone files — do not move them without updating paths here.
#
# Steps:
#   0  Tool and file check
#   1  Salmon index (one-time build, skipped if already exists)
#   2  Parse samplesheet
#   3  Salmon quantification (one job per sample)
#   4  DESeq2 (one job per donor × comparison, calls run_deseq2_salmon.R)
#   5  GSEA PreRanked (one job per donor × comparison, calls make_rnk.py)
#   6  Summary table
#
# Usage:
#   source config/config.sh
#   bash pipeline/rnaseq_pipeline_salmon.sh 2>&1 | tee results/pipeline.log
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Locate repo root so all relative paths work from any working directory
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export REPO_ROOT

R_SCRIPT="${REPO_ROOT}/r_scripts/run_deseq2_salmon.R"
PY_RNK="${REPO_ROOT}/python_scripts/make_rnk.py"

# ---------------------------------------------------------------------------
# Load config if environment variables are not already set
# ---------------------------------------------------------------------------
if [ -z "${FASTQ_DIR:-}" ]; then
    CONFIG="${REPO_ROOT}/config/config.sh"
    if [ -f "${CONFIG}" ]; then
        # shellcheck source=/dev/null
        source "${CONFIG}"
        echo "  Loaded config: ${CONFIG}"
    else
        echo "ERROR: config/config.sh not found."
        echo "  Run:  cp config/config.template.sh config/config.sh"
        echo "  Then fill in your paths and comparisons."
        exit 1
    fi
fi

# Apply defaults for anything not set in config
FASTQ_DIR="${FASTQ_DIR:-${REPO_ROOT}/tests/data/fastq}"
OUTDIR="${OUTDIR:-${REPO_ROOT}/results}"
SALMON_INDEX="${SALMON_INDEX:-${HOME}/ref/salmon_index}"
GTF="${GTF:-${HOME}/ref/genome.gtf}"
TRANSCRIPTOME_FA="${TRANSCRIPTOME_FA:-${HOME}/ref/cdna.all.fa.gz}"
GSEA_JAR="${GSEA_JAR:-${HOME}/tools/gsea/gsea-cli.jar}"
HALLMARK_GMT="${HALLMARK_GMT:-${HOME}/ref/h.all.v2023.2.Hs.symbols.gmt}"
SAMPLESHEET="${SAMPLESHEET:-${REPO_ROOT}/samplesheets/diversity_samplesheet.csv}"
THREADS="${THREADS:-8}"
PAIRED_END="${PAIRED_END:-true}"
PADJ_THRESHOLD="${PADJ_THRESHOLD:-0.05}"
LFC_THRESHOLD="${LFC_THRESHOLD:-1}"
GSEA_SET_MIN="${GSEA_SET_MIN:-15}"
GSEA_SET_MAX="${GSEA_SET_MAX:-500}"
GSEA_NPERM="${GSEA_NPERM:-1000}"

export FASTQ_DIR OUTDIR SALMON_INDEX GTF TRANSCRIPTOME_FA
export GSEA_JAR HALLMARK_GMT SAMPLESHEET THREADS PAIRED_END
export PADJ_THRESHOLD LFC_THRESHOLD GSEA_SET_MIN GSEA_SET_MAX GSEA_NPERM

# Convert COMPARISONS array to a delimited string for passing to Python
# Build COMPARISONS_STR from array only if not already set directly in config
if [ -z "${COMPARISONS_STR:-}" ]; then
    COMPARISONS_STR="$(IFS='|||'; echo "${COMPARISONS[*]:-}")"
fi
export COMPARISONS_STR

# ---------------------------------------------------------------------------
# 0. Check all required tools and scripts are present
# ---------------------------------------------------------------------------
echo "=============================================================="
echo " RNA-Seq Pipeline"
echo " $(date)"
echo "=============================================================="
echo ""
echo "[0/6] Checking tools and scripts..."

ALL_OK=true

for tool in salmon Rscript python3 java wget; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "  ✓ ${tool}"
    else
        echo "  ✗ ${tool}  NOT FOUND"
        ALL_OK=false
    fi
done

for f in "${R_SCRIPT}" "${PY_RNK}"; do
    if [ -f "$f" ]; then
        echo "  ✓ $(basename "$f")"
    else
        echo "  ✗ $(basename "$f")  NOT FOUND at ${f}"
        ALL_OK=false
    fi
done

if [ "${ALL_OK}" = false ]; then
    echo ""
    echo "ERROR: Missing tools or scripts."
    echo "  Run:  bash .devcontainer/setup.sh"
    exit 1
fi

# Warn if COMPARISONS is empty
if [ -z "${COMPARISONS_STR}" ]; then
    echo ""
    echo "ERROR: COMPARISONS is empty. Set it in config/config.sh."
    echo "  Example:"
    echo '    export COMPARISONS=("DrugA,Vehicle,DrugA_vs_Vehicle")'
    exit 1
fi

echo ""
echo "  FASTQ_DIR:   ${FASTQ_DIR}"
echo "  OUTDIR:      ${OUTDIR}"
echo "  SAMPLESHEET: ${SAMPLESHEET}"
echo "  THREADS:     ${THREADS}"
echo "  PAIRED_END:  ${PAIRED_END}"
echo "  COMPARISONS:"
IFS='|||' read -ra _COMP_DISPLAY <<< "$COMPARISONS_STR"
for c in "${_COMP_DISPLAY[@]}"; do echo "    ${c}"; done
echo ""

mkdir -p "${OUTDIR}"/{salmon,deseq2,gsea,logs}

# ---------------------------------------------------------------------------
# 1. Build Salmon index (skipped if already present)
# ---------------------------------------------------------------------------
echo "[1/6] Salmon index..."

if [ -d "${SALMON_INDEX}" ] && [ -f "${SALMON_INDEX}/info.json" ]; then
    echo "  [SKIP] Index already exists: ${SALMON_INDEX}"
else
    if [ ! -f "${TRANSCRIPTOME_FA}" ]; then
        echo "ERROR: TRANSCRIPTOME_FA not found: ${TRANSCRIPTOME_FA}"
        echo "  Run:  bash .devcontainer/setup.sh"
        exit 1
    fi
    echo "  Building from ${TRANSCRIPTOME_FA} — takes ~10 minutes..."
    mkdir -p "${SALMON_INDEX}"
    salmon index \
        --transcripts "${TRANSCRIPTOME_FA}" \
        --index       "${SALMON_INDEX}" \
        --threads     "${THREADS}" \
        2>&1 | tee "${OUTDIR}/logs/salmon_index.log"
    echo "  Index built: ${SALMON_INDEX}"
fi

# ---------------------------------------------------------------------------
# 2. Parse samplesheet → /tmp/rnaseq_sample_meta.json
# ---------------------------------------------------------------------------
echo ""
echo "[2/6] Parsing samplesheet..."

python3 - <<'PYEOF'
import csv, json, os, sys

ss = os.environ["SAMPLESHEET"]
if not os.path.exists(ss):
    sys.exit(f"ERROR: Samplesheet not found: {ss}")

rows = []
with open(ss) as fh:
    for line in fh:
        if line.strip().startswith("#"):
            continue
    fh.seek(0)
    reader = csv.DictReader(
        (line for line in fh if not line.strip().startswith("#"))
    )
    for r in reader:
        row = {k.strip(): v.strip() for k, v in r.items()}
        donor = row.get("Donor_ID") or row.get("Donor_IDs") or "unknown"
        rows.append({
            "sample_id": row["Sample_ID"],
            "treatment": row["Treatment"],
            "condition": row["Condition"],
            "donor":     donor,
        })

with open("/tmp/rnaseq_sample_meta.json", "w") as fh:
    json.dump(rows, fh, indent=2)

donors = sorted(set(r["donor"] for r in rows))
print(f"  Samples: {len(rows)}  |  Donors: {len(donors)}")
for r in sorted(rows, key=lambda x: (x["donor"], x["treatment"])):
    print(f"    Donor {r['donor']:6s}  |  {r['treatment']:16s}  |  {r['sample_id']}")
PYEOF

# ---------------------------------------------------------------------------
# 3. Salmon quantification (per sample)
# ---------------------------------------------------------------------------
echo ""
echo "[3/6] Salmon quantification..."

python3 - <<'PYEOF'
import json, subprocess, os, pathlib, sys, glob

meta   = json.load(open("/tmp/rnaseq_sample_meta.json"))
fq_dir = os.environ["FASTQ_DIR"]
out_d  = os.environ["OUTDIR"]
idx    = os.environ["SALMON_INDEX"]
t      = os.environ["THREADS"]
paired = os.environ.get("PAIRED_END", "true").lower() == "true"
sal_d  = f"{out_d}/salmon"
log_d  = f"{out_d}/logs"

def find_fq(fq_dir, sid, read):
    for ext in [f"_R{read}.fastq.gz", f"_R{read}.fastq"]:
        p = os.path.join(fq_dir, f"{sid}{ext}")
        if os.path.exists(p):
            return p
    sid_under = sid.replace("-", "_")
    for pattern in [
        f"{sid_under}_S*_L*_R{read}_001.fastq",
        f"{sid_under}_S*_L*_R{read}_001.fastq.gz",
        f"{sid_under}*_R{read}*.fastq",
        f"{sid_under}*_R{read}*.fastq.gz",
    ]:
        matches = sorted(glob.glob(os.path.join(fq_dir, pattern)))
        if matches:
            return matches[0]
    return None

print(f"  Mode: {'paired-end (R1+R2)' if paired else 'single-end (R1 only)'}")

for s in meta:
    sid  = s["sample_id"]
    pfx  = f"{sal_d}/{sid}"
    done = f"{pfx}/quant.sf"

    if os.path.exists(done):
        print(f"  [SKIP] {sid}")
        continue

    r1 = find_fq(fq_dir, sid, 1)
    r2 = find_fq(fq_dir, sid, 2) if paired else None

    if not r1:
        print(f"  ERROR: R1 not found for {sid} in {fq_dir}")
        sys.exit(1)
    if paired and not r2:
        print(f"  ERROR: R2 not found for {sid} in {fq_dir}")
        sys.exit(1)

    print(f"  Quantifying {sid}...")
    print(f"    R1: {os.path.basename(r1)}")
    if paired:
        print(f"    R2: {os.path.basename(r2)}")

    pathlib.Path(pfx).mkdir(parents=True, exist_ok=True)

    cmd = [
        "salmon", "quant",
        "--index",          idx,
        "--libType",        "A",
        "--output",         pfx,
        "--threads",        t,
        "--gcBias",
        "--seqBias",
        "--validateMappings",
        "--numBootstraps",  "10",
    ]
    cmd += (["-1", r1, "-2", r2] if paired else ["-r", r1])

    log_f = f"{log_d}/{sid}_salmon.log"
    with open(log_f, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  ERROR: Salmon failed. See {log_f}")
        sys.exit(1)
    print(f"    ✓ {done}")

print("  All samples complete.")
PYEOF

# ---------------------------------------------------------------------------
# 4. DESeq2 — per donor × per comparison
# ---------------------------------------------------------------------------
echo ""
echo "[4/6] DESeq2..."

python3 - <<'PYEOF'
import json, subprocess, os, pathlib, sys

meta      = json.load(open("/tmp/rnaseq_sample_meta.json"))
out_dir   = os.environ["OUTDIR"]
gtf       = os.environ["GTF"]
repo_root = os.environ.get("REPO_ROOT", os.path.expanduser("~/rnaseq-pipeline"))
r_script  = f"{repo_root}/r_scripts/run_deseq2_salmon.R"
salmon_d  = f"{out_dir}/salmon"
deseq_d   = f"{out_dir}/deseq2"

comp_str  = os.environ.get("COMPARISONS_STR", "")
comps     = [c.strip().split(",") for c in comp_str.split("|||")
             if c.strip() and len(c.strip().split(",")) == 3]

donors    = sorted(set(s["donor"] for s in meta))
total     = len(donors) * len(comps)
print(f"  {len(donors)} donors × {len(comps)} comparisons = {total} DESeq2 jobs")

for donor in donors:
    by_trt = {s["treatment"]: s["sample_id"]
              for s in meta if s["donor"] == donor}
    out_d  = f"{deseq_d}/donor_{donor}"
    pathlib.Path(out_d).mkdir(parents=True, exist_ok=True)

    for trt, ctrl, label in comps:
        treat_s = by_trt.get(trt)
        ctrl_s  = by_trt.get(ctrl)

        if not treat_s or not ctrl_s:
            print(f"  [SKIP] Donor {donor} | {label} — sample(s) not found in samplesheet")
            continue

        out_csv = f"{out_d}/{label}_donor{donor}_deseq2_out.csv"
        if os.path.exists(out_csv):
            print(f"  [SKIP] {label} donor {donor} — already done")
            continue

        pm = {"treat_sample": treat_s, "ctrl_sample": ctrl_s,
              "donor": donor, "group_label": label,
              "salmon_dir": salmon_d, "gtf_path": gtf}
        pm_f = f"/tmp/pm_{donor}_{label}.json"
        with open(pm_f, "w") as fh:
            json.dump(pm, fh)

        print(f"\n  DESeq2  Donor {donor}  |  {label}")
        print(f"    treatment={treat_s}   control={ctrl_s}")
        env = {**os.environ,
               "SALMON_INDEX":     os.environ.get("SALMON_INDEX",""),
               "TRANSCRIPTOME_FA": os.environ.get("TRANSCRIPTOME_FA",""),
               "PADJ_THRESHOLD":   os.environ.get("PADJ_THRESHOLD","0.05"),
               "LFC_THRESHOLD":    os.environ.get("LFC_THRESHOLD","1")}
        r = subprocess.run(["Rscript", r_script, pm_f, out_d], env=env)
        if r.returncode != 0:
            print(f"  ERROR: DESeq2 failed. Check {out_dir}/logs/")
            sys.exit(1)

print("\n  DESeq2 all done.")
PYEOF

# ---------------------------------------------------------------------------
# 5. GSEA PreRanked — per donor × per comparison
# ---------------------------------------------------------------------------
echo ""
echo "[5/6] GSEA PreRanked..."

python3 - <<'PYEOF'
import json, subprocess, os, pathlib, sys

meta      = json.load(open("/tmp/rnaseq_sample_meta.json"))
out_dir   = os.environ["OUTDIR"]
gsea_jar  = os.environ["GSEA_JAR"]
gmt       = os.environ["HALLMARK_GMT"]
repo_root = os.environ.get("REPO_ROOT", os.path.expanduser("~/rnaseq-pipeline"))
py_rnk    = f"{repo_root}/python_scripts/make_rnk.py"
deseq_d   = f"{out_dir}/deseq2"
gsea_d    = f"{out_dir}/gsea"

comp_str  = os.environ.get("COMPARISONS_STR", "")
comp_lbls = [c.strip().split(",")[2].strip()
             for c in comp_str.split("|||")
             if c.strip() and len(c.strip().split(",")) == 3]

donors = sorted(set(s["donor"] for s in meta))

for tool_f, name in [(gsea_jar, "GSEA jar"), (gmt, "GMT file")]:
    if not os.path.exists(tool_f):
        print(f"  ERROR: {name} not found: {tool_f}")
        print("         Run:  bash .devcontainer/setup.sh")
        sys.exit(1)

for donor in donors:
    for label in comp_lbls:
        deseq_csv  = f"{deseq_d}/donor_{donor}/{label}_donor{donor}_deseq2_out.csv"
        run_label  = f"{label}_donor{donor}"
        gsea_run_d = f"{gsea_d}/donor_{donor}/{label}"
        rnk_file   = f"{gsea_run_d}/{run_label}.rnk"
        done_mark  = f"{gsea_run_d}/.gsea_done"

        if not os.path.exists(deseq_csv):
            print(f"  [SKIP] DESeq2 output missing: donor {donor} | {label}")
            continue

        pathlib.Path(gsea_run_d).mkdir(parents=True, exist_ok=True)

        print(f"\n  Building .rnk: donor {donor} | {label}")
        subprocess.run(["python3", py_rnk, deseq_csv, rnk_file], check=True)

        if os.path.exists(done_mark):
            print(f"  [SKIP] GSEA already complete: {run_label}")
            continue

        print(f"  Running GSEA: {run_label}")
        cmd = ["java", "-Xmx8g", "-cp", gsea_jar,
               "xtools.gsea.GseaPreranked",
               "-gmx",            gmt,
               "-rnk",            rnk_file,
               "-out",            gsea_run_d,
               "-rpt_label",      run_label,
               "-nperm",          os.environ.get("GSEA_NPERM","1000"),
               "-scoring_scheme", "weighted",
               "-set_min",        os.environ.get("GSEA_SET_MIN","15"),
               "-set_max",        os.environ.get("GSEA_SET_MAX","500"),
               "-create_svgs",    "false",
               "-make_sets",      "true",
               "-plot_top_x",     "20",
               "-norm",           "meandiv",
               "-zip_report",     "false",
               "-gui",            "false"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            open(done_mark, "w").close()
            print(f"    ✓ {gsea_run_d}")
        else:
            print(f"  ERROR: GSEA failed for {run_label}")
            print(r.stderr[-2000:])

print("\n  GSEA all done.")
PYEOF

# ---------------------------------------------------------------------------
# 6. Summary table
# ---------------------------------------------------------------------------
echo ""
echo "[6/6] Summary..."

python3 - <<'PYEOF'
import json, os
import pandas as pd

out_dir   = os.environ["OUTDIR"]
deseq_d   = f"{out_dir}/deseq2"
meta      = json.load(open("/tmp/rnaseq_sample_meta.json"))
donors    = sorted(set(s["donor"] for s in meta))
comp_str  = os.environ.get("COMPARISONS_STR","")
comp_lbls = [c.strip().split(",")[2].strip()
             for c in comp_str.split("|||")
             if c.strip() and len(c.strip().split(",")) == 3]
padj_t    = float(os.environ.get("PADJ_THRESHOLD","0.05"))

rows = []
for donor in donors:
    for label in comp_lbls:
        f = f"{deseq_d}/donor_{donor}/{label}_donor{donor}_deseq2_out.csv"
        if os.path.exists(f):
            df  = pd.read_csv(f)
            sig = df[df["significant"] == True] if "significant" in df.columns \
                  else df[df["padj"].fillna(1) < padj_t]
            rows.append({"donor": donor, "comparison": label,
                         "total_genes": len(df), "sig_genes": len(sig),
                         "up": len(sig[sig["log2FoldChange"] > 0]),
                         "down": len(sig[sig["log2FoldChange"] < 0]),
                         "status": "DONE"})
        else:
            rows.append({"donor": donor, "comparison": label,
                         "total_genes": 0, "sig_genes": 0,
                         "up": 0, "down": 0, "status": "MISSING"})

sumdf = pd.DataFrame(rows)
out_f = f"{deseq_d}/DESeq2_summary.csv"
sumdf.to_csv(out_f, index=False)
print(sumdf.to_string(index=False))
print(f"\n  Summary: {out_f}")
PYEOF

echo ""
echo "=============================================================="
echo " Pipeline complete — $(date)"
echo ""
echo " Outputs:"
echo "   ${OUTDIR}/salmon/      quant.sf per sample"
echo "   ${OUTDIR}/deseq2/      DESeq2 CSVs + DESeq2_summary.csv"
echo "   ${OUTDIR}/gsea/        GSEA HTML reports + .rnk files"
echo "   ${OUTDIR}/pipeline.log full log"
echo "=============================================================="
