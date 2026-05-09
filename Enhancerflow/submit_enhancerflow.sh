#!/bin/bash
# ------------------------------------------------------------------------------
# SLURM driver for Enhancerflow (super-enhancer calling on chipseq output).
#
# Pipeline order:
#   1. module load apptainer/system nextflow/25.10.4 encodefetch/0.5.0
#   2. make_enhancerflow_samplesheet.py
#        chipseq samplesheet + chipseq out/ → enhancerflow_samplesheet.csv
#   3. nextflow run khan-lab/enhancerflow ...
#
# Layout (everything stays under Enhancerflow/):
#   logs/driver/      — driver job stdout/err
#   logs/samplesheet/ — make_enhancerflow_samplesheet.py logs
#   logs/nextflow/    — Nextflow .nextflow.log
#   logs/reports/     — trace/report/timeline/dag (set by nextflow.slurm.config)
#   out/              — Enhancerflow --outdir
#
# Usage:
#   sbatch submit_enhancerflow.sh                       # member1 (default)
#   sbatch submit_enhancerflow.sh --member 2
#   sbatch submit_enhancerflow.sh --member all
#   sbatch submit_enhancerflow.sh --member 1 --skip-samplesheet
#
# Driver sizing rationale: Nextflow's executor is 'slurm' (see
# nextflow.slurm.config), so all real compute is dispatched as child sbatch
# jobs on compute nodes. The driver only runs the Nextflow JVM + the python
# samplesheet step → 4 CPUs, 8 GB is plenty (validated empirically — JVM RES
# typically <2 GB, 4 GB peak).
# ------------------------------------------------------------------------------

#SBATCH --job-name=nf-enhancerflow
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH -o logs/driver/driver-%x-%j.log
#SBATCH -e logs/driver/driver-%x-%j.err

set -euo pipefail

# ── parse args ────────────────────────────────────────────────────────────────
MEMBER="1"
SKIP_SAMPLESHEET=false

usage() {
  echo "Usage: sbatch $(basename "$0") [--member <1-5|all>] [--skip-samplesheet]"
  echo ""
  echo "  --member            which chipseq output to consume (default: 1)"
  echo "                      maps to ../out/member<N> (or ../out/all)"
  echo "  --skip-samplesheet  reuse the existing enhancerflow_samplesheet.csv"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --member)            MEMBER="$2";       shift 2 ;;
    --skip-samplesheet)  SKIP_SAMPLESHEET=true; shift ;;
    -h|--help)           usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# ── project layout ────────────────────────────────────────────────────────────
# Under sbatch, $SLURM_SUBMIT_DIR is the dir you ran sbatch from. We expect
# that to be Enhancerflow/. Falls back to BASH_SOURCE for manual `bash` runs.
ENHANCERFLOW_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ENHANCERFLOW_DIR"

if [[ "$MEMBER" == "all" ]]; then
  CHIPSEQ_OUTDIR="../out/all"
else
  CHIPSEQ_OUTDIR="../out/member${MEMBER}"
fi
CHIPSEQ_SAMPLESHEET="../raw/encode_results/nfcore_chipseq_samplesheet.local.csv"
SAMPLESHEET="enhancerflow_samplesheet.csv"

TS="$(date +%Y%m%d-%H%M%S)"
mkdir -p logs/driver logs/samplesheet logs/nextflow logs/reports out

SAMPLESHEET_LOG="logs/samplesheet/samplesheet-${MEMBER}-${TS}.log"

echo "[$(date)] === enhancerflow driver start ==="
echo "  enhancerflow_dir : $ENHANCERFLOW_DIR"
echo "  member           : $MEMBER"
echo "  chipseq_outdir   : $CHIPSEQ_OUTDIR"
echo "  samplesheet_out  : $SAMPLESHEET"
echo "  skip-samplesheet : $SKIP_SAMPLESHEET"
echo "  timestamp        : $TS"

# ── 1. load environment modules ───────────────────────────────────────────────
# Karakoram uses Lmod; sbatch jobs run in non-login shells where the `module`
# function isn't defined, so source the init scripts explicitly. set +u around
# the calls because Lmod's bash function reads unset variables.
#
#   apptainer/system   → singularity/apptainer container runtime
#   nextflow/25.10.4   → Nextflow CLI + bundled OpenJDK 23
#   encodefetch/0.5.0  → Python 3.11 with pandas (used by
#                        make_enhancerflow_samplesheet.py)
#
# Load order matters: apptainer prepends /usr/bin to PATH (which would shadow
# the encodefetch env's python3), so it must come BEFORE encodefetch.
echo "[$(date)] loading environment modules"
set +u
source /etc/profile.d/lmod.sh
source /etc/profile.d/modules.sh
module load apptainer/system nextflow/25.10.4 encodefetch/0.5.0
set -u

# ── 2. make_enhancerflow_samplesheet.py ──────────────────────────────────────
if $SKIP_SAMPLESHEET; then
  echo "[$(date)] --skip-samplesheet given; reusing existing $SAMPLESHEET"
  [[ -f "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET not found"; exit 1; }
else
  [[ -f "$CHIPSEQ_SAMPLESHEET" ]] || {
    echo "ERROR: chipseq samplesheet not found at $CHIPSEQ_SAMPLESHEET"; exit 1; }
  [[ -d "$CHIPSEQ_OUTDIR" ]] || {
    echo "ERROR: chipseq outdir not found at $CHIPSEQ_OUTDIR"; exit 1; }

  echo "[$(date)] running make_enhancerflow_samplesheet.py -> $SAMPLESHEET_LOG"
  python3 make_enhancerflow_samplesheet.py \
    --chipseq-samplesheet "$CHIPSEQ_SAMPLESHEET" \
    --chipseq-outdir      "$CHIPSEQ_OUTDIR" \
    --absolute \
    >"$SAMPLESHEET_LOG" 2>&1

  [[ -s "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET empty"; exit 1; }
fi

# ── 3. nextflow run ───────────────────────────────────────────────────────────
CONFIG="nextflow.slurm.config"
[[ -f "$CONFIG" ]] || { echo "ERROR: missing $CONFIG"; exit 1; }

echo "[$(date)] launching khan-lab/enhancerflow"
echo "  samplesheet : $SAMPLESHEET"
echo "  config      : $CONFIG"
echo "  outdir      : out/"

nextflow -log logs/nextflow/.nextflow.log \
  run khan-lab/enhancerflow -r main \
  -profile apptainer \
  -c "$CONFIG" \
  -resume \
  --input          "$SAMPLESHEET" \
  --genome         hg38 \
  --outdir         out/ \
  --skip_motifs \
  --skip_comparison

echo "[$(date)] driver exit status: $?"
