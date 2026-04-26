#!/bin/bash
# ------------------------------------------------------------------------------
# SLURM driver for the dbsuper_pipeline (member4 variant).
#
# Pipeline order inside the driver:
#   1. init.sh                                (once, if not already done)
#   2. conda activate encodefetch
#   3. fetch_member_accessions.sh             (manifest [+ optional FASTQ])
#   4. make_nfcode_samplesheet.py             (only if --download was given)
#   5. nextflow run nf-core/chipseq ...       (only if --download was given)
#
# Usage:
#   sbatch submit_driver.sh --member <1-5> \
#     [--accessions ACC1,ACC2,...] [--download] [--skip-fetch]
#
# Without --download the driver stops after step 3 (manifest/samplesheet only).
# With --skip-fetch, steps 3 and 4 are skipped entirely -- go straight to
# nextflow -resume. Use this when nothing has changed and you just want
# nextflow to pick up where it left off.
#
# Each step's stdout+stderr is captured in its own logs/<step>/ subdir with a
# timestamp so multiple submissions don't overwrite each other.
#
# Idempotency / re-run safety:
#   - init.sh is guarded by a .init_done marker (skipped on re-run).
#   - Before encodefetch runs, any *.fastq.gz under raw/encode_results/files/
#     is renamed back to *.fastq so encodefetch's built-in skip-if-present
#     check (which keys on .fastq) finds existing files and skips them.
#     make_nfcode_samplesheet.py renames them back to .fastq.gz afterwards.
# ------------------------------------------------------------------------------

#SBATCH -J nf-chipseq-m4
#SBATCH -p cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH -o logs/driver/driver-%x-%j.log
#SBATCH -e logs/driver/driver-%x-%j.err
#SBATCH --nodelist=cn-08

set -euo pipefail

# ── parse args ────────────────────────────────────────────────────────────────
MEMBER=""
ACCESSIONS=""
DOWNLOAD=false
SKIP_FETCH=false

usage() {
  echo "Usage: sbatch $(basename "$0") [--member <1-5|all>] \\"
  echo "         [--accessions ACC1,ACC2,...] [--download] [--skip-fetch]"
  echo ""
  echo "  Omitting --member (or passing --member all) fetches the union of all"
  echo "  five member files = the full 913-accession dataset."
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --member)      MEMBER="$2";     shift 2 ;;
    --accessions)  ACCESSIONS="$2"; shift 2 ;;
    --download)    DOWNLOAD=true;   shift   ;;
    --skip-fetch)  SKIP_FETCH=true; shift   ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Default: no --member means run on the full dataset (union of all 5 members).
[[ -z "$MEMBER" ]] && MEMBER="all"

# ── project layout ────────────────────────────────────────────────────────────
# Under sbatch, SLURM copies this script to a staging dir on the compute node,
# so BASH_SOURCE[0] points somewhere read-only. SLURM_SUBMIT_DIR is the dir you
# ran sbatch from (the project root, as long as you sbatch from here).
# Falls back to BASH_SOURCE for manual `bash submit_driver.sh` runs.
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJECT_ROOT"

TS="$(date +%Y%m%d-%H%M%S)"
mkdir -p logs/init logs/fetch logs/samplesheet logs/driver logs/nextflow logs/reports

INIT_LOG="logs/init/init-${TS}.log"
FETCH_LOG="logs/fetch/fetch-${MEMBER}-${TS}.log"
SAMPLESHEET_LOG="logs/samplesheet/samplesheet-${TS}.log"

echo "[$(date)] === driver start ==="
echo "  project    : $PROJECT_ROOT"
echo "  member     : $MEMBER"
echo "  accessions : ${ACCESSIONS:-<default list for member ${MEMBER}>}"
echo "  download   : $DOWNLOAD"
echo "  timestamp  : $TS"

# ── 1. init.sh (once) ─────────────────────────────────────────────────────────
if [[ -f .init_done ]]; then
  echo "[$(date)] init.sh already ran (.init_done present); skipping"
else
  echo "[$(date)] running init.sh -> $INIT_LOG"
  bash init.sh >"$INIT_LOG" 2>&1
  touch .init_done
fi

# ── 2. conda activate encodefetch ─────────────────────────────────────────────
# CIAI cluster: conda at /apps/local/anaconda3. Disable -u around activation
# because the openjdk activate hook reads $JAVA_HOME before setting it.
echo "[$(date)] activating conda env 'encodefetch'"
set +u
source /apps/local/conda_init.sh 2>/dev/null \
  || source /apps/local/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || { echo "Could not find conda init; adjust path."; exit 1; }
conda activate encodefetch
set -u

SAMPLESHEET="raw/encode_results/nfcore_chipseq_samplesheet.local.csv"

if $SKIP_FETCH; then
  echo "[$(date)] --skip-fetch given; skipping fetch + samplesheet steps."
  [[ -f "$SAMPLESHEET" ]] || { echo "ERROR: --skip-fetch but $SAMPLESHEET not found"; exit 1; }
else
  # ── 3a. un-rename .fastq.gz -> .fastq so encodefetch skip-existing works ────
  # First run: no *.fastq.gz present, this is a no-op.
  # Later runs: make_nfcode_samplesheet.py renamed everything to .fastq.gz;
  # encodefetch's skip check keys on the .fastq suffix, so we flip back first.
  if [[ -d raw/encode_results/files ]]; then
    GZ_COUNT=$(find raw/encode_results/files -name '*.fastq.gz' -printf . 2>/dev/null | wc -c)
    if (( GZ_COUNT > 0 )); then
      echo "[$(date)] un-renaming $GZ_COUNT *.fastq.gz -> *.fastq (so encodefetch can skip-existing)"
      find raw/encode_results/files -name '*.fastq.gz' -exec bash -c '
        for f; do mv "$f" "${f%.gz}"; done
      ' _ {} +
    fi
  fi

  # ── 3b. fetch_member_accessions.sh ────────────────────────────────────────
  FETCH_ARGS=(--member "$MEMBER")
  [[ -n "$ACCESSIONS" ]] && FETCH_ARGS+=(--accessions "$ACCESSIONS")
  $DOWNLOAD && FETCH_ARGS+=(--download)

  echo "[$(date)] running fetch_member_accessions.sh ${FETCH_ARGS[*]} -> $FETCH_LOG"
  bash fetch_member_accessions.sh "${FETCH_ARGS[@]}" >"$FETCH_LOG" 2>&1

  # ── gate: stop here if --download was not given ──────────────────────────
  if ! $DOWNLOAD; then
    echo "[$(date)] --download not given; stopping after fetch step."
    echo "  manifest + URL samplesheet available under raw/encode_results/"
    echo "  to continue: resubmit with --download"
    exit 0
  fi

  # ── 4. make_nfcode_samplesheet.py ────────────────────────────────────────
  # Script assumes cwd == subset dir (it picks up raw/encode_results/* from '.').
  # It re-renames *.fastq -> *.fastq.gz and rebuilds the .local.csv sheet.
  echo "[$(date)] running make_nfcode_samplesheet.py -> $SAMPLESHEET_LOG"
  (
    cd raw/encode_results
    python3 make_nfcode_samplesheet.py
  ) >"$SAMPLESHEET_LOG" 2>&1

  [[ -f "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET not produced"; exit 1; }
fi

# ── 5. nextflow run ───────────────────────────────────────────────────────────
CONFIG="nfcore_chipseq.config"
if [[ "$MEMBER" == "all" ]]; then
  OUTDIR="out/all"
else
  OUTDIR="out/member${MEMBER}"
fi

[[ -f "$CONFIG" ]] || { echo "ERROR: missing $CONFIG"; exit 1; }

echo "[$(date)] launching nf-core/chipseq"
echo "  samplesheet : $SAMPLESHEET"
echo "  config      : $CONFIG"
echo "  outdir      : $OUTDIR"

nextflow -log logs/nextflow/.nextflow.log \
  run nf-core/chipseq -r 2.1.0 \
  -profile singularity \
  -c "$CONFIG" \
  -resume \
  --input   "$SAMPLESHEET" \
  --outdir  "$OUTDIR" \
  --genome  GRCh38 \
  --read_length 50 \
  --narrow_peak

echo "[$(date)] driver exit status: $?"
