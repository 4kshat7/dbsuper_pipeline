#!/bin/bash
# ------------------------------------------------------------------------------
# submit_staging.sh — single SLURM entry point for Enhancerflow staging.
#
# What `sbatch submit_staging.sh` does:
#   1. Activate the dbsuper_pipeline micromamba env (for Python).
#   2. Hand off to bin/stage_orchestrator.py, which:
#        - parses the nf-core/chipseq samplesheet
#        - folds _se / _pe ENCODE entries with the same ENCSR id into ONE
#          biological sample
#        - submits per-sample stage_sample.sh sbatch jobs (throttled to
#          --max-jobs concurrent in the queue)
#        - waits for them all to finish
#        - generates enhancerflow_samplesheet.csv
#
# Run it like:
#     sbatch submit_staging.sh                       # member4 (default)
#     sbatch submit_staging.sh --member 1
#     sbatch submit_staging.sh --skip-stage          # only regen samplesheet
#     sbatch submit_staging.sh --dry-run             # preview without submitting
#     sbatch submit_staging.sh --force               # re-stage even done samples
#     sbatch submit_staging.sh --clean-stale         # delete legacy _se/_pe dirs
#
# Anything after the recognised flags is passed straight through to the
# orchestrator, so `--max-jobs 30` (etc.) works as you'd expect.
#
# Driver sizing: the orchestrator runs the Nextflow JVM... no it doesn't,
# it just submits sbatch jobs and polls squeue. So 2 CPU / 4 GB is enough.
# We give it 1 day; if you have a huge backlog needing throttled submission
# over many hours, bump --time.
# ------------------------------------------------------------------------------

#SBATCH --job-name=stage-orchestrator
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH -o logs/staging-driver.log
#SBATCH -e logs/staging-driver.log

set -euo pipefail

# ── parse args ────────────────────────────────────────────────────────────────
MEMBER="4"
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --member) MEMBER="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: sbatch $(basename "$0") [--member N] [orchestrator args ...]"
      echo "Orchestrator args (passed through):"
      echo "  --skip-stage           regen samplesheet only"
      echo "  --dry-run              preview without submitting"
      echo "  --force                re-stage even completed samples"
      echo "  --clean-stale          delete staging dirs that don't match any base sample"
      echo "  --max-jobs N           cap concurrent stage-* jobs (default 50)"
      echo "  --poll-seconds N       queue poll interval (default 30)"
      exit 0
      ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

# ── project layout ────────────────────────────────────────────────────────────
ENHANCERFLOW_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ENHANCERFLOW_DIR"

CHIPSEQ_SAMPLESHEET="../raw/encode_results/nfcore_chipseq_samplesheet.local.csv"
CHIPSEQ_OUTDIR="../out/member${MEMBER}"
STAGING_DIR="staging"
LOGS_DIR="logs/staging"
OUTPUT_SAMPLESHEET="enhancerflow_samplesheet.csv"

mkdir -p "$LOGS_DIR" "$STAGING_DIR" logs

echo "=== submit_staging start  $(date '+%F %T') ==="
echo "  Enhancerflow dir   : $ENHANCERFLOW_DIR"
echo "  member             : $MEMBER"
echo "  chipseq samplesheet: $CHIPSEQ_SAMPLESHEET"
echo "  chipseq outdir     : $CHIPSEQ_OUTDIR"
echo "  staging dir        : $STAGING_DIR"
echo "  per-sample logs    : $LOGS_DIR/<SAMPLE>.log"
echo "  output samplesheet : $OUTPUT_SAMPLESHEET"
echo "  extra args         : ${EXTRA_ARGS[*]:-(none)}"

# ── activate env (just for python) ───────────────────────────────────────────
set +u
export MAMBA_ROOT_PREFIX="/storage/software/micromamba"
eval "$(/usr/local/bin/micromamba shell hook --shell bash)" \
  || { echo "could not init micromamba"; exit 1; }
micromamba activate dbsuper_pipeline
set -u

# ── run orchestrator ──────────────────────────────────────────────────────────
python3 bin/stage_orchestrator.py \
    --chipseq-samplesheet "$CHIPSEQ_SAMPLESHEET" \
    --chipseq-outdir      "$CHIPSEQ_OUTDIR" \
    --staging-dir         "$STAGING_DIR" \
    --logs-dir            "$LOGS_DIR" \
    --output-samplesheet  "$OUTPUT_SAMPLESHEET" \
    --condition           H3K27ac \
    --max-jobs            50 \
    --poll-seconds        30 \
    "${EXTRA_ARGS[@]}"

echo "=== submit_staging exit $? at $(date '+%F %T') ==="
