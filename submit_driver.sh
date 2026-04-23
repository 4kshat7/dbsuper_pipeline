#!/bin/bash
# ------------------------------------------------------------------------------
# SLURM driver script for nf-core/chipseq v2.1.0
#
# Runs the FULL samplesheet (all members combined) in a single Nextflow driver.
# Do not split by member — this script processes
#   raw/encode_results/nfcore_chipseq_samplesheet.csv in one shot.
#
# Usage:
#   sbatch submit_driver.sh
#
# Real compute happens in child SLURM jobs that Nextflow submits via
# executor='slurm' configured in nfcore_chipseq.config.
# ------------------------------------------------------------------------------

#SBATCH -J nf-chipseq-full
#SBATCH -p cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH -o /nfs-stor/akshat.mistry/cluster/dbsuper_pipeline/logs/driver/driver-%x-%j.log
#SBATCH -e /nfs-stor/akshat.mistry/cluster/dbsuper_pipeline/logs/driver/driver-%x-%j.err
# NOTE: do NOT set --nodelist here. Driver can run anywhere; children spread freely.

set -euo pipefail

PROJECT_ROOT="/nfs-stor/$USER/cluster/dbsuper_pipeline"

cd "$PROJECT_ROOT"
mkdir -p logs/driver logs/nextflow logs/reports

# ---- Environment activation ------------------------------------------------
# CIAI cluster: conda is installed at /apps/local/anaconda3; login shells pick
# it up via /apps/local/conda_init.sh (same file referenced from ~/.bashrc).
# Disable -u around activation: conda hooks (e.g. openjdk_activate.sh) read
# $JAVA_HOME before it is set, which trips nounset.
set +u
source /apps/local/conda_init.sh 2>/dev/null \
  || source /apps/local/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || { echo "Could not find conda init; adjust path."; exit 1; }
conda activate encodefetch
set -u

# ---- Sanity checks ----------------------------------------------------------
SAMPLESHEET="raw/encode_results/nfcore_chipseq_samplesheet.local.csv"
CONFIG="nfcore_chipseq.config"
OUTDIR="out/full"

[[ -f "$SAMPLESHEET" ]] || { echo "Missing $SAMPLESHEET"; exit 1; }
[[ -f "$CONFIG"      ]] || { echo "Missing $CONFIG";      exit 1; }

echo "[$(date)] Launching nf-core/chipseq driver (full samplesheet)"
echo "  project     : $PROJECT_ROOT"
echo "  samplesheet : $SAMPLESHEET"
echo "  config      : $CONFIG"
echo "  outdir      : $OUTDIR"

# ---- Launch -----------------------------------------------------------------
# -log redirects Nextflow's own .nextflow.log into logs/nextflow/
#       (Nextflow auto-rotates to .log.1, .log.2, ... on each run)
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

echo "[$(date)] Driver exit status: $?"
