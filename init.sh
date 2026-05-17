#!/bin/bash

# Initialize directory structure for dbsuper_pipeline.
#
# Project tree:
#   out/, raw/, ref/, logs/, Enhancerflow/   — Nextflow inputs/outputs/logs
#
# Enhancerflow sub-tree (read by submit_staging.sh and submit_enhancerflow.sh —
# we pre-create them here so SLURM's `-o logs/.../X.log` directive in those
# sbatch scripts has a directory to write into BEFORE the script body runs):
#   Enhancerflow/staging/         — per-sample case/control/peaks/bigwig bundles
#   Enhancerflow/logs/staging/    — per-sample stage_sample.sh logs
#   Enhancerflow/logs/driver/     — submit_enhancerflow.sh driver logs
#   Enhancerflow/logs/samplesheet/— samplesheet generation logs
#   Enhancerflow/logs/nextflow/   — nextflow .log
#   Enhancerflow/logs/reports/    — nextflow trace/report/timeline/dag
#
# Per-task scratch + container paths (referenced by nfcore_chipseq.config and
# Enhancerflow/nextflow.*.config via env.TMPDIR / env.APPTAINER_TMPDIR /
# env.APPTAINER_CACHEDIR):
#   /scratch/$USER/tmp            — TMPDIR for nextflow tasks
#   /scratch/$USER/apptainer_tmp  — APPTAINER_TMPDIR for container scratch
#   ~/.singularity_cache          — APPTAINER_CACHEDIR (kept under the legacy

set -euo pipefail

mkdir -p out raw ref logs Enhancerflow
mkdir -p \
    Enhancerflow/staging \
    Enhancerflow/logs/staging \
    Enhancerflow/logs/driver \
    Enhancerflow/logs/samplesheet \
    Enhancerflow/logs/nextflow \
    Enhancerflow/logs/reports
mkdir -p "/scratch/$USER/tmp" "/scratch/$USER/apptainer_tmp"
mkdir -p "$HOME/.singularity_cache"

echo "Directories created:"
echo "  out, raw, ref, logs, Enhancerflow"
echo "  Enhancerflow/staging"
echo "  Enhancerflow/logs/{staging,driver,samplesheet,nextflow,reports}"
echo "  /scratch/$USER/tmp"
echo "  /scratch/$USER/apptainer_tmp"
echo "  $HOME/.singularity_cache"
