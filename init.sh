#!/bin/bash

# Initialize directory structure for dbsuper_pipeline.
#
# Project tree:
#   out/, raw/, ref/, logs/, Enhancerflow/   — Nextflow inputs/outputs/logs
#
# Per-task scratch + container paths (referenced by nfcore_chipseq.config and
# Enhancerflow/nextflow.*.config via env.TMPDIR / env.APPTAINER_TMPDIR /
# env.APPTAINER_CACHEDIR):
#   /scratch/$USER/tmp            — TMPDIR for nextflow tasks
#   /scratch/$USER/apptainer_tmp  — APPTAINER_TMPDIR for container scratch
#   ~/.singularity_cache          — APPTAINER_CACHEDIR (kept under the legacy

set -euo pipefail

mkdir -p out raw ref logs Enhancerflow
mkdir -p "/scratch/$USER/tmp" "/scratch/$USER/apptainer_tmp"
mkdir -p "$HOME/.singularity_cache"

echo "Directories created:"
echo "  out, raw, ref, logs, Enhancerflow"
echo "  /scratch/$USER/tmp"
echo "  /scratch/$USER/apptainer_tmp"
echo "  $HOME/.singularity_cache"
