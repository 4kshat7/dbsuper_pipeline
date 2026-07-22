#!/bin/bash
# ------------------------------------------------------------------------------
# SLURM driver for the dbsuper_pipeline.
#
# Pipeline order inside the driver:
#   1. init.sh                                (once, if not already done)
#   2. module load nextflow/25.10.4 encodefetch/0.5.0 apptainer/system
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

#SBATCH --job-name=nf-chipseq
#SBATCH --partition=cpu
#SBATCH --account=$USER
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH -o logs/driver/driver-%x-%j.log
#SBATCH -e logs/driver/driver-%x-%j.err

set -euo pipefail

# ── parse args ────────────────────────────────────────────────────────────────
MEMBER=""
ACCESSIONS=""
DOWNLOAD=false
SKIP_FETCH=false
SAMPLESHEET_OVERRIDE=""              # override the hardcoded input samplesheet path
OUTDIR_OVERRIDE=""                   # override the member-derived nextflow --outdir
BWA_INDEX_OVERRIDE=""                # pass --bwa_index (needed for the flat local mm10 index)

# ── genome selection (EDIT HERE) ────────────────────────────────────────────
# Which iGenomes build to align against, and its MACS2 effective genome size.
# Override at submit time with --genome <key> (and optionally --macs-gsize N).
# To support a new build, add a row to GSIZE_BY_GENOME below.
#
# gsize values are deepTools "effective genome size" (non-N bases) — the same
# convention as the original GRCh38 value this config shipped with:
#   https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GENOME="mm10"                       # default: mouse (UCSC mm10 == Ensembl GRCm38)
MACS_GSIZE=""                       # blank = auto-pick from GSIZE_BY_GENOME[$GENOME]
declare -A GSIZE_BY_GENOME=(
  [mm10]=2652783500                 # mouse — UCSC mm10 (== GRCm38)
  [GRCm38]=2652783500               # mouse — Ensembl GRCm38
  [GRCm37]=2620345972               # mouse — GRCm37 / mm9
  [GRCh38]=2913022398               # human — GRCh38
  [GRCh37]=2864785220               # human — GRCh37 / hg19
)

usage() {
  echo "Usage: sbatch $(basename "$0") [--member <1-5|all>] \\"
  echo "         [--accessions ACC1,ACC2,...] [--genome KEY] [--macs-gsize N] \\"
  echo "         [--samplesheet PATH] [--outdir NAME] [--bwa-index DIR] \\"
  echo "         [--download] [--skip-fetch]"
  echo ""
  echo "  Omitting --member (or passing --member all) fetches the union of all"
  echo "  five member files = the full 913-accession dataset."
  echo ""
  echo "  --genome      iGenomes build to align against (default: ${GENOME})."
  echo "                Known keys: ${!GSIZE_BY_GENOME[*]}"
  echo "  --macs-gsize  Override MACS2 effective genome size (default: auto from --genome)."
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --member)      MEMBER="$2";     shift 2 ;;
    --accessions)  ACCESSIONS="$2"; shift 2 ;;
    --genome)      GENOME="$2";     shift 2 ;;
    --macs-gsize)  MACS_GSIZE="$2"; shift 2 ;;
    --download)    DOWNLOAD=true;   shift   ;;
    --skip-fetch)  SKIP_FETCH=true; shift   ;;
    --samplesheet) SAMPLESHEET_OVERRIDE="$2"; shift 2 ;;
    --outdir)      OUTDIR_OVERRIDE="$2";      shift 2 ;;
    --bwa-index)   BWA_INDEX_OVERRIDE="$2";   shift 2 ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Default: no --member means run on the full dataset (union of all 5 members).
[[ -z "$MEMBER" ]] && MEMBER="all"

# MACS2 effective genome size. Source of truth is params.macs_gsize in
# nfcore_chipseq.config -- edit it there. The driver leaves it alone unless you
# explicitly pass --macs-gsize N, in which case that value overrides the config.
# (GSIZE_BY_GENOME above is kept as a reference table / set of valid --genome keys.)

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
echo "  genome     : $GENOME (macs_gsize=${MACS_GSIZE:-from nfcore_chipseq.config})"
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

# ── 2. load environment modules ───────────────────────────────────────────────
# Karakoram uses Lmod (modules under /storage/software/modules). The `module`
# function isn't defined in non-login shells (sbatch jobs), so source the init
# scripts explicitly. Disable -u around `module` because Lmod's bash function
# reads unset variables internally.
#
#   apptainer/system   → singularity/apptainer runtime for nf-core containers
#   nextflow/25.10.4   → Nextflow CLI + bundled OpenJDK 23 (matches old env)
#   encodefetch/0.5.0  → encodefetch CLI + Python 3.11 with pandas (used by
#                        make_nfcode_samplesheet.py)
#
# Load order matters: apptainer prepends /usr/bin to PATH (shadows env pythons),
# so it must be loaded BEFORE encodefetch — otherwise `python3` resolves to
# /usr/bin/python3 which lacks pandas.
echo "[$(date)] loading environment modules"
set +u
source /etc/profile.d/lmod.sh
source /etc/profile.d/modules.sh
module load apptainer/system nextflow/25.10.4 encodefetch/0.5.0
set -u

SAMPLESHEET="raw/encode_results/nfcore_chipseq_samplesheet.local.csv"
[[ -n "$SAMPLESHEET_OVERRIDE" ]] && SAMPLESHEET="$SAMPLESHEET_OVERRIDE"

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
[[ -n "$OUTDIR_OVERRIDE" ]] && OUTDIR="$OUTDIR_OVERRIDE"

[[ -f "$CONFIG" ]] || { echo "ERROR: missing $CONFIG"; exit 1; }

echo "[$(date)] launching nf-core/chipseq"
echo "  samplesheet : $SAMPLESHEET"
echo "  config      : $CONFIG"
echo "  outdir      : $OUTDIR"
echo "  genome      : $GENOME (macs_gsize=${MACS_GSIZE:-from nfcore_chipseq.config})"
echo "  bwa_index   : ${BWA_INDEX_OVERRIDE:-from igenomes ($GENOME key)}"

# --skip_spp: some GEO/legacy samples align at ~0% (bad chemistry/library, not a
# pipeline bug) and leave phantompeakqualtools nothing to read, which crashes it
# and (since its errorStrategy isn't 'retry'-eligible) aborts the whole run. Its
# NSC/RSC/correlation output feeds MultiQC only -- not MACS3 or anything
# downstream -- so skipping it is safe; those samples still yield ~0 peaks and
# get dropped later by make_enhancerflow_samplesheet.py's --min-peaks check.
nextflow -log logs/nextflow/.nextflow.log \
  run nf-core/chipseq -r 2.1.0 \
  -profile apptainer \
  -c "$CONFIG" \
  -resume \
  --input      "$SAMPLESHEET" \
  --outdir     "$OUTDIR" \
  --genome     "$GENOME" \
  ${MACS_GSIZE:+--macs_gsize "$MACS_GSIZE"} \
  ${BWA_INDEX_OVERRIDE:+--bwa_index "$BWA_INDEX_OVERRIDE"} \
  --narrow_peak \
  --skip_spp

echo "[$(date)] driver exit status: $?"
