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
#   sbatch submit_driver.sh --member <1-5> --macs-gsize N \
#     [--accessions ACC1,ACC2,...] [--gtf PATH] [--blacklist PATH] \
#     [--download] [--skip-fetch]
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
GTF_OVERRIDE=""                      # pass a local --gtf instead of the iGenomes GTF
FASTA_OVERRIDE=""                    # pass a local --fasta instead of the iGenomes (S3) FASTA
BLACKLIST_OVERRIDE=""                # pass --blacklist (the CHM13 iGenomes key defines none)
WORK_DIR_OVERRIDE=""                 # pass -w (separate nf_work/ task cache per dataset)

# ── genome selection ────────────────────────────────────────────────────────
# Override the iGenomes build at submit time with --genome <key>.
# The effective genome size must be supplied explicitly with --macs-gsize N.
GENOME="mm10"                       # default: mouse (UCSC mm10 == Ensembl GRCm38)
MACS_GSIZE=""

usage() {
  echo "Usage: sbatch $(basename "$0") [--member <1-5|all>] \\"
  echo "         --macs-gsize N [--accessions ACC1,ACC2,...] [--genome KEY] \\"
  echo "         [--samplesheet PATH] [--outdir NAME] [--bwa-index DIR] \\"
  echo "         [--fasta PATH] [--gtf PATH] [--blacklist PATH] \\"
  echo "         [--work-dir DIR] [--download] [--skip-fetch]"
  echo ""
  echo "  Omitting --member (or passing --member all) fetches the union of all"
  echo "  five member files = the full 913-accession dataset."
  echo ""
  echo "  --genome      iGenomes build to align against (default: ${GENOME})."
  echo "  --macs-gsize  Required positive integer: MACS2 effective genome size."
  echo "  --gtf         Optional local GTF; overrides the GTF from iGenomes."
  echo "  --fasta       Optional local genome FASTA; overrides the iGenomes one."
  echo "                STRONGLY RECOMMENDED: setting params.igenomes_base in the"
  echo "                -c config does NOT redirect iGenomes to local disk (the"
  echo "                paths are baked in at includeConfig time), so without this"
  echo "                the FASTA is pulled from s3://ngi-igenomes/ over the network."
  echo "  --blacklist   Optional exclusion-set BED. REQUIRED for --genome CHM13:"
  echo "                the CHM13 iGenomes key defines no blacklist, so without"
  echo "                this the run is unfiltered. See T2T note below."
  echo "  --work-dir    Optional Nextflow work dir (-w); overrides the workDir in"
  echo "                nfcore_chipseq.config so a dataset gets its own task cache."
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
    --gtf)         GTF_OVERRIDE="$2";         shift 2 ;;
    --fasta)       FASTA_OVERRIDE="$2";       shift 2 ;;
    --blacklist)   BLACKLIST_OVERRIDE="$2";   shift 2 ;;
    --work-dir)    WORK_DIR_OVERRIDE="$2";    shift 2 ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Require an explicit MACS2 effective genome size for every driver invocation.
if [[ -z "$MACS_GSIZE" ]]; then
  echo "ERROR: --macs-gsize is required." >&2
  exit 2
fi
if [[ ! "$MACS_GSIZE" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --macs-gsize must be a positive integer; got: $MACS_GSIZE" >&2
  exit 2
fi

# Default: no --member means run on the full dataset (union of all 5 members).
[[ -z "$MEMBER" ]] && MEMBER="all"

# ── project layout ────────────────────────────────────────────────────────────
# Under sbatch, SLURM copies this script to a staging dir on the compute node,
# so BASH_SOURCE[0] points somewhere read-only. SLURM_SUBMIT_DIR is the dir you
# ran sbatch from (the project root, as long as you sbatch from here).
# Falls back to BASH_SOURCE for manual `bash submit_driver.sh` runs.
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJECT_ROOT"

if [[ -n "$GTF_OVERRIDE" ]]; then
  [[ -f "$GTF_OVERRIDE" ]] || {
    echo "ERROR: --gtf file not found: $GTF_OVERRIDE" >&2
    exit 2
  }
  GTF_OVERRIDE="$(realpath "$GTF_OVERRIDE")"
fi

if [[ -n "$FASTA_OVERRIDE" ]]; then
  [[ -f "$FASTA_OVERRIDE" ]] || {
    echo "ERROR: --fasta file not found: $FASTA_OVERRIDE" >&2
    exit 2
  }
  FASTA_OVERRIDE="$(realpath "$FASTA_OVERRIDE")"
fi

if [[ -n "$BLACKLIST_OVERRIDE" ]]; then
  [[ -f "$BLACKLIST_OVERRIDE" ]] || {
    echo "ERROR: --blacklist file not found: $BLACKLIST_OVERRIDE" >&2
    exit 2
  }
  BLACKLIST_OVERRIDE="$(realpath "$BLACKLIST_OVERRIDE")"
fi

# The CHM13 iGenomes key carries no blacklist (unlike GRCh38/hg38, which default
# to ENCODE hg38-blacklist.v3). T2T assembles the centromeric alpha-satellite
# arrays and acrocentric short arms that hg38 left out, and those become
# high-signal multi-mapping sinks -- exactly the rank-1 artefact that inflates
# ROSE2's super-enhancer cutoff downstream. Warn loudly rather than fail, so
# deliberate unfiltered runs stay possible.
if [[ "$GENOME" == "CHM13" && -z "$BLACKLIST_OVERRIDE" ]]; then
  echo "WARNING: --genome CHM13 without --blacklist; no exclusion set will be applied." >&2
  echo "         Suggested: --blacklist /storage/datasets/igenomes/Homo_sapiens/NCBI/CHM13/Annotation/Blacklist/T2T.excluderanges.bed" >&2
fi

TS="$(date +%Y%m%d-%H%M%S)"
mkdir -p logs/init logs/fetch logs/samplesheet logs/driver logs/nextflow logs/reports

INIT_LOG="logs/init/init-${TS}.log"
FETCH_LOG="logs/fetch/fetch-${MEMBER}-${TS}.log"
SAMPLESHEET_LOG="logs/samplesheet/samplesheet-${TS}.log"

echo "[$(date)] === driver start ==="
echo "  project    : $PROJECT_ROOT"
echo "  member     : $MEMBER"
echo "  accessions : ${ACCESSIONS:-<default list for member ${MEMBER}>}"
echo "  genome     : $GENOME (macs_gsize=$MACS_GSIZE)"
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
echo "  genome      : $GENOME (macs_gsize=$MACS_GSIZE)"
echo "  fasta       : ${FASTA_OVERRIDE:-from igenomes ($GENOME key) -- WARNING: S3 download}"
echo "  gtf         : ${GTF_OVERRIDE:-from igenomes ($GENOME key)}"
echo "  bwa_index   : ${BWA_INDEX_OVERRIDE:-from igenomes ($GENOME key)}"
echo "  blacklist   : ${BLACKLIST_OVERRIDE:-from igenomes ($GENOME key; none for CHM13)}"
echo "  work_dir    : ${WORK_DIR_OVERRIDE:-from $CONFIG (workDir)}"

# phantompeakqualtools (SPP) runs: it gives every sample its NSC, RSC and
# fragment length, which dbSUPER stores in experiment.more_info.qc. It used to
# be turned off here with --skip_spp (93b627f, for the GEO batch), because a
# sample aligning at ~0% leaves SPP nothing to read and its crash aborted the
# whole run; every run since then has null NSC/RSC. A failing SPP task is now
# dropped for that sample alone -- see PHANTOMPEAKQUALTOOLS in
# nfcore_chipseq.config. SPP feeds MultiQC only, never MACS3 or anything
# downstream, so turning it on changes no peaks or super-enhancers; ~0%-aligned
# samples still yield ~0 peaks and are dropped by make_enhancerflow_samplesheet.py's
# --min-peaks check.
nextflow -log logs/nextflow/.nextflow.log \
  run nf-core/chipseq -r 2.1.0 \
  -profile apptainer \
  -c "$CONFIG" \
  ${WORK_DIR_OVERRIDE:+-w "$WORK_DIR_OVERRIDE"} \
  -resume \
  --input      "$SAMPLESHEET" \
  --outdir     "$OUTDIR" \
  --genome     "$GENOME" \
  --macs_gsize "$MACS_GSIZE" \
  ${FASTA_OVERRIDE:+--fasta "$FASTA_OVERRIDE"} \
  ${GTF_OVERRIDE:+--gtf "$GTF_OVERRIDE"} \
  ${BWA_INDEX_OVERRIDE:+--bwa_index "$BWA_INDEX_OVERRIDE"} \
  ${BLACKLIST_OVERRIDE:+--blacklist "$BLACKLIST_OVERRIDE"} \
  --narrow_peak

echo "[$(date)] driver exit status: $?"
