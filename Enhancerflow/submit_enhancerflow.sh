#!/bin/bash
# ------------------------------------------------------------------------------
# SLURM driver for Enhancerflow (super-enhancer calling on chipseq output).
#
# Pipeline order:
#   1. activate micromamba env 'dbsuper_pipeline'
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
MEMBER="all"
SKIP_SAMPLESHEET=false
GENOME="mm10"
CHIPSEQ_OUTDIR_OVERRIDE=""           # override chipseq out dir (e.g. ../out/human)
STAGING_DIR_OVERRIDE=""              # override staging dir (e.g. staging_human)
SAMPLESHEET_OVERRIDE=""              # override enhancerflow samplesheet name
SKIP_GREAT=false                     # pass --skip_great (rgreat lacks the mouse TxDb)
OUTDIR="out/"                        # enhancerflow --outdir; override for reruns
CUSTOM_GENOME=""                     # pass --custom_genome (ROSE2 annotation for T2T/CHM13)
GTF=""                               # pass a local --gtf instead of the iGenomes GTF

usage() {
  echo "Usage: sbatch $(basename "$0") [--member <1-5|all>] [--genome KEY] \\"
  echo "         [--chipseq-outdir DIR] [--staging-dir DIR] [--samplesheet NAME] \\"
  echo "         [--outdir DIR] [--custom-genome PATH] [--gtf PATH] \\"
  echo "         [--skip-great] [--skip-samplesheet]"
  echo ""
  echo "  --member            which chipseq output to consume (default: all)"
  echo "                      maps to ../out/member<N> (or ../out/all)"
  echo "  --genome            assembly passed to ROSE2 + rGREAT (default: mm10)"
  echo "  --custom-genome     ROSE2 gene annotation table (refseq.ucsc) for an"
  echo "                      assembly ROSE2 ships no built-in annotation for"
  echo "                      (T2T/CHM13). When set, ROSE2 uses --custom INSTEAD"
  echo "                      of -g, so --genome no longer reaches ROSE2."
  echo "  --gtf               local GTF, overriding the one --genome resolves to"
  echo "                      from iGenomes. Required for CHM13: the iGenomes"
  echo "                      CHM13 GTF is an S3 path that does not exist, and"
  echo "                      nf-schema rejects it at launch."
  echo "  --outdir            enhancerflow --outdir (default: out/). Use a"
  echo "                      separate dir for artefact reruns so the existing"
  echo "                      run's results are left intact."
  echo "  --chipseq-outdir    override chipseq out dir (e.g. ../out/human)"
  echo "  --staging-dir       override staging dir (e.g. staging_human)"
  echo "  --samplesheet       override enhancerflow samplesheet name"
  echo "  --skip-great        pass --skip_great to enhancerflow (use for mouse)"
  echo "  --skip-samplesheet  reuse the existing enhancerflow samplesheet"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --member)            MEMBER="$2";       shift 2 ;;
    --genome)            GENOME="$2";       shift 2 ;;
    --outdir)            OUTDIR="$2";       shift 2 ;;
    --chipseq-outdir)    CHIPSEQ_OUTDIR_OVERRIDE="$2"; shift 2 ;;
    --staging-dir)       STAGING_DIR_OVERRIDE="$2";    shift 2 ;;
    --samplesheet)       SAMPLESHEET_OVERRIDE="$2";    shift 2 ;;
    --custom-genome)     CUSTOM_GENOME="$2";           shift 2 ;;
    --gtf)               GTF="$2";                     shift 2 ;;
    --skip-great)        SKIP_GREAT=true;   shift ;;
    --skip-samplesheet)  SKIP_SAMPLESHEET=true; shift ;;
    -h|--help)           usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# The ROSE2 module interpolates params.custom_genome straight into the command
# line — Nextflow does NOT stage it as a path input. So it has to be an absolute
# path that is visible *inside* the container, which on Karakoram means it must
# live under a root that nextflow.slurm.config binds (`--bind /storage`).
# Resolve it here so a relative path on the sbatch line can't silently break the
# ROSE2 tasks half an hour into the run.
if [[ -n "$CUSTOM_GENOME" ]]; then
  [[ -f "$CUSTOM_GENOME" ]] || {
    echo "ERROR: --custom-genome file not found: $CUSTOM_GENOME" >&2
    exit 2
  }
  CUSTOM_GENOME="$(realpath "$CUSTOM_GENOME")"
fi

# --genome resolves gtf/fasta from iGenomes, and nf-schema validates both with
# `exists: true` BEFORE any task runs. The iGenomes CHM13 entry points its GTF at
# s3://ngi-igenomes/.../Homo_sapiens/NCBI/CHM13/Annotation/Genes/genes.gtf, which
# is not actually in the public bucket — the run aborts at parameter validation
# with "the file or directory ... does not exist". So CHM13 needs a local --gtf,
# the same override submit_driver.sh needs for the chipseq run. (hg38/mm10 GTFs
# do exist on S3, hence no default here: leaving GTF empty preserves the
# existing behaviour for those assemblies.)
if [[ -n "$GTF" ]]; then
  [[ -f "$GTF" ]] || {
    echo "ERROR: --gtf file not found: $GTF" >&2
    exit 2
  }
  GTF="$(realpath "$GTF")"
fi

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

# per-organism overrides (default = above, i.e. ENCODE behaviour unchanged)
[[ -n "$CHIPSEQ_OUTDIR_OVERRIDE" ]] && CHIPSEQ_OUTDIR="$CHIPSEQ_OUTDIR_OVERRIDE"
[[ -n "$SAMPLESHEET_OVERRIDE" ]]    && SAMPLESHEET="$SAMPLESHEET_OVERRIDE"

TS="$(date +%Y%m%d-%H%M%S)"
mkdir -p logs/driver logs/samplesheet logs/nextflow logs/reports "$OUTDIR"

SAMPLESHEET_LOG="logs/samplesheet/samplesheet-${MEMBER}-${TS}.log"

echo "[$(date)] === enhancerflow driver start ==="
echo "  enhancerflow_dir : $ENHANCERFLOW_DIR"
echo "  member           : $MEMBER"
echo "  genome           : $GENOME"
echo "  custom_genome    : ${CUSTOM_GENOME:-(none — ROSE2 uses its built-in -g annotation)}"
echo "  gtf              : ${GTF:-(none — resolved from iGenomes via --genome)}"
echo "  chipseq_outdir   : $CHIPSEQ_OUTDIR"
echo "  samplesheet_out  : $SAMPLESHEET"
echo "  outdir           : $OUTDIR"
echo "  skip-samplesheet : $SKIP_SAMPLESHEET"
echo "  timestamp        : $TS"

# ── 1. activate dbsuper_pipeline environment ──────────────────────────────────
echo "[$(date)] activating micromamba env 'dbsuper_pipeline'"
set +u
export MAMBA_ROOT_PREFIX="/storage/software/micromamba"
eval "$(/usr/local/bin/micromamba shell hook --shell bash)" \
  || { echo "Could not init micromamba; check /usr/local/bin/micromamba"; exit 1; }
micromamba activate dbsuper_pipeline
set -u

# ── 2. make_enhancerflow_samplesheet.py ──────────────────────────────────────
# The samplesheet is now generated by walking the per-sample staging dir
# produced upstream by `sbatch submit_staging.sh`. Each subdir of staging/
# becomes one row.
STAGING_DIR="staging"
[[ -n "$STAGING_DIR_OVERRIDE" ]] && STAGING_DIR="$STAGING_DIR_OVERRIDE"
if $SKIP_SAMPLESHEET; then
  echo "[$(date)] --skip-samplesheet given; reusing existing $SAMPLESHEET"
  [[ -f "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET not found"; exit 1; }
else
  [[ -d "$STAGING_DIR" ]] || {
    echo "ERROR: staging dir not found at $STAGING_DIR — run "
    echo "       'sbatch submit_staging.sh' first to populate it"; exit 1; }

  echo "[$(date)] running make_enhancerflow_samplesheet.py -> $SAMPLESHEET_LOG"
  python3 make_enhancerflow_samplesheet.py \
    --staging-dir "$STAGING_DIR" \
    --condition   H3K27ac \
    --output      "$SAMPLESHEET" \
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
echo "  outdir      : $OUTDIR"

# Track upstream main. `nextflow pull` (run beforehand or via -latest) refreshes
# the local asset cache; resolved SHA is logged below for traceability.
ENHANCERFLOW_REV="main"
echo "  enhancerflow rev: $ENHANCERFLOW_REV (resolving to current main)"

# --genome ($GENOME, default mm10) is the assembly name passed to ROSE2 (-g)
# and (when enabled) rGREAT. --fasta overrides the iGenomes S3 URL (compute
# nodes can't reach it) with the local FASTA that nf-core/chipseq published —
# the same reference the staged BAMs were aligned to, so coordinates match.
#
# --custom_genome (--custom-genome here) is for assemblies ROSE2 has no built-in
# annotation for. ROSE2 only accepts -g {MM8,MM9,MM10,HG18,HG19,HG38}, and its
# argparse puts -g and --custom in a mutually exclusive group, so --custom has to
# REPLACE -g rather than accompany it. Requires khan-lab/enhancerflow >= 612018b
# ("update rose2, to allow custom_genome ..."); on older revisions the ROSE2
# module always emits -g and every ROSE2 task aborts with "invalid choice".
# Skips:
#   --skip_homer  HOMER off (FIMO/SEA motif analysis still runs)
#   --skip_great  GREAT off: the rgreat container ships only the human TxDb,
#                 not TxDb.Mmusculus.UCSC.mm10.knownGene, so GREAT local mode
#                 fails on mouse. Drop this flag once the container bundles the
#                 mouse TxDb (+ org.Mm.eg.db).
GENOME_FASTA="$CHIPSEQ_OUTDIR/genome/genome.fa"
[[ -f "$GENOME_FASTA" ]] || { echo "ERROR: genome FASTA not at $GENOME_FASTA"; exit 1; }

# GREAT: off only when requested (--skip-great). The rgreat container ships just
# the human TxDb, so mouse runs must pass --skip-great; human leaves it on.
SKIP_GREAT_FLAG=""
if $SKIP_GREAT; then SKIP_GREAT_FLAG="--skip_great"; fi
echo "  skip_great      : $SKIP_GREAT"

nextflow -log logs/nextflow/.nextflow.log \
  run khan-lab/enhancerflow -r "$ENHANCERFLOW_REV" \
  -profile apptainer \
  -c "$CONFIG" \
  -resume \
  --input          "$SAMPLESHEET" \
  --genome         "$GENOME" \
  --fasta          "$GENOME_FASTA" \
  --outdir         "$OUTDIR" \
  ${CUSTOM_GENOME:+--custom_genome "$CUSTOM_GENOME"} \
  ${GTF:+--gtf "$GTF"} \
  --skip_crc \
  --skip_comparison \
  --skip_homer \
  ${SKIP_GREAT_FLAG}

echo "[$(date)] driver exit status: $?"
