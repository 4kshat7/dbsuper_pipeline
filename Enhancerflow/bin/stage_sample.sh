#!/bin/bash
# ------------------------------------------------------------------------------
# stage_sample.sh — per-sample staging job for Enhancerflow.
#
# Given one biological sample and its per-replicate inputs from nf-core/chipseq,
# produce a tidy <SAMPLE>/ directory that Enhancerflow can ingest as one row:
#
#   <staging_dir>/<SAMPLE>/<SAMPLE>_case.bam        (+ .bai)   merged H3K27ac IP BAM
#   <staging_dir>/<SAMPLE>/<SAMPLE>_control.bam     (+ .bai)   merged input BAM
#   <staging_dir>/<SAMPLE>/<SAMPLE>_peaks.bed                  bedtools-merge of rep narrowPeaks
#   <staging_dir>/<SAMPLE>/<SAMPLE>.bigWig                     CPM-normalised merged-case BigWig
#
# Single-rep inputs are symlinked (saves disk, avoids needless reindexing).
# Multi-rep inputs go through samtools merge / bedtools merge.
#
# Note on the <SAMPLE> name: it can be a "merged" name such as ENCSR597UDW
# that combines reads from BOTH _se (single-end) and _pe (paired-end) ENCODE
# entries for the same biological sample. The orchestrator decides what
# replicate files feed into each <SAMPLE>; this script just merges what it's
# given.
#
# Why merge at all: ROSE has no replicate concept (Whyte 2013 / dbSUPER pooled
# rep BAMs before calling SE). Per-replicate calling gives shallow signal and
# an unstable elbow on the ranked-enhancer curve.
#
# Args (positional, all required):
#   $1  SAMPLE          base sample name (no _REP / _se / _pe suffix)
#   $2  CASE_BAMS       comma-separated absolute paths to per-replicate IP BAMs
#   $3  CONTROL_BAMS    comma-separated absolute paths to per-replicate input BAMs
#   $4  PEAKS           comma-separated absolute paths to per-replicate narrowPeaks
#   $5  STAGING_DIR     absolute path where <SAMPLE>/ will be created
# ------------------------------------------------------------------------------

#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 SAMPLE CASE_BAMS CONTROL_BAMS PEAKS STAGING_DIR" >&2
  echo "  CASE_BAMS / CONTROL_BAMS / PEAKS are comma-separated path lists" >&2
  exit 2
fi

SAMPLE="$1"
IFS=',' read -ra CASE_BAMS    <<< "$2"
IFS=',' read -ra CONTROL_BAMS <<< "$3"
IFS=',' read -ra PEAKS        <<< "$4"
STAGING_DIR="$5"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
OUT="${STAGING_DIR}/${SAMPLE}"
CASE_BAM="${OUT}/${SAMPLE}_case.bam"
CTRL_BAM="${OUT}/${SAMPLE}_control.bam"
PEAKS_BED="${OUT}/${SAMPLE}_peaks.bed"
BIGWIG="${OUT}/${SAMPLE}.bigWig"

# Send stderr to the same place as stdout so a single -o log captures everything.
exec 2>&1

ts() { date '+%F %T'; }

echo "[$(ts)] === stage_sample START — ${SAMPLE} ==="
echo "  case reps     : ${#CASE_BAMS[@]}"
for b in "${CASE_BAMS[@]}"; do echo "                    $b"; done
echo "  control reps  : ${#CONTROL_BAMS[@]}"
for b in "${CONTROL_BAMS[@]}"; do echo "                    $b"; done
echo "  peak reps     : ${#PEAKS[@]}"
for b in "${PEAKS[@]}"; do echo "                    $b"; done
echo "  staging_dir   : $STAGING_DIR"
echo "  threads       : $THREADS"
echo "  output dir    : $OUT"

# ── load tools ────────────────────────────────────────────────────────────────
# Versions pinned to the current Karakoram defaults (from `ml -d av` on
# 2026-05-13). If a default changes upstream we'll notice when this fails
# rather than silently picking up a new version.
set +u
source /etc/profile.d/lmod.sh 2>/dev/null || source /etc/profile.d/modules.sh 2>/dev/null || true
module load samtools/1.23.1 bedtools/2.31.1 deeptools/3.5.6
set -u
echo "  samtools      : $(command -v samtools)  ($(samtools --version | head -1))"
echo "  bedtools      : $(command -v bedtools)  ($(bedtools --version))"
echo "  bamCoverage   : $(command -v bamCoverage)  ($(bamCoverage --version 2>&1))"

mkdir -p "$OUT"

# ── deduplicate control list ──────────────────────────────────────────────────
# Same input control often serves multiple IP reps — merging it twice would
# double-count reads. Sort+uniq before merging.
mapfile -t CONTROL_BAMS_UNIQ < <(printf '%s\n' "${CONTROL_BAMS[@]}" | sort -u)
echo "  control reps (deduped): ${#CONTROL_BAMS_UNIQ[@]}"

# ── helpers ───────────────────────────────────────────────────────────────────
merge_or_link_bam() {
  # $1 = output path, remaining args = input BAM paths
  local out="$1"; shift
  local inputs=("$@")
  if [[ ${#inputs[@]} -eq 1 ]]; then
    echo "  [link] $out -> ${inputs[0]}"
    ln -sfn "${inputs[0]}" "$out"
    if [[ -f "${inputs[0]}.bai" ]]; then
      ln -sfn "${inputs[0]}.bai" "${out}.bai"
    else
      echo "  [index] ${out}.bai (source .bai missing — building)"
      samtools index -@ "$THREADS" "$out"
    fi
  else
    echo "  [merge] $out  <-  ${#inputs[@]} BAMs"
    samtools merge -@ "$THREADS" -f "$out" "${inputs[@]}"
    echo "  [index] ${out}.bai"
    samtools index -@ "$THREADS" "$out"
  fi
}

# ── case BAM ──────────────────────────────────────────────────────────────────
echo "[$(ts)] CASE BAM ..."
merge_or_link_bam "$CASE_BAM" "${CASE_BAMS[@]}"

# ── control BAM ──────────────────────────────────────────────────────────────
echo "[$(ts)] CONTROL BAM ..."
merge_or_link_bam "$CTRL_BAM" "${CONTROL_BAMS_UNIQ[@]}"

# ── per-sample consensus peaks ────────────────────────────────────────────────
# Union (>=1 rep): bedtools merge of all replicate narrowPeak intervals.
# Matches nf-core --min_reps_consensus 1 default and Whyte/dbSUPER convention.
#
# Output: 6-column BED (chrom start end name score strand) — the metadata
# columns are CARRIED OVER from the source narrowPeak files, not synthesised:
#   col 4 (name)   ← `first`: the first source peak's MACS-assigned name
#                    (e.g. ENCSR600TOW_REP1_peak_42) for traceability.
#   col 5 (score)  ← `max`:   the strongest source peak's score across the
#                    merge group, preserving the best MACS q-value-derived
#                    quality among overlapping rep peaks.
#   col 6 (strand) ← `first`: inherited (narrowPeak strand is always `.`).
#
# Why 6 cols (and not 3): ROSE2's bed_to_gff reads line[3] and line[4] when
# converting the BED into the GFF it uses internally. 3-col BEDs throw an
# IndexError. ROSE2 doesn't use the score for SE ranking (it re-derives signal
# from the BAMs) but it must be parseable as a numeric column.
#
# Single-rep samples also go through bedtools merge — MACS3 peaks are usually
# already non-overlapping so the merge is a no-op there, but treating both
# cases identically keeps the output schema uniform.
echo "[$(ts)] PEAKS (bedtools merge, carrying name/score/strand from source) ..."
cat "${PEAKS[@]}" | sort -k1,1 -k2,2n | \
  bedtools merge -c 4,5,6 -o first,max,first -i - > "$PEAKS_BED"
PEAK_COUNT=$(wc -l < "$PEAKS_BED")
echo "  intervals: $PEAK_COUNT  (6-col BED)"

# ── merged-replicate BigWig (CPM) ─────────────────────────────────────────────
echo "[$(ts)] BIGWIG (bamCoverage --normalizeUsing CPM) ..."
bamCoverage \
  --bam "$CASE_BAM" \
  --outFileName "$BIGWIG" \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --binSize 25 \
  --numberOfProcessors "$THREADS"

echo "[$(ts)] === stage_sample DONE — ${SAMPLE} ==="
echo "  outputs:"
echo "    $CASE_BAM ($(stat -c%s "$CASE_BAM" 2>/dev/null || echo '?') bytes)"
echo "    $CTRL_BAM ($(stat -c%s "$CTRL_BAM" 2>/dev/null || echo '?') bytes)"
echo "    $PEAKS_BED ($PEAK_COUNT intervals)"
echo "    $BIGWIG"
