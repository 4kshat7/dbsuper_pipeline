#!/bin/bash
# ------------------------------------------------------------------------------
# filter_staged_peaks.sh — retro-fit the primary-contig filter onto an EXISTING
# staging directory, without re-running stage_sample.sh.
#
# Why this exists: stage_sample.sh now filters non-primary contigs out of the
# peaks BED before ROSE2 sees them (see its header for the full rationale). For
# samples that were already staged under the old behaviour, re-running the whole
# staging step would needlessly re-merge BAMs and regenerate BigWigs — neither
# of which the filter touches. This script rewrites only the peaks BED.
#
# For each <SAMPLE>:
#   <SAMPLE>_peaks.bed             -> preserved as <SAMPLE>_peaks.allcontigs.bed
#   <SAMPLE>_peaks.bed             -> rewritten, primary chromosomes only
#
# IDEMPOTENT: if <SAMPLE>_peaks.allcontigs.bed already exists it is treated as
# the authoritative pre-filter source and is NOT overwritten, so re-running this
# script cannot progressively eat the BED.
#
# Usage:
#   bin/filter_staged_peaks.sh --staging-dir staging
#   bin/filter_staged_peaks.sh --staging-dir staging --samples rerun_samples.txt
#   bin/filter_staged_peaks.sh --staging-dir staging --samples rerun_samples.txt --dry-run
#
# Options:
#   --staging-dir DIR   staging root containing one subdir per sample (required)
#   --samples FILE      restrict to sample names in FILE, one per line;
#                       '#' comments and blank lines ignored. Default: all.
#   --dry-run           report what would change, write nothing
#
# Env overrides:
#   PRIMARY_CONTIG_RE   same semantics as in stage_sample.sh. Keep the two in
#                       sync — the default below is deliberately identical.
# ------------------------------------------------------------------------------

set -euo pipefail

STAGING_DIR=""
SAMPLES_FILE=""
DRY_RUN=false

# MUST stay identical to the default in stage_sample.sh.
PRIMARY_CONTIG_RE="${PRIMARY_CONTIG_RE:-^chr([0-9]|1[0-9]|2[0-2]|X|Y)$}"

usage() {
  echo "Usage: $0 --staging-dir DIR [--samples FILE] [--dry-run]" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --staging-dir) STAGING_DIR="$2"; shift 2 ;;
    --samples)     SAMPLES_FILE="$2"; shift 2 ;;
    --dry-run)     DRY_RUN=true; shift ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "$STAGING_DIR" ]] || usage
[[ -d "$STAGING_DIR" ]] || { echo "ERROR: staging dir not found: $STAGING_DIR" >&2; exit 1; }

# ── build the sample list ─────────────────────────────────────────────────────
declare -a SAMPLES=()
if [[ -n "$SAMPLES_FILE" ]]; then
  [[ -f "$SAMPLES_FILE" ]] || { echo "ERROR: samples file not found: $SAMPLES_FILE" >&2; exit 1; }
  while read -r line; do
    line="${line%%#*}"                       # strip trailing comment
    line="$(printf '%s' "$line" | tr -d '[:space:]')"
    [[ -n "$line" ]] && SAMPLES+=("$line")
  done < "$SAMPLES_FILE"
else
  while IFS= read -r d; do
    SAMPLES+=("$(basename "$d")")
  done < <(find "$STAGING_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
fi

[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no samples to process" >&2; exit 1; }

echo "=== filter_staged_peaks ==="
echo "  staging_dir : $STAGING_DIR"
echo "  samples     : ${#SAMPLES[@]}$([[ -n "$SAMPLES_FILE" ]] && echo " (from $SAMPLES_FILE)")"
echo "  primary re  : $PRIMARY_CONTIG_RE"
echo "  dry_run     : $DRY_RUN"
echo

CHANGED=0; UNCHANGED=0; MISSING=0
printf "%-16s %10s %10s %9s  %s\n" SAMPLE BEFORE AFTER DROPPED CONTIGS

for S in "${SAMPLES[@]}"; do
  BED="${STAGING_DIR}/${S}/${S}_peaks.bed"
  ALL="${STAGING_DIR}/${S}/${S}_peaks.allcontigs.bed"

  if [[ ! -f "$BED" && ! -f "$ALL" ]]; then
    printf "%-16s %10s %10s %9s  %s\n" "$S" - - - "MISSING peaks BED"
    MISSING=$((MISSING + 1))
    continue
  fi

  # Prefer an existing pre-filter snapshot as the source of truth.
  SRC="$ALL"
  [[ -f "$ALL" ]] || SRC="$BED"

  BEFORE=$(wc -l < "$SRC")
  AFTER=$(awk -v re="$PRIMARY_CONTIG_RE" 'BEGIN{FS="\t"} $1 ~ re' "$SRC" | wc -l)
  DROPPED=$(( BEFORE - AFTER ))
  CONTIGS=$(awk -v re="$PRIMARY_CONTIG_RE" 'BEGIN{FS="\t"} $1 !~ re {print $1}' "$SRC" \
            | sort | uniq -c | sort -rn | awk '{printf "%s(%s) ",$2,$1}')

  printf "%-16s %10d %10d %9d  %s\n" "$S" "$BEFORE" "$AFTER" "$DROPPED" "${CONTIGS:-none}"

  if [[ "$DROPPED" -eq 0 && -f "$BED" && ! -f "$ALL" ]]; then
    UNCHANGED=$((UNCHANGED + 1))
    continue
  fi

  if [[ "$AFTER" -eq 0 ]]; then
    echo "  ERROR: $S — filter would empty the BED; refusing. Check contig naming." >&2
    exit 1
  fi

  if ! $DRY_RUN; then
    # Snapshot the pre-filter BED once, then rewrite in place via a temp file so
    # an interrupted run can never leave a truncated peaks.bed behind.
    [[ -f "$ALL" ]] || cp -p "$BED" "$ALL"
    TMP="$(mktemp "${STAGING_DIR}/${S}/.${S}_peaks.XXXXXX")"
    awk -v re="$PRIMARY_CONTIG_RE" 'BEGIN{FS=OFS="\t"} $1 ~ re' "$ALL" > "$TMP"
    mv -f "$TMP" "$BED"
  fi
  CHANGED=$((CHANGED + 1))
done

echo
echo "  changed   : $CHANGED"
echo "  unchanged : $UNCHANGED"
echo "  missing   : $MISSING"
$DRY_RUN && echo "  (dry run — nothing written)"
echo "=== done ==="
