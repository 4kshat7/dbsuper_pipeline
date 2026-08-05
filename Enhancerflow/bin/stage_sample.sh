#!/bin/bash
# ------------------------------------------------------------------------------
# stage_sample.sh — per-sample staging job for Enhancerflow.
#
# Given one biological sample and its per-replicate inputs from nf-core/chipseq,
# produce a tidy <SAMPLE>/ directory that Enhancerflow can ingest as one row:
#
#   <staging_dir>/<SAMPLE>/<SAMPLE>_case.bam        (+ .bai)   merged H3K27ac IP BAM
#   <staging_dir>/<SAMPLE>/<SAMPLE>_control.bam     (+ .bai)   merged input BAM
#   <staging_dir>/<SAMPLE>/<SAMPLE>_peaks.bed                  bedtools-merge of rep narrowPeaks,
#                                                              PRIMARY CHROMOSOMES ONLY (see below)
#   <staging_dir>/<SAMPLE>/<SAMPLE>_peaks.allcontigs.bed       pre-filter BED, kept for audit
#   <staging_dir>/<SAMPLE>/<SAMPLE>_contig_readcounts.tsv      samtools idxstats of the case BAM
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
# Contig filtering — why the peaks BED is restricted to primary chromosomes:
#
# ROSE2 derives its super-enhancer cutoff geometrically. It sorts all stitched
# enhancers by background-subtracted signal and finds where a line of slope
#     (max(signal) - min(signal)) / n_stitched
# is tangent to that curve. The slope is therefore set by the SINGLE highest
# signal value, so one artefactual region at rank 1 steepens the line, pushes
# the tangent point right, and RAISES the cutoff — silently discarding real
# super-enhancers on real chromosomes.
#
# Non-primary GRCh38/GRCm38 contigs are exactly such artefacts:
#   chrUn_* / *_random  unplaced + unlocalized scaffolds. Alpha-satellite and
#                       segmental-duplication rich, so they act as multi-mapping
#                       read sinks. chrUn_KI270438v1 and chrUn_KI270467v1 are
#                       the worst offenders in this dataset.
#   chrM                mtDNA is packaged by TFAM into nucleoids, NOT into
#                       nucleosomes — there is no histone H3 on chrM, therefore
#                       no H3K27ac. The signal is pure copy-number background
#                       (100s-1000s of mtDNA copies per cell vs 2 nuclear).
#                       ROSE2 stitched all 16,569 bp of it into one "enhancer"
#                       at rank 1 in 14 human samples.
#   chrEBV              viral decoy; real EBV chromatin in LCLs, but not a
#                       human super-enhancer.
#
# Measured cost of NOT filtering: ENCSR405ESP was run twice by accident. With
# chrM at rank 1 the cutoff was 5358.29 and 194 SEs were called; with chr22 at
# rank 1 the cutoff was 2685.54 and 586 SEs were called. Same sample, same
# pipeline — the artefact cost 392 genuine super-enhancers.
#
# We filter INTERVALS, not READS. chrM stays in the FASTA and in both BAMs, so
# the reads-per-million denominator ROSE2 normalises against is unchanged and
# scores stay comparable with earlier runs. Mitochondrial read fraction is
# still captured, as QC, in <SAMPLE>_contig_readcounts.tsv.
#
# Args (positional, all required):
#   $1  SAMPLE          base sample name (no _REP / _se / _pe suffix)
#   $2  CASE_BAMS       comma-separated absolute paths to per-replicate IP BAMs
#   $3  CONTROL_BAMS    comma-separated absolute paths to per-replicate input BAMs
#   $4  PEAKS           comma-separated absolute paths to per-replicate narrowPeaks
#   $5  STAGING_DIR     absolute path where <SAMPLE>/ will be created
#
# Env overrides:
#   PRIMARY_CONTIG_RE   ERE matched against BED col 1; non-matching intervals
#                       are dropped. Default covers hg38 (chr1-22,X,Y) and
#                       mm10 (chr1-19,X,Y) with one pattern. Set to '.' to
#                       disable filtering entirely.
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
PEAKS_ALL="${OUT}/${SAMPLE}_peaks.allcontigs.bed"
CONTIG_QC="${OUT}/${SAMPLE}_contig_readcounts.tsv"
BIGWIG="${OUT}/${SAMPLE}.bigWig"

# Primary-chromosome allowlist. Deliberately written WITHOUT {n,m} interval
# expressions: Ubuntu's default awk is mawk, whose interval support has been
# inconsistent across versions. This pattern is plain ERE and safe everywhere.
#   hg38  -> chr1..chr22, chrX, chrY
#   mm10  -> chr1..chr19, chrX, chrY
# Rejects: chrM, chrEBV, chrUn_*, *_random, *_alt, *_fix, and the bare Ensembl
# accessions (GL456216.1, JH584304.1, ...) that the mouse reference uses for
# unplaced scaffolds.
PRIMARY_CONTIG_RE="${PRIMARY_CONTIG_RE:-^chr([0-9]|1[0-9]|2[0-2]|X|Y)$}"

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
echo "  primary re    : $PRIMARY_CONTIG_RE"

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

# ── contig read-count QC ──────────────────────────────────────────────────────
# We drop chrM *intervals* from the peak set, but the mitochondrial read
# fraction is a genuine signal-to-noise metric (high mito % => poor chromatin
# prep), so capture it here rather than losing it. idxstats is O(1) on an
# indexed BAM — this costs nothing.
echo "[$(ts)] CONTIG READ COUNTS (samtools idxstats) ..."
samtools idxstats "$CASE_BAM" > "$CONTIG_QC"
MITO_READS=$(awk -F'\t' '$1=="chrM" || $1=="MT" {s+=$3} END {print s+0}' "$CONTIG_QC")
MAPPED_READS=$(awk -F'\t' '$1!="*" {s+=$3} END {print s+0}' "$CONTIG_QC")
if [[ "$MAPPED_READS" -gt 0 ]]; then
  echo "  mito reads    : $MITO_READS / $MAPPED_READS mapped" \
       "($(awk -v m="$MITO_READS" -v t="$MAPPED_READS" 'BEGIN{printf "%.2f", 100*m/t}')%)"
else
  echo "  mito reads    : $MITO_READS / 0 mapped (WARNING: no mapped reads)"
fi

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
  bedtools merge -c 4,5,6 -o first,max,first -i - > "$PEAKS_ALL"
ALL_COUNT=$(wc -l < "$PEAKS_ALL")
echo "  intervals (all contigs): $ALL_COUNT  (6-col BED)"

# ── restrict to primary chromosomes ───────────────────────────────────────────
# This is what ROSE2 actually consumes (-i). See the contig-filtering rationale
# in the header block. Filtering here rather than post-hoc on the ROSE2 tables
# is essential: post-filtering would delete the artefact rows but leave the
# inflated cutoff — and the real super-enhancers it suppressed — unrecovered.
echo "[$(ts)] PEAKS filter -> primary chromosomes only ..."
awk -v re="$PRIMARY_CONTIG_RE" 'BEGIN{FS=OFS="\t"} $1 ~ re' "$PEAKS_ALL" > "$PEAKS_BED"
PEAK_COUNT=$(wc -l < "$PEAKS_BED")
DROPPED=$(( ALL_COUNT - PEAK_COUNT ))
echo "  intervals kept   : $PEAK_COUNT"
echo "  intervals dropped: $DROPPED"
if [[ "$DROPPED" -gt 0 ]]; then
  echo "  dropped by contig:"
  awk -v re="$PRIMARY_CONTIG_RE" 'BEGIN{FS="\t"} $1 !~ re {print $1}' "$PEAKS_ALL" \
    | sort | uniq -c | sort -rn | sed 's/^/    /'
fi
if [[ "$PEAK_COUNT" -eq 0 ]]; then
  echo "ERROR: no intervals survived the primary-contig filter." >&2
  echo "       PRIMARY_CONTIG_RE='$PRIMARY_CONTIG_RE' matched nothing in $PEAKS_ALL." >&2
  echo "       Check the reference's chromosome naming (Ensembl '1' vs UCSC 'chr1')." >&2
  exit 1
fi

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
echo "    $PEAKS_BED ($PEAK_COUNT intervals, primary chromosomes only)"
echo "    $PEAKS_ALL ($ALL_COUNT intervals, pre-filter)"
echo "    $CONTIG_QC (mito: $MITO_READS / $MAPPED_READS mapped reads)"
echo "    $BIGWIG"
