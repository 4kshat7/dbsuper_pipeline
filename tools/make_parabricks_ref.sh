#!/bin/bash
#
# Build the reference directory used by the hybrid pipeline's Parabricks GPU
# lane (--parabricks_ref). pbrun fq2bam resolves --ref to its real path and
# expects the BWA index next to that real file, named <fasta>.{amb,ann,bwt,pac,sa},
# plus <fasta>.fai. The iGenomes layout keeps FASTA and index in different
# directories (and mm10's index is named genome.* not genome.fa.*), so we copy
# real files into one private directory.
#
# Sources are only READ. Never write into /storage/datasets/igenomes: the BWA
# index directory is an input of every BWA_MEM task, and any new file there
# changes its fingerprint and re-runs all CPU alignments.
#
# Usage:
#   bash tools/make_parabricks_ref.sh GRCh38|mm10 [OUT_DIR]
#   bash tools/make_parabricks_ref.sh --fasta F --bwa-index DIR --out OUT_DIR
# Default OUT_DIR: <dbsuper_pipeline>/ref/parabricks/<genome>  (ref/ is gitignored)

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IG=/storage/datasets/igenomes

FASTA=""; BWA_INDEX=""; OUT=""; NAME=""
case "${1:-}" in
  GRCh38)
    NAME=GRCh38
    FASTA=$IG/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
    BWA_INDEX=$IG/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0
    OUT="${2:-}" ;;
  mm10)
    NAME=mm10
    FASTA=$IG/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
    BWA_INDEX=$IG/Mus_musculus/UCSC/mm10/Sequence/BWAIndex
    OUT="${2:-}" ;;
  *)
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --fasta)     FASTA="$2";     shift 2 ;;
        --bwa-index) BWA_INDEX="$2"; shift 2 ;;
        --out)       OUT="$2";       shift 2 ;;
        -h|--help)   sed -n '2,19p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
      esac
    done
    NAME="custom" ;;
esac
[[ -n "$FASTA" && -n "$BWA_INDEX" ]] || { sed -n '2,19p' "$0"; exit 1; }
OUT="${OUT:-$REPO_DIR/ref/parabricks/$NAME}"

[[ -s "$FASTA" && -s "$FASTA.fai" ]] || { echo "ERROR: need $FASTA and $FASTA.fai" >&2; exit 1; }

# Index prefix = the *.amb file in the index dir (genome.fa.amb or genome.amb)
mapfile -t AMB < <(find -L "$BWA_INDEX" -maxdepth 1 -name '*.amb')
[[ ${#AMB[@]} -eq 1 ]] || { echo "ERROR: expected exactly one *.amb in $BWA_INDEX, found ${#AMB[@]}" >&2; exit 1; }
PREFIX="${AMB[0]%.amb}"
for ext in amb ann bwt pac sa; do
  [[ -s "$PREFIX.$ext" ]] || { echo "ERROR: missing $PREFIX.$ext" >&2; exit 1; }
done

if [[ -e "$OUT" ]]; then
  echo "ERROR: $OUT already exists (it may be an input of a resumable run). Remove it deliberately first." >&2
  exit 1
fi

# The index and the FASTA must describe the same contigs (name, length, order)
python3 - "$PREFIX.ann" "$FASTA.fai" <<'PY'
import sys
ann, fai = sys.argv[1], sys.argv[2]
lines = open(ann).read().split('\n')
n = int(lines[0].split()[1])
ann_contigs = []
for i in range(n):
    name = lines[1 + 2 * i].split()[1]
    length = int(lines[2 + 2 * i].split()[1])
    ann_contigs.append((name, length))
fai_contigs = [(l.split('\t')[0], int(l.split('\t')[1])) for l in open(fai) if l.strip()]
if ann_contigs != fai_contigs:
    sys.exit(f"ERROR: BWA index contigs ({len(ann_contigs)}) differ from FASTA .fai ({len(fai_contigs)})")
print(f"contigs match: {len(ann_contigs)} (index == fasta.fai)")
PY

TMP="$OUT.partial.$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' ERR
echo "Copying into $OUT (reads only from the sources)"
cp "$FASTA"     "$TMP/genome.fa"
cp "$FASTA.fai" "$TMP/genome.fa.fai"
for ext in amb ann bwt pac sa; do
  cp "$PREFIX.$ext" "$TMP/genome.fa.$ext"
done
cmp -s "$FASTA" "$TMP/genome.fa" || { echo "ERROR: FASTA copy differs" >&2; exit 1; }

{
  echo "Parabricks reference built $(date -u +%Y-%m-%dT%H:%M:%SZ) by ${USER:-unknown}"
  echo "fasta:     $FASTA"
  echo "bwa index: $PREFIX.{amb,ann,bwt,pac,sa}"
  ls -la "$TMP"
} > "$TMP/SOURCES.txt"

mv "$TMP" "$OUT"
trap - ERR
# Read-only: nothing (pbrun included) may add or change files here, since the
# directory is a fingerprinted input of every GPU-lane task.
chmod -R a-w "$OUT"
echo "Built: $OUT"
ls -la "$OUT"
