#!/bin/bash
#
# Build the patched nf-core/chipseq 2.1.0 with the CPU+GPU hybrid alignment
# lane (pipelines/chipseq-2.1.0-hybrid.patch) and deploy it to a fixed
# directory OUTSIDE this git working tree, for `submit_driver.sh --pipeline`.
#
# Why a fixed, never-touched directory: when Nextflow runs a pipeline from a
# local path it fingerprints the pipeline's bin/ scripts and assets/ files by
# path + size + modification time (not content). Re-creating, moving or
# editing the deployed copy during a run's lifetime would therefore re-run
# every downstream step that uses those files. Deploy once per session.
#
# Usage:
#   bash tools/install_chipseq_hybrid.sh [DEST] [--force]
#     DEST     default: <dbsuper_pipeline>/../pipelines/chipseq-2.1.0-hybrid
#     --force  replace an existing deployment (only between sessions!)
#   CHIPSEQ_SRC=<git url or path> overrides the source repository
#   (default: the local Nextflow asset clone if present, else GitHub).

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PATCH="$REPO_DIR/pipelines/chipseq-2.1.0-hybrid.patch"
TAG="2.1.0"
TAG_COMMIT="76e2382b6d443db4dc2396e6831d1243256d80b0"

DEST=""
FORCE=false
for arg in "$@"; do
  case "$arg" in
    --force) FORCE=true ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) DEST="$arg" ;;
  esac
done
DEST="${DEST:-$REPO_DIR/../pipelines/chipseq-2.1.0-hybrid}"

LOCAL_ASSET="${NXF_HOME:-$HOME/.nextflow}/assets/nf-core/chipseq"
if [[ -z "${CHIPSEQ_SRC:-}" ]]; then
  if [[ -d "$LOCAL_ASSET/.git" ]]; then CHIPSEQ_SRC="$LOCAL_ASSET"
  else CHIPSEQ_SRC="https://github.com/nf-core/chipseq.git"; fi
fi

[[ -f "$PATCH" ]] || { echo "ERROR: patch not found: $PATCH" >&2; exit 1; }

if [[ -e "$DEST" ]]; then
  if ! $FORCE; then
    echo "ERROR: $DEST already exists. It may be in use by a running/resumable session;" >&2
    echo "       re-deploying would re-run downstream steps. Use --force only between sessions." >&2
    exit 1
  fi
  echo "Removing existing deployment $DEST (--force)"
  rm -rf "$DEST"
fi

mkdir -p "$(dirname "$DEST")"
echo "Cloning nf-core/chipseq $TAG from $CHIPSEQ_SRC"
git clone --quiet "$CHIPSEQ_SRC" "$DEST"
git -C "$DEST" -c advice.detachedHead=false checkout --quiet "$TAG"

HEAD_COMMIT="$(git -C "$DEST" rev-parse HEAD)"
if [[ "$HEAD_COMMIT" != "$TAG_COMMIT" ]]; then
  echo "ERROR: $TAG resolved to $HEAD_COMMIT, expected $TAG_COMMIT" >&2
  exit 1
fi

git -C "$DEST" apply --whitespace=nowarn "$PATCH"
echo "Applied $(basename "$PATCH"):"
git -C "$DEST" status --short | sed 's/^/  /'

cat > "$DEST/DBSUPER_HYBRID_DEPLOYMENT.txt" <<EOF
nf-core/chipseq $TAG ($TAG_COMMIT) + $(basename "$PATCH")
patch sha256: $(sha256sum "$PATCH" | cut -d' ' -f1)
deployed: $(date -u +%Y-%m-%dT%H:%M:%SZ) by ${USER:-unknown}
Do not edit, move or re-create this directory while a session that uses it
may still be resumed.
EOF
echo "Deployed: $DEST"
