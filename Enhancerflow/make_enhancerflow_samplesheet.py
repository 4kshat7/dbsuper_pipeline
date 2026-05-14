#!/usr/bin/env python3
"""
Build the Enhancerflow samplesheet from a per-sample staging directory.

This is the sample-level rewrite of the old replicate-level generator. It no
longer walks nf-core/chipseq output directly — instead it reads the staging
directory produced by bin/stage_all_samples.py, where each <SAMPLE>/ already
contains the merged-replicate BAMs, sample-level consensus peak BED, and a
merged BigWig.

Why this switched from per-replicate to per-sample: ROSE has no replicate
concept. Whyte 2013 / dbSUPER pool replicate IP reads (and input reads)
before SE calling. Running ROSE per-replicate gave shallow signal and an
unstable elbow on the ranked-enhancer curve.

Input layout:
    <staging_dir>/<SAMPLE>/<SAMPLE>_case.bam        (+ .bai)
    <staging_dir>/<SAMPLE>/<SAMPLE>_control.bam     (+ .bai)
    <staging_dir>/<SAMPLE>/<SAMPLE>_peaks.bed
    <staging_dir>/<SAMPLE>/<SAMPLE>.bigWig          (informational, not in sheet)

Output columns (Enhancerflow schema — unchanged):
    sample,condition,timepoint,bam,peaks,control_bam

Usage:
  python3 Enhancerflow/make_enhancerflow_samplesheet.py \
      --staging-dir Enhancerflow/staging \
      --condition   H3K27ac \
      --absolute
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


ENHANCERFLOW_COLS = ("sample", "condition", "timepoint",
                     "bam", "peaks", "control_bam")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--staging-dir", required=True, type=Path,
                   help="directory containing per-sample subdirs "
                        "(<SAMPLE>/<SAMPLE>_case.bam etc.)")
    p.add_argument("--condition", default="H3K27ac",
                   help="condition / antibody label written into every row "
                        "(default: H3K27ac)")
    p.add_argument("--output", type=Path,
                   default=Path(__file__).resolve().parent / "enhancerflow_samplesheet.csv",
                   help="path to write Enhancerflow samplesheet CSV "
                        "(default: <script_dir>/enhancerflow_samplesheet.csv)")
    p.add_argument("--absolute", action="store_true",
                   help="write absolute paths (default: paths as given, "
                        "i.e. relative if --staging-dir was relative)")
    p.add_argument("--strict", action="store_true",
                   help="fail if any expected file is missing for a sample "
                        "(default: skip the sample with a warning)")
    p.add_argument("--min-peaks", type=int, default=1,
                   help="drop samples whose peaks.bed has fewer than N lines "
                        "(default: 1 — skip empty peak files). ROSE2 fails with "
                        "'xmin not less than xmax' on empty input.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.staging_dir.is_dir():
        sys.exit(f"ERROR: staging dir not found: {args.staging_dir}")

    root = args.staging_dir.resolve() if args.absolute else args.staging_dir
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows_written, rows_skipped = 0, 0
    with args.output.open("w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=ENHANCERFLOW_COLS)
        writer.writeheader()

        for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            sample = sample_dir.name
            bam = sample_dir / f"{sample}_case.bam"
            bai = sample_dir / f"{sample}_case.bam.bai"
            ctl_bam = sample_dir / f"{sample}_control.bam"
            ctl_bai = sample_dir / f"{sample}_control.bam.bai"
            peaks = sample_dir / f"{sample}_peaks.bed"

            missing = [str(p) for p in (bam, bai, ctl_bam, ctl_bai, peaks)
                       if not p.exists()]
            if missing:
                msg = (f"[warn] {sample}: missing files, skipping: "
                       + ", ".join(Path(m).name for m in missing))
                if args.strict:
                    sys.exit("ERROR: " + msg)
                print(msg, file=sys.stderr)
                rows_skipped += 1
                continue

            with peaks.open() as pf:
                peak_count = sum(1 for _ in pf)
            if peak_count < args.min_peaks:
                msg = (f"[warn] {sample}: only {peak_count} peaks "
                       f"(< --min-peaks={args.min_peaks}), skipping")
                if args.strict:
                    sys.exit("ERROR: " + msg)
                print(msg, file=sys.stderr)
                rows_skipped += 1
                continue

            writer.writerow({
                "sample":      sample,
                "condition":   args.condition,
                "timepoint":   "",
                "bam":         str(bam),
                "peaks":       str(peaks),
                "control_bam": str(ctl_bam),
            })
            rows_written += 1

    print(f"[done] wrote {rows_written} rows to {args.output} "
          f"(skipped {rows_skipped})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
