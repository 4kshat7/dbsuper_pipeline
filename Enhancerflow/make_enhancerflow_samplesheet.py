#!/usr/bin/env python3
"""
Build the Enhancerflow samplesheet from the nf-core/chipseq output directory.

Input:
  - chipseq samplesheet (sample,fastq_1,fastq_2,replicate,antibody,control,
    control_replicate) — the same .local.csv produced by
    raw/encode_results/make_nfcode_samplesheet.py
  - chipseq output directory (the --outdir given to nextflow run)

Output (Enhancerflow expects these columns, in order):
    sample,condition,timepoint,bam,peaks,control_bam

For every ChIP row (antibody non-empty) we emit one Enhancerflow row using:
  sample      = <experiment>_REP<replicate>
  condition   = antibody           (e.g. H3K27ac, BRD4, SMAD3 ...)
  timepoint   = ''                 (left blank — schema requires the column)
  bam         = bwa/merged_library/<sample>_REP<rep>.mLb.clN.sorted.bam
  peaks       = bwa/merged_library/macs3/<narrow|broad>_peak/
                  <sample>_REP<rep>_peaks.<narrowPeak|broadPeak>
  control_bam = bwa/merged_library/<control>_REP<control_replicate>.mLb.clN.sorted.bam

Peak-mode (narrow vs broad) is auto-detected from which subdirectory exists
under macs3/. Pass --peak-mode to force one or the other.

Rows are skipped (with a warning) if any of the three required files is missing
from the chipseq output. Use --strict to fail instead of skipping.

By default the output is written next to this script as
`enhancerflow_samplesheet.csv`. Override with --output if needed.

Usage (from project root):
  python3 Enhancerflow/make_enhancerflow_samplesheet.py \
      --chipseq-samplesheet raw/encode_results/nfcore_chipseq_samplesheet.local.csv \
      --chipseq-outdir      out/member1 \
      --absolute
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


CHIPSEQ_COLS = ("sample", "fastq_1", "fastq_2", "replicate",
                "antibody", "control", "control_replicate")
ENHANCERFLOW_COLS = ("sample", "condition", "timepoint",
                     "bam", "peaks", "control_bam")


def detect_peak_mode(merged_lib: Path) -> str:
    narrow = merged_lib / "macs3" / "narrow_peak"
    broad  = merged_lib / "macs3" / "broad_peak"
    has_narrow = narrow.is_dir() and any(narrow.glob("*_peaks.narrowPeak"))
    has_broad  = broad.is_dir()  and any(broad.glob("*_peaks.broadPeak"))
    if has_narrow and not has_broad:
        return "narrow"
    if has_broad and not has_narrow:
        return "broad"
    if has_narrow and has_broad:
        sys.exit("ERROR: both narrow_peak/ and broad_peak/ found under "
                 f"{merged_lib}/macs3 — pass --peak-mode to choose one")
    sys.exit(f"ERROR: no peak files found under {merged_lib}/macs3 "
             "(expected narrow_peak/ or broad_peak/)")


def peak_paths(peak_mode: str, merged_lib: Path, sample_rep: str) -> Path:
    if peak_mode == "narrow":
        return (merged_lib / "macs3" / "narrow_peak"
                / f"{sample_rep}_peaks.narrowPeak")
    return (merged_lib / "macs3" / "broad_peak"
            / f"{sample_rep}_peaks.broadPeak")


def bam_path(merged_lib: Path, sample_rep: str) -> Path:
    return merged_lib / f"{sample_rep}.mLb.clN.sorted.bam"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chipseq-samplesheet", required=True, type=Path,
                   help="nf-core/chipseq input CSV (.local.csv)")
    p.add_argument("--chipseq-outdir", required=True, type=Path,
                   help="nf-core/chipseq --outdir (contains bwa/merged_library/...)")
    p.add_argument("--output", type=Path,
                   default=Path(__file__).resolve().parent / "enhancerflow_samplesheet.csv",
                   help="path to write Enhancerflow samplesheet CSV "
                        "(default: <script_dir>/enhancerflow_samplesheet.csv)")
    p.add_argument("--peak-mode", choices=("narrow", "broad"), default=None,
                   help="force narrow or broad peaks (default: auto-detect)")
    p.add_argument("--absolute", action="store_true",
                   help="write absolute paths (default: paths as given, "
                        "i.e. relative if --chipseq-outdir was relative)")
    p.add_argument("--strict", action="store_true",
                   help="fail if any expected BAM/peak file is missing "
                        "(default: skip the row with a warning)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.chipseq_samplesheet.is_file():
        sys.exit(f"ERROR: samplesheet not found: {args.chipseq_samplesheet}")
    if not args.chipseq_outdir.is_dir():
        sys.exit(f"ERROR: chipseq outdir not found: {args.chipseq_outdir}")

    out_root = args.chipseq_outdir.resolve() if args.absolute else args.chipseq_outdir
    merged_lib = out_root / "bwa" / "merged_library"
    if not merged_lib.is_dir():
        sys.exit(f"ERROR: expected {merged_lib} to exist (chipseq --outdir layout)")

    peak_mode = args.peak_mode or detect_peak_mode(merged_lib)
    print(f"[info] peak-mode = {peak_mode}", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows_written, rows_skipped = 0, 0
    with args.chipseq_samplesheet.open() as fh, \
         args.output.open("w", newline="") as out_fh:
        reader = csv.DictReader(fh)
        missing = [c for c in CHIPSEQ_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: chipseq samplesheet missing columns: {missing}")

        writer = csv.DictWriter(out_fh, fieldnames=ENHANCERFLOW_COLS)
        writer.writeheader()

        seen: set[str] = set()
        for row in reader:
            antibody = (row["antibody"] or "").strip()
            if not antibody:
                continue  # control rows: no ChIP target, skip

            sample = row["sample"].strip()
            replicate = row["replicate"].strip()
            control = (row["control"] or "").strip()
            control_rep = (row["control_replicate"] or "").strip()
            if not control or not control_rep:
                msg = f"[warn] {sample}_REP{replicate}: missing control mapping"
                if args.strict:
                    sys.exit("ERROR: " + msg)
                print(msg, file=sys.stderr)
                rows_skipped += 1
                continue

            sample_rep  = f"{sample}_REP{replicate}"
            control_rep_id = f"{control}_REP{control_rep}"

            # Multiple chipseq rows can produce the same Enhancerflow row
            # (one per FASTQ). Dedupe on sample_rep.
            if sample_rep in seen:
                continue
            seen.add(sample_rep)

            bam = bam_path(merged_lib, sample_rep)
            peaks = peak_paths(peak_mode, merged_lib, sample_rep)
            ctl_bam = bam_path(merged_lib, control_rep_id)

            missing_files = [str(p) for p in (bam, peaks, ctl_bam) if not p.is_file()]
            if missing_files:
                msg = (f"[warn] {sample_rep}: missing files, skipping: "
                       + ", ".join(missing_files))
                if args.strict:
                    sys.exit("ERROR: " + msg)
                print(msg, file=sys.stderr)
                rows_skipped += 1
                continue

            writer.writerow({
                "sample":      sample_rep,
                "condition":   antibody,
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
