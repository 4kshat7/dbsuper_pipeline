#!/usr/bin/env python3
"""
Build a URL-based nf-core/chipseq v2.1.0 samplesheet.

Identical sanitization pipeline to make_nfcode_samplesheet.py but keeps the
ENCODE HTTPS URLs in fastq_1 / fastq_2 instead of rewriting them to local
paths. nf-core/chipseq v2.1.0 accepts URLs directly — each Nextflow process
fetches its file into workDir when it runs, so this sidesteps any
pre-download permission issues.

Stages:
  1. Drop the 'single_end' column (not in v2.1.0 schema).
  2. sanitize_pairings    — collapse comma-joined control IDs, drop orphans
                            and unused-control rows.
  3. renumber_replicates  — contiguous 1..N per sample, propagate + clamp
                            control_replicate.
  4. drop_duplicate_rows  — collapse encodefetch's double-emitted control rows.

What is NOT done (compared to the local-path script):
  * rename_misnamed_fastqs   — no local files to rename
  * make_local_samplesheet   — URLs stay as URLs
  * drop_missing_fastqs      — URLs are never 'missing'; nf-core will fetch them

Usage:
  python3 make_nfcode_samplesheet_urls.py \\
      --subset-dir . \\
      --url-samplesheet nfcore_chipseq_samplesheet.csv \\
      --output nfcore_chipseq_samplesheet.urls.csv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

# Reuse the cleanup functions already defined in the local-path script.
from make_nfcode_samplesheet import (
    sanitize_pairings,
    renumber_replicates,
    drop_duplicate_rows,
)


def main():
    ap = argparse.ArgumentParser(
        description="Prepare a URL-based nf-core/chipseq v2.1.0 samplesheet "
                    "(no local-path mapping; URLs kept verbatim)."
    )
    ap.add_argument("--subset-dir", default=".",
                    help="Directory containing the URL samplesheet (default: current dir).")
    ap.add_argument("--url-samplesheet", default="nfcore_chipseq_samplesheet.csv",
                    help="ENCODEfetch-produced URL samplesheet filename "
                         "(default: nfcore_chipseq_samplesheet.csv).")
    ap.add_argument("--output", default="nfcore_chipseq_samplesheet.urls.csv",
                    help="Output filename (default: nfcore_chipseq_samplesheet.urls.csv).")
    args = ap.parse_args()

    subset = Path(args.subset_dir).resolve()
    url_sheet = (subset / args.url_samplesheet).resolve()
    out_path = (subset / args.output).resolve()

    if not url_sheet.exists():
        sys.exit(f"[ERROR] URL samplesheet not found: {url_sheet}")

    print(f"[INFO] Reading URL samplesheet: {url_sheet}")
    df = pd.read_csv(url_sheet)

    # Project to v2.1.0 schema (drop single_end, keep required columns)
    if "single_end" in df.columns:
        df = df.drop(columns=["single_end"])

    required = ["sample", "fastq_1", "fastq_2", "replicate",
                "antibody", "control", "control_replicate"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] Missing required columns for nf-core v2.1.0: {missing}")
    df = df[required].copy()

    n_url = len(df)
    print(f"[INFO] Rows in URL samplesheet: {n_url}")

    # Sanity check: fastq_1 should be all URLs
    bad_fq = df["fastq_1"].dropna().astype(str)
    not_url = bad_fq[~bad_fq.str.startswith(("http://", "https://", "ftp://"))]
    if len(not_url):
        print(f"[WARN] {len(not_url)} fastq_1 values are not URLs — this script is "
              f"intended for URL-based input. Example: {not_url.iloc[0][:80]}")

    # Sanitization pipeline (no drop_missing_fastqs needed for URL mode)
    df, _ = sanitize_pairings(df)
    df = renumber_replicates(df)
    df, _ = drop_duplicate_rows(df)

    df.to_csv(out_path, index=False)

    print()
    print("[SUMMARY]")
    print(f"  Input URL samplesheet : {url_sheet}")
    print(f"  Output URL samplesheet: {out_path}")
    print(f"  Rows                  : {len(df)}  (started with {n_url})")
    chip = df[df["antibody"].notna() & df["antibody"].astype(str).str.strip().ne("")]
    ctrl = df.loc[~df.index.isin(chip.index)]
    print(f"  ChIP rows             : {len(chip)} (samples: {chip['sample'].nunique()})")
    print(f"  Control rows          : {len(ctrl)} (samples: {ctrl['sample'].nunique()})")
    print("  Ready for nf-core/chipseq v2.1.0 — point --input at this CSV; "
          "nf-core will fetch each FASTQ URL into its workDir at runtime.")


if __name__ == "__main__":
    main()
