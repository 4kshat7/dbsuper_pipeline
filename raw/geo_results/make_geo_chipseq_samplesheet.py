#!/usr/bin/env python3
"""
make_geo_chipseq_samplesheet.py

Build an nf-core/chipseq v2.1.0 samplesheet for the GEO (GSM) to-process set,
from FASTQs already downloaded locally under the geo_data segregation tree.

There is no `encodefetch --nfcore` equivalent for GEO, so this builds the sheet
directly, then reuses the SAME post-processing functions the ENCODE path uses
(`../encode_results/make_nfcode_samplesheet.py`) so both sheets are cleaned
identically:
  handle_mixed_datatypes -> sanitize_pairings -> renumber_replicates -> drop_duplicate_rows

Target schema (7 cols, exact order):
  sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate

Sample naming = accession (the GSM), matching the pipeline convention: the
`sample` value propagates verbatim through staging -> Enhancerflow -> output.

Inputs (per organism):
  - <geo-root>/<org>/to_process.tsv        (segregation output; GSM rows used)
  - <geo-root>/{downloaded,data}/<org>/<GSM>/fastq/SRR*.fastq.gz   (local FASTQ)

A catalog row is usable only if BOTH its case GSM and its control (Input_ID_Name)
GSM have local FASTQ — availability is RE-SCANNED from disk each run (the snapshot
`pair_ready` column is not trusted), so re-running after more downloads picks up
newly-available pairs.

Requires pandas — run under `module load encodefetch/0.5.0` or
`micromamba activate dbsuper_pipeline`.
"""

import argparse
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENCODE_DIR = HERE.parent / "encode_results"
sys.path.insert(0, str(ENCODE_DIR))

try:
    import pandas as pd
    from make_nfcode_samplesheet import (
        handle_mixed_datatypes,
        sanitize_pairings,
        renumber_replicates,
        drop_duplicate_rows,
        convert_to_nfcore_v2_1_0,
    )
except ImportError as e:
    sys.exit(
        f"[ERROR] needs pandas + {ENCODE_DIR}/make_nfcode_samplesheet.py: {e}\n"
        f"        run under: module load encodefetch/0.5.0  "
        f"(or: micromamba activate dbsuper_pipeline)"
    )

GEO_ROOT_DEFAULT = Path("/storage/projects/akshat.mistry/geo_data")
NFCORE_COLS = ["sample", "fastq_1", "fastq_2", "replicate",
               "antibody", "control", "control_replicate"]
ANTIBODY = "H3K27ac"

# SRA run fastq: SRR/ERR/DRR + digits, optional _1/_2 mate, .fastq.gz
SRR_RE = re.compile(r"^((?:SRR|ERR|DRR)\d+)(?:_([12]))?\.fastq\.gz$")


def resolve_fastq_dir(geo_root: Path, org: str, gsm: str) -> Path | None:
    """Return the GSM's fastq/ dir from downloaded/ (symlink) or data/ (real)."""
    for sub in ("downloaded", "data"):
        fq = geo_root / sub / org / gsm / "fastq"
        if fq.is_dir():
            return fq
    return None


def fastq_pairs_for_gsm(fastq_dir: Path) -> list[tuple[str, str]]:
    """One (fastq_1, fastq_2) per SRR run. Single-end -> fastq_2=''.

    Paths are resolved to the real file (realpath) so symlinks in downloaded/
    don't leak into the samplesheet.
    """
    runs: dict[str, dict[str, str]] = {}
    for name in sorted(os.listdir(fastq_dir)):
        m = SRR_RE.match(name)
        if not m:
            continue
        srr, mate = m.group(1), m.group(2)
        path = str((fastq_dir / name).resolve())
        r = runs.setdefault(srr, {})
        if mate == "1":
            r["r1"] = path
        elif mate == "2":
            r["r2"] = path
        else:
            r["single"] = path

    pairs = []
    for srr in sorted(runs):
        r = runs[srr]
        if "r1" in r and "r2" in r:
            pairs.append((r["r1"], r["r2"]))       # paired-end
        elif "r1" in r:
            pairs.append((r["r1"], ""))            # R1 only -> single-end
        elif "single" in r:
            pairs.append((r["single"], ""))        # single-end
    return pairs


def read_to_process_gsm_pairs(to_process: Path):
    """Yield (case_gsm, control_gsm) from GSM rows of to_process.tsv.

    First control wins if a case ever pairs with more than one (mirrors the
    ENCODE builder's `_pick_first`).
    """
    import csv
    seen_case = {}
    order = []
    with to_process.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("accession_namespace") != "GSM":
                continue
            case = (row.get("Case_ID_Name") or "").strip()
            ctrl = (row.get("Input_ID_Name") or "").strip()
            if not case.startswith("GSM") or not ctrl.startswith("GSM"):
                continue
            if case not in seen_case:
                seen_case[case] = ctrl
                order.append(case)
    return [(c, seen_case[c]) for c in order]


def build_rows(geo_root: Path, org: str, pairs) -> tuple[list[dict], dict]:
    """Emit samplesheet rows for pairs whose case AND control FASTQ are on disk."""
    rows = []
    emitted_controls = set()
    stats = {"pairs_total": len(pairs), "pairs_usable": 0,
             "case_missing_fastq": 0, "control_missing_fastq": 0}

    for case, ctrl in pairs:
        case_fq = resolve_fastq_dir(geo_root, org, case)
        ctrl_fq = resolve_fastq_dir(geo_root, org, ctrl)
        if case_fq is None:
            stats["case_missing_fastq"] += 1
            continue
        if ctrl_fq is None:
            stats["control_missing_fastq"] += 1
            continue
        case_pairs = fastq_pairs_for_gsm(case_fq)
        ctrl_pairs = fastq_pairs_for_gsm(ctrl_fq)
        if not case_pairs or not ctrl_pairs:
            # a dir exists but holds no usable SRR fastq -> treat as missing
            if not case_pairs:
                stats["case_missing_fastq"] += 1
            else:
                stats["control_missing_fastq"] += 1
            continue

        stats["pairs_usable"] += 1

        # ChIP (case) rows: one per SRR run, all replicate 1 (nf-core merges them)
        for fq1, fq2 in case_pairs:
            rows.append({"sample": case, "fastq_1": fq1, "fastq_2": fq2,
                         "replicate": 1, "antibody": ANTIBODY,
                         "control": ctrl, "control_replicate": 1})
        # Control rows: emit once per distinct control GSM
        if ctrl not in emitted_controls:
            emitted_controls.add(ctrl)
            for fq1, fq2 in ctrl_pairs:
                rows.append({"sample": ctrl, "fastq_1": fq1, "fastq_2": fq2,
                             "replicate": 1, "antibody": "",
                             "control": "", "control_replicate": ""})
    return rows, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", required=True, choices=["human", "mouse"])
    ap.add_argument("--geo-root", type=Path, default=GEO_ROOT_DEFAULT,
                    help="segregation root holding <org>/, downloaded/, data/")
    ap.add_argument("--to-process", type=Path, default=None,
                    help="override path to <org>/to_process.tsv")
    ap.add_argument("--append-encode", type=Path, default=None,
                    help="an nf-core .local.csv (e.g. the 34 ENCSRs) to fold in "
                         "before sanitizing, producing a combined sheet")
    ap.add_argument("--out", type=Path, default=None,
                    help="output samplesheet (default: <script_dir>/<org>_chipseq_samplesheet.local.csv)")
    args = ap.parse_args()

    org = args.organism
    to_process = args.to_process or (args.geo_root / org / "to_process.tsv")
    out = args.out or (HERE / f"{org}_chipseq_samplesheet.local.csv")
    if not to_process.exists():
        sys.exit(f"[ERROR] to_process.tsv not found: {to_process}")

    pairs = read_to_process_gsm_pairs(to_process)
    rows, stats = build_rows(args.geo_root, org, pairs)
    print(f"[INFO] {org}: {stats['pairs_total']} GSM case/control pairs in catalog; "
          f"{stats['pairs_usable']} usable now "
          f"(case-missing {stats['case_missing_fastq']}, "
          f"control-missing {stats['control_missing_fastq']}).")

    df = pd.DataFrame(rows, columns=NFCORE_COLS) if rows \
        else pd.DataFrame(columns=NFCORE_COLS)

    if args.append_encode:
        if not args.append_encode.exists():
            sys.exit(f"[ERROR] --append-encode file not found: {args.append_encode}")
        enc = pd.read_csv(args.append_encode, dtype=str).fillna("")
        missing = [c for c in NFCORE_COLS if c not in enc.columns]
        if missing:
            sys.exit(f"[ERROR] --append-encode missing columns {missing}")
        enc = enc[NFCORE_COLS]
        print(f"[INFO] appending {len(enc)} ENCODE rows from {args.append_encode}")
        df = pd.concat([df, enc], ignore_index=True)

    if df.empty:
        print("[WARN] no usable rows — writing header-only samplesheet.")
        df[NFCORE_COLS].to_csv(out, index=False)
        return

    # Reuse the ENCODE post-processing verbatim.
    df, _ = handle_mixed_datatypes(df)
    df, _ = sanitize_pairings(df)
    df = renumber_replicates(df)
    df, _ = drop_duplicate_rows(df)
    convert_to_nfcore_v2_1_0(df, out)

    n_chip = int((df["antibody"].astype(str).str.strip() != "").sum())
    n_ctrl = len(df) - n_chip
    print(f"[INFO] wrote {out}: {len(df)} rows "
          f"({n_chip} ChIP + {n_ctrl} control), "
          f"{df['sample'].nunique()} distinct samples.")


if __name__ == "__main__":
    main()
