#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

ENCFF_RE = re.compile(r"(ENCFF[0-9A-Z]+)")


def is_gzip_file(path: Path) -> bool:
    """Return True if file looks like gzip (by magic bytes)."""
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def rename_misnamed_fastqs(files_root: Path) -> None:
    """
    Rename *.fastq -> *.fastq.gz if it is actually gzip content.
    Leaves non-gzip *.fastq unchanged (prints a warning).
    """
    fastq_files = list(files_root.rglob("*.fastq"))
    if not fastq_files:
        print("[INFO] No *.fastq files found to rename.")
        return

    renamed = 0
    skipped = 0
    for f in fastq_files:
        if is_gzip_file(f):
            gz = f.with_name(f.name + ".gz")  # ENCFF...fastq.gz
            if gz.exists():
                print(f"[WARN] Target exists, skipping rename: {gz}")
                skipped += 1
                continue
            f.rename(gz)
            renamed += 1
        else:
            print(f"[WARN] Not gzip (or corrupt), leaving as-is: {f}")
            skipped += 1

    print(f"[INFO] Renamed {renamed} file(s) from .fastq -> .fastq.gz; skipped {skipped}.")


def url_to_encff(url: str) -> str:
    m = ENCFF_RE.search(str(url))
    if not m:
        raise ValueError(f"Could not parse ENCFF accession from: {url}")
    return m.group(1)


def _pick_existing_path(p: Path, files_root: Path) -> Optional[Path]:
    """
    Resolve a path that may be absolute or relative to subset dir, and ensure it exists.
    """
    # Absolute path
    if p.is_absolute() and p.exists():
        return p.resolve()

    # Relative path (common in manifests)
    cand = (files_root.parent / p).resolve()
    if cand.exists():
        return cand

    # Sometimes manifest stores path relative to files/
    cand2 = (files_root / p).resolve()
    if cand2.exists():
        return cand2

    return None


def build_encff_index_from_manifest(manifest_tsv: Path, files_root: Path) -> dict:
    """
    Build ENCFF -> local path using manifest.tsv (preferred).
    Works even if filenames don't match *.fastq.gz patterns.
    """
    # dtype=str avoids a numpy 2.4.x dtype-inference segfault on this wide TSV
    df = pd.read_csv(manifest_tsv, sep="\t", dtype=str)

    # Try to guess the accession column
    acc_cols = [c for c in df.columns if c.lower() in ("file_accession", "accession", "encff", "file")]
    if not acc_cols:
        # fallback: any column that contains "encff" in its values
        for c in df.columns:
            if df[c].astype(str).str.contains("ENCFF", regex=False).any():
                acc_cols = [c]
                break
    if not acc_cols:
        raise RuntimeError(f"Could not find an accession column in {manifest_tsv}. Columns: {list(df.columns)}")
    acc_col = acc_cols[0]

    # Try to guess the local path column
    path_cols = [c for c in df.columns if c.lower() in ("path", "local_path", "file_path", "filepath", "downloaded_path")]
    if not path_cols:
        # fallback: look for any column that looks like a path under files/
        for c in df.columns:
            s = df[c].astype(str)
            if s.str.contains("/files/", regex=False).any() or s.str.contains("files/", regex=False).any():
                path_cols = [c]
                break

    index = {}

    if path_cols:
        path_col = path_cols[0]
        for _, row in df.iterrows():
            acc = str(row[acc_col]).strip()
            if not acc or acc == "nan":
                continue
            m = ENCFF_RE.search(acc)
            if m:
                acc = m.group(1)

            raw_path = str(row[path_col]).strip()
            if not raw_path or raw_path == "nan":
                continue

            p = _pick_existing_path(Path(raw_path), files_root)
            if p is not None:
                index.setdefault(acc, p)
    else:
        # If no path column exists, we can still use manifest to get accessions,
        # but we must locate files on disk by searching.
        accs = []
        for v in df[acc_col].astype(str).tolist():
            m = ENCFF_RE.search(v)
            if m:
                accs.append(m.group(1))
        accs = sorted(set(accs))
        for acc in accs:
            hits = list(files_root.rglob(f"*{acc}*"))
            for h in hits:
                if h.is_file():
                    index.setdefault(acc, h.resolve())
                    break

    return index


def build_encff_index_by_scanning(files_root: Path) -> dict:
    """
    Aggressive fallback: scan ALL files under files_root and map any that contain ENCFF in
    the filename or any parent directory name.
    """
    index = {}

    # Fast path: common extensions first (prefer gz)
    for ext in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq", "*.gz"):
        for p in files_root.rglob(ext):
            m = ENCFF_RE.search(str(p))  # search full path, not just name
            if not m:
                continue
            acc = m.group(1)
            index.setdefault(acc, p.resolve())

    # Super-aggressive: any file anywhere that has ENCFF in its path
    if not index:
        for p in files_root.rglob("*"):
            if not p.is_file():
                continue
            m = ENCFF_RE.search(str(p))
            if not m:
                continue
            acc = m.group(1)
            index.setdefault(acc, p.resolve())

    return index


def make_local_samplesheet(url_sheet: Path, files_root: Path, out_local: Path, manifest_tsv: Optional[Path]) -> pd.DataFrame:
    """
    Read ENCODEfetch nf-core samplesheet (URL-based) and replace fastq_1/fastq_2 URLs
    with local paths found under files_root (prefer manifest.tsv mapping).
    """
    df = pd.read_csv(url_sheet)

    encff_index = {}
    if manifest_tsv and manifest_tsv.exists():
        encff_index = build_encff_index_from_manifest(manifest_tsv, files_root)
        print(f"[INFO] Built ENCFF index from manifest.tsv: {len(encff_index)} file(s) indexed.")
    if not encff_index:
        encff_index = build_encff_index_by_scanning(files_root)
        print(f"[INFO] Built ENCFF index by scanning files/: {len(encff_index)} file(s) indexed.")

    if not encff_index:
        raise RuntimeError(f"No ENCFF files found under: {files_root}")

    def map_fastq(v: str) -> str:
        v = str(v).strip()
        if not v or v == "nan":
            return ""
        acc = url_to_encff(v)

        if acc in encff_index:
            return str(encff_index[acc])

        # Last-chance: search on disk for any file containing the accession
        hits = [p for p in files_root.rglob(f"*{acc}*") if p.is_file()]
        if hits:
            encff_index[acc] = hits[0].resolve()
            return str(encff_index[acc])

        # File not found - print warning and return empty string
        print(f"[WARN] Local file not found for {acc} (from {v}), skipping this file")
        return ""

    for col in ("fastq_1", "fastq_2"):
        if col in df.columns:
            df[col] = df[col].fillna("").map(map_fastq)

    df.to_csv(out_local, index=False)
    print(f"[INFO] Wrote local-path samplesheet: {out_local}")
    return df


def convert_to_nfcore_v2_1_0(df: pd.DataFrame, out_v2: Path) -> pd.DataFrame:
    """
    nf-core/chipseq v2.1.0 expects header EXACTLY:
    sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
    """
    df = df.copy()
    if "single_end" in df.columns:
        df = df.drop(columns=["single_end"])

    required = ["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns for nf-core v2.1.0: {missing}")

    df = df[required].copy()
    df.to_csv(out_v2, index=False)
    print(f"[INFO] Wrote nf-core v2.1.0-format samplesheet: {out_v2}")
    return df


def fix_control_replicates(df: pd.DataFrame, out_fixed: Path) -> pd.DataFrame:
    """
    Renumber control sample replicates to 1..N (per control sample),
    then set control_replicate for each case row to match its replicate when possible.
    """
    df = df.copy()

    control_samples = sorted(set(df["control"].dropna().astype(str)) - {""})

    for s in control_samples:
        idx = df.index[df["sample"] == s].tolist()
        for n, i in enumerate(idx, start=1):
            df.loc[i, "replicate"] = n

    avail = df.groupby("sample")["replicate"].apply(lambda x: set(map(int, x.tolist()))).to_dict()

    def fix_crep(row):
        ctrl = str(row.get("control", "")).strip()
        if not ctrl:
            return row.get("control_replicate")
        rep = int(row["replicate"])
        if rep in avail.get(ctrl, set()):
            return rep
        a = avail.get(ctrl, {1})
        return min(a) if a else 1

    df["control_replicate"] = df.apply(fix_crep, axis=1).astype(int)

    df.to_csv(out_fixed, index=False)
    print(f"[INFO] Wrote control-fixed samplesheet: {out_fixed}")
    return df


def sanitize_pairings(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean ChIP/control pairings in an nf-core/chipseq v2.1.0 samplesheet.

    Fixes three classes of problems that ENCODEfetch leaves in the sheet:
      1. Comma-joined control IDs (e.g. 'ENCSR001BSB,ENCSR704GTT' in a single
         'control' cell) are collapsed to the first ID so nf-core can resolve
         it against the 'sample' column.
      2. ChIP rows whose control does not resolve to any sample in the sheet
         (true orphans) are dropped. Leaving them would either fail schema
         validation or run MACS2 without input, producing noisy peaks.
      3. Control-style rows (antibody blank) that no ChIP row references are
         dropped. Leaving them wastes alignment compute with no downstream use.

    Prints a before/after summary. Returns the cleaned dataframe.
    """
    df = df.copy()
    n0 = len(df)

    def _pick_first(v):
        s = str(v).strip()
        if not s or s.lower() == "nan":
            return pd.NA
        return s.split(",")[0].strip() or pd.NA

    had_comma = int(df["control"].astype(str).str.contains(",", na=False).sum())
    df["control"] = df["control"].map(_pick_first)

    def _chip_mask(frame):
        return frame["antibody"].notna() & frame["antibody"].astype(str).str.strip().ne("")

    chip_mask = _chip_mask(df)
    ctrl_mask = ~chip_mask
    ctrl_sample_ids = set(df.loc[ctrl_mask, "sample"].astype(str))

    ctrl_valid = df["control"].notna() & df["control"].astype(str).isin(ctrl_sample_ids)
    drop_chip_mask = chip_mask & ~ctrl_valid
    n_drop_chip_rows = int(drop_chip_mask.sum())
    n_drop_chip_samples = int(df.loc[drop_chip_mask, "sample"].nunique())
    df = df.loc[~drop_chip_mask].copy()

    chip_mask = _chip_mask(df)
    ctrl_mask = ~chip_mask
    referenced = set(df.loc[chip_mask, "control"].dropna().astype(str))
    unused_ctrl_mask = ctrl_mask & ~df["sample"].astype(str).isin(referenced)
    n_drop_ctrl_rows = int(unused_ctrl_mask.sum())
    n_drop_ctrl_samples = int(df.loc[unused_ctrl_mask, "sample"].nunique())
    df = df.loc[~unused_ctrl_mask].copy()

    print("[SANITIZE] Pairing cleanup:")
    print(f"  comma-joined control cells collapsed : {had_comma}")
    print(f"  orphan ChIP rows dropped             : {n_drop_chip_rows} "
          f"({n_drop_chip_samples} unique ChIP samples)")
    print(f"  unused control rows dropped          : {n_drop_ctrl_rows} "
          f"({n_drop_ctrl_samples} unique control samples)")
    print(f"  rows before -> after                 : {n0} -> {len(df)}")

    return df.reset_index(drop=True)


def drop_missing_fastqs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop rows whose FASTQ didn't resolve to a local file.

    make_local_samplesheet() replaces each URL with the local path on disk,
    but when a file wasn't downloaded (permission error, network failure,
    still-in-flight) the URL becomes an empty string. nf-core/chipseq schema
    requires a non-empty fastq_1, so such rows must be dropped.

    Also drops broken paired-end rows: if any row of a sample has a non-empty
    fastq_2 (→ sample is paired-end), then every row of that sample needs its
    fastq_2 populated. A missing fastq_2 in a paired-end sample silently turns
    into single-end mid-sample, which nf-core rejects.
    """
    df = df.copy()
    n0 = len(df)

    def _blank(v):
        if pd.isna(v):
            return True
        s = str(v).strip()
        return s in ("", "nan", "NaN", "None")

    fq1_blank = df["fastq_1"].apply(_blank)
    fq2_blank = df["fastq_2"].apply(_blank)

    # Samples with at least one non-blank fastq_2 are paired-end; every row of
    # such a sample must have both fastqs.
    paired_samples = set(df.loc[~fq2_blank, "sample"].astype(str))
    broken_paired = df["sample"].astype(str).isin(paired_samples) & fq2_blank

    drop_mask = fq1_blank | broken_paired
    dropped = df[drop_mask]
    df = df[~drop_mask].copy()

    n_fq1 = int(fq1_blank.sum())
    n_pe_broken = int((broken_paired & ~fq1_blank).sum())
    affected = dropped["sample"].astype(str).unique()
    kept_samples = set(df["sample"].astype(str))
    fully_dropped = [s for s in affected if s not in kept_samples]

    print("[DROP-MISSING-FASTQ] Rows whose local FASTQ path was not found:")
    print(f"  rows dropped (blank fastq_1)         : {n_fq1}")
    print(f"  rows dropped (broken paired-end)     : {n_pe_broken}")
    print(f"  samples affected                     : {len(affected)}")
    print(f"    fully dropped (no replicates left) : {len(fully_dropped)}")
    print(f"    partially dropped (some reps kept) : {len(affected) - len(fully_dropped)}")
    print(f"  rows before -> after                 : {n0} -> {len(df)}")
    if fully_dropped:
        shown = fully_dropped[:10]
        suffix = "..." if len(fully_dropped) > 10 else ""
        print(f"  fully-dropped samples: {shown}{suffix}")

    return df.reset_index(drop=True)


def renumber_replicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize replicate numbering (replaces the separate
    fix_samplesheet_replicates.py step).

    Pass A — per sample, rewrite 'replicate' so values are contiguous 1..N in
             first-seen order. Propagate the renumbering into 'control_replicate'
             of any row pointing to a renumbered value.

    Pass B — for each ChIP row, make sure control_replicate points to a
             replicate that actually exists for the named control. If not,
             clamp to the control's lowest available replicate.

    Final: format control_replicate as an int (or empty string for
           control-style rows, which legitimately have no control of their own).
    """
    df = df.copy()

    # ---- Pass A: build old -> new replicate map per sample -----------------
    rep_map: dict = {}
    for sample, grp in df.groupby("sample", sort=False):
        seen_order = []
        for v in grp["replicate"].astype(str).tolist():
            if v not in seen_order:
                seen_order.append(v)
        for new_idx, old_rep in enumerate(seen_order, start=1):
            rep_map[(str(sample), old_rep)] = new_idx

    def _new_rep(row):
        return rep_map.get((str(row["sample"]), str(row["replicate"])), row["replicate"])

    def _prop_crep(row):
        ctrl = row.get("control")
        crep = row.get("control_replicate")
        if pd.isna(ctrl) or not str(ctrl).strip() or pd.isna(crep):
            return crep
        # control_replicate may come back as 1.0 from CSV parsing
        try:
            crep_key = str(int(float(crep)))
        except (ValueError, TypeError):
            crep_key = str(crep)
        return rep_map.get((str(ctrl).strip(), crep_key), crep)

    df["replicate"] = df.apply(_new_rep, axis=1).astype(int)
    df["control_replicate"] = df.apply(_prop_crep, axis=1)

    # ---- Pass B: clamp control_replicate to an existing rep of the control -
    avail = (df.groupby("sample")["replicate"]
               .apply(lambda s: set(int(v) for v in s.tolist()))
               .to_dict())
    n_clamped = 0

    def _clamp_crep(row):
        nonlocal n_clamped
        ctrl = row.get("control")
        crep = row.get("control_replicate")
        if pd.isna(ctrl) or not str(ctrl).strip() or pd.isna(crep):
            return crep
        pool = avail.get(str(ctrl).strip(), set())
        if not pool:
            return crep
        try:
            cv = int(float(crep))
        except (ValueError, TypeError):
            return crep
        if cv in pool:
            return cv
        n_clamped += 1
        return min(pool)

    df["control_replicate"] = df.apply(_clamp_crep, axis=1)

    # ---- Final formatting: empty string for control rows, int for ChIP rows
    def _fmt(v):
        if pd.isna(v):
            return ""
        if isinstance(v, str) and v.strip() == "":
            return ""
        try:
            return int(float(v))
        except (ValueError, TypeError):
            return ""
    df["control_replicate"] = df["control_replicate"].apply(_fmt)

    print("[RENUMBER] Replicate normalization:")
    print(f"  samples renumbered               : {len(set(k[0] for k in rep_map.keys()))}")
    print(f"  control_replicate values clamped : {n_clamped}")

    return df


def drop_duplicate_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows nf-core/chipseq's check_samplesheet.py rejects as duplicates.

    encodefetch emits every control experiment twice in the URL samplesheet:
      - once as '<ctrl_exp> referenced by a ChIP at its replicate N'
        (control_replicate = N)
      - once as '<ctrl_exp> standalone control row'
        (control_replicate = blank)

    Both rows have identical sample / fastq_1 / fastq_2 / replicate and blank
    antibody + blank control, differing only in control_replicate. For a
    control-style row (antibody empty) the control / control_replicate fields
    are meaningless, so nf-core treats these two rows as duplicates and aborts.

    Fix: for control rows, normalize control and control_replicate to empty,
    then drop exact duplicates on the key tuple. ChIP rows are deduped on the
    full tuple too (same row literally listed twice = duplicate regardless).
    """
    df = df.copy()
    n0 = len(df)

    ctrl_mask = df["antibody"].isna() | df["antibody"].astype(str).str.strip().eq("")

    # For control rows, clear the meaningless control / control_replicate fields
    # so duplicates collapse.
    df.loc[ctrl_mask, "control"] = ""
    df.loc[ctrl_mask, "control_replicate"] = ""

    key_cols = ["sample", "fastq_1", "fastq_2", "replicate",
                "antibody", "control", "control_replicate"]
    # Stable comparison: cast everything to stripped strings (handles 1 vs 1.0 vs NaN).
    key_df = df[key_cols].astype(str).apply(lambda c: c.str.strip())
    dup_mask = key_df.duplicated(keep="first")

    n_dropped_ctrl = int((dup_mask & ctrl_mask).sum())
    n_dropped_chip = int((dup_mask & ~ctrl_mask).sum())
    df = df.loc[~dup_mask].copy()

    print("[DEDUP] Duplicate-row cleanup:")
    print(f"  duplicate control rows dropped : {n_dropped_ctrl}")
    print(f"  duplicate ChIP rows dropped    : {n_dropped_chip}")
    print(f"  rows before -> after           : {n0} -> {len(df)}")

    return df.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(
        description="Prepare ENCODEfetch subset for nf-core/chipseq: rename gzip FASTQs and generate local nf-core samplesheets."
    )
    ap.add_argument("--subset-dir", default=".", help="Path to encode_subset directory (default: current dir).")
    ap.add_argument("--files-dir", default="files", help="Relative path to files/ under subset-dir (default: files).")
    ap.add_argument("--url-samplesheet", default="nfcore_chipseq_samplesheet.csv",
                    help="ENCODEfetch-produced URL samplesheet filename (default: nfcore_chipseq_samplesheet.csv).")
    ap.add_argument("--manifest", default="manifest.tsv",
                    help="ENCODEfetch manifest.tsv filename (default: manifest.tsv).")
    args = ap.parse_args()

    subset = Path(args.subset_dir).resolve()
    files_root = (subset / args.files_dir).resolve()
    url_sheet = (subset / args.url_samplesheet).resolve()
    manifest_tsv = (subset / args.manifest).resolve()

    if not files_root.exists():
        sys.exit(f"[ERROR] files directory not found: {files_root}")
    if not url_sheet.exists():
        sys.exit(f"[ERROR] URL samplesheet not found: {url_sheet}")

    rename_misnamed_fastqs(files_root)

    # Create temporary local samplesheet with URLs replaced
    out_local = subset / "nfcore_chipseq_samplesheet.local.csv"
    df_local = make_local_samplesheet(url_sheet, files_root, out_local, manifest_tsv if manifest_tsv.exists() else None)

    # Convert to nf-core v2.1.0 format (removes single_end column)
    df_v2 = df_local.copy()
    if "single_end" in df_v2.columns:
        df_v2 = df_v2.drop(columns=["single_end"])
    
    required = ["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"]
    missing = [c for c in required if c not in df_v2.columns]
    if missing:
        raise RuntimeError(f"Missing required columns for nf-core v2.1.0: {missing}")
    df_v2 = df_v2[required].copy()

    # 1) Drop rows whose FASTQ didn't download (blank fastq_1 or broken paired).
    #    Must run BEFORE sanitize_pairings so any control that loses all its
    #    replicates cascades into ChIP-orphan detection below.
    df_v2 = drop_missing_fastqs(df_v2)

    # 2) Clean ChIP/control pairings: comma-joined IDs, orphan ChIPs whose
    #    control is no longer in the sheet, and control rows nothing references.
    df_v2 = sanitize_pairings(df_v2)

    # 3) Renumber replicates 1..N per sample and propagate to control_replicate.
    #    Replaces the separate fix_samplesheet_replicates.py step.
    df_v2 = renumber_replicates(df_v2)

    # 4) Drop the duplicate control rows encodefetch emits (same fastq/sample
    #    listed both as 'standalone control' and 'referenced by ChIP at rep N').
    #    nf-core's check_samplesheet.py rejects these as duplicates.
    df_v2 = drop_duplicate_rows(df_v2)

    # Overwrite .local.csv with the fully cleaned, pipeline-ready samplesheet.
    df_v2.to_csv(out_local, index=False)

    print("\n[SUMMARY]")
    print(f"  Final samplesheet : {out_local}")
    print(f"  Rows              : {len(df_v2)}")
    print(f"  ChIP rows         : {int(df_v2['antibody'].notna().sum())} "
          f"(samples: {df_v2.loc[df_v2['antibody'].notna(), 'sample'].nunique()})")
    print(f"  Control rows      : {int(df_v2['antibody'].isna().sum())} "
          f"(samples: {df_v2.loc[df_v2['antibody'].isna(), 'sample'].nunique()})")
    print("  Ready for nf-core/chipseq v2.1.0 — running fix_samplesheet_replicates.py")
    print("  is no longer necessary; this file is already sanitized + renumbered.")


if __name__ == "__main__":
    main()
