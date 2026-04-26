#!/usr/bin/env python3
"""
Per-member control alignment.

Keeps each member's existing ChIP assignments unchanged. The only allowed
modifications are to the *control* set within each member's list:

  * ADD any matched control referenced by a ChIP in the member's list
    that isn't already there.
  * REMOVE any control currently in the member's list that no ChIP in the
    same member references (orphan).

This eliminates the "wasted orphan-control download" problem the original
hardcoded split caused, while keeping per-member ChIP ownership stable.

Inputs
------
  raw/encode_results/manifest.tsv

  members/member1.txt … member5.txt   (lists of ENCSR IDs, one per line)

  If the member files don't exist yet, this script bootstraps them by
  parsing the original ACCESSIONS_N strings out of fetch_member_accessions.sh.

Outputs
-------
  members/member1.txt … member5.txt   (rewritten with aligned controls)

Use --check to re-validate the existing files without rewriting them.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

NUM_MEMBERS = 5
PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ── manifest helpers ──────────────────────────────────────────────────────────
def load_experiment_index(manifest_tsv: Path) -> pd.DataFrame:
    """One row per experiment_accession. dtype=str dodges numpy 2.4 segfault."""
    df = pd.read_csv(manifest_tsv, sep="\t", dtype=str)
    cols = ["experiment_accession", "is_control", "matched_control_experiments"]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit(f"manifest.tsv missing columns: {missing}")
    df = df[cols].drop_duplicates(subset=["experiment_accession"]).reset_index(drop=True)
    df["is_control"] = df["is_control"].astype(str).str.strip().str.lower().eq("true")
    return df


def parse_controls(cell) -> list[str]:
    if cell is None:
        return []
    s = str(cell).strip()
    if not s or s.lower() == "nan":
        return []
    return [t.strip() for t in s.split(",") if t.strip().startswith("ENCSR")]


def chip_to_controls_map(df: pd.DataFrame) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        if row["is_control"]:
            continue
        out[row["experiment_accession"]] = parse_controls(row["matched_control_experiments"])
    return out


# ── member file i/o ───────────────────────────────────────────────────────────
def members_dir() -> Path:
    return PROJECT_ROOT / "members"


def read_member_files() -> list[list[str]]:
    out_dir = members_dir()
    return [
        [ln.strip() for ln in (out_dir / f"member{i}.txt").read_text().splitlines() if ln.strip()]
        for i in range(1, NUM_MEMBERS + 1)
    ]


def write_member_files(member_lists: list[list[str]]) -> None:
    out_dir = members_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    for i, members in enumerate(member_lists, start=1):
        (out_dir / f"member{i}.txt").write_text("\n".join(members) + "\n")
    print(f"\n  wrote {NUM_MEMBERS} files to {out_dir}/")


def bootstrap_from_bash(fetch_script: Path) -> list[list[str]]:
    """Extract ACCESSIONS_1..5='...' from fetch_member_accessions.sh."""
    text = fetch_script.read_text()
    lists: list[list[str]] = []
    for n in range(1, NUM_MEMBERS + 1):
        m = re.search(rf'^ACCESSIONS_{n}="([^"]+)"', text, flags=re.MULTILINE)
        if not m:
            sys.exit(f"could not find ACCESSIONS_{n} in {fetch_script}")
        lists.append([acc.strip() for acc in m.group(1).split(",") if acc.strip()])
    print(f"bootstrapped 5 member lists from {fetch_script}")
    return lists


# ── alignment ─────────────────────────────────────────────────────────────────
def align_controls(member_chips: set[str],
                   member_controls_original: set[str],
                   chip_ctrl_map: dict[str, list[str]]) -> tuple[set[str], set[str], set[str]]:
    """Return (final_set, added_controls, removed_orphan_controls)."""
    required_controls: set[str] = set()
    for chip in member_chips:
        required_controls.update(chip_ctrl_map.get(chip, []))

    added = required_controls - member_controls_original
    removed = member_controls_original - required_controls
    final = member_chips | required_controls
    return final, added, removed


# ── reporting ─────────────────────────────────────────────────────────────────
def summarise(label: str, member_lists, df, edges_per_member=None):
    is_ctrl = dict(zip(df["experiment_accession"], df["is_control"]))
    print(f"\n[{label}] per-member breakdown:")
    print("    member  total  ChIP  control")
    for i, members in enumerate(member_lists, start=1):
        tot = len(members)
        ctrl = sum(1 for a in members if is_ctrl.get(a, False))
        chip = tot - ctrl
        print(f"      {i:>3}    {tot:>3}   {chip:>3}    {ctrl:>3}")


def check_invariants(member_lists, df, chip_ctrl_map) -> int:
    """Return number of ChIPs whose matched controls aren't in the same member."""
    is_ctrl = dict(zip(df["experiment_accession"], df["is_control"]))
    failures = 0
    orphans = 0
    for members in member_lists:
        s = set(members)
        chips = [a for a in members if not is_ctrl.get(a, False)]
        for chip in chips:
            for ctrl in chip_ctrl_map.get(chip, []):
                if ctrl not in s:
                    failures += 1
        ctrls = [a for a in members if is_ctrl.get(a, False)]
        chip_set = set(chips)
        for ctrl in ctrls:
            referenced_here = any(
                ctrl in chip_ctrl_map.get(c, []) for c in chip_set
            )
            if not referenced_here:
                orphans += 1
    print(f"\n  unmatched ChIP→control links across members : {failures}  (must be 0)")
    print(f"  orphan controls (no ChIP in same member ref): {orphans}  (must be 0)")
    return failures + orphans


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--manifest", type=Path,
                    default=PROJECT_ROOT / "raw/encode_results/manifest.tsv")
    ap.add_argument("--fetch-script", type=Path,
                    default=PROJECT_ROOT / "fetch_member_accessions.sh",
                    help="bash script to bootstrap from when members/ doesn't exist yet")
    ap.add_argument("--check", action="store_true",
                    help="Verify existing memberN.txt files, don't write.")
    args = ap.parse_args()

    if not args.manifest.exists():
        sys.exit(f"manifest not found: {args.manifest}")

    print(f"reading: {args.manifest}")
    df = load_experiment_index(args.manifest)
    is_ctrl = dict(zip(df["experiment_accession"], df["is_control"]))
    chip_ctrl_map = chip_to_controls_map(df)
    print(f"  unique experiments: {len(df)}  (controls={int(df['is_control'].sum())})")
    print(f"  ChIP→control links : {sum(len(v) for v in chip_ctrl_map.values())}")

    if (members_dir() / "member1.txt").exists():
        member_lists_in = read_member_files()
        print(f"loaded existing member files from {members_dir()}/")
    else:
        member_lists_in = bootstrap_from_bash(args.fetch_script)

    summarise("BEFORE", member_lists_in, df)

    if args.check:
        rc = check_invariants(member_lists_in, df, chip_ctrl_map)
        sys.exit(0 if rc == 0 else 1)

    # ── apply control-only alignment, ChIPs untouched ─────────────────────────
    new_lists: list[list[str]] = []
    print("\n[ALIGNMENT] per-member changes:")
    for i, members in enumerate(member_lists_in, start=1):
        chips = {a for a in members if not is_ctrl.get(a, False)}
        ctrls = {a for a in members if is_ctrl.get(a, False)}
        unknown = set(members) - chips - ctrls
        if unknown:
            print(f"  member {i}: WARNING — {len(unknown)} accessions not found in manifest: "
                  f"{sorted(unknown)[:5]}{'...' if len(unknown) > 5 else ''}")

        final, added, removed = align_controls(chips, ctrls, chip_ctrl_map)
        new_lists.append(sorted(final))
        print(f"  member {i}: kept {len(chips)} ChIPs; "
              f"+{len(added)} controls added, -{len(removed)} orphans removed "
              f"(final size {len(final)})")

    summarise("AFTER", new_lists, df)
    check_invariants(new_lists, df, chip_ctrl_map)
    write_member_files(new_lists)


if __name__ == "__main__":
    main()
