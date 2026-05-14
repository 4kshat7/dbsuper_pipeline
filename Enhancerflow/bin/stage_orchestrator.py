#!/usr/bin/env python3
"""
Enhancerflow staging orchestrator — single entry point.

What this does, top to bottom:

  1. Read the nf-core/chipseq samplesheet.
  2. Group its rows into BIOLOGICAL SAMPLES. Two ENCODE entries that share an
     ENCSR id but differ only in their _se / _pe suffix (single-end vs paired-
     end sequencing of the same biological sample) are folded into ONE merge.
     E.g. ENCSR597UDW_pe_REP1, ENCSR597UDW_se_REP1, ENCSR597UDW_se_REP2 all
     merge into a single ENCSR597UDW staging dir.
  3. Resolve each sample's per-replicate BAM and narrowPeak paths under
     <chipseq_outdir>/bwa/merged_library/...
  4. Decide which samples need staging:
       - has staging/<SAMPLE>/<SAMPLE>.bigWig  →  done, skip
       - has stage-<SAMPLE> already in SLURM   →  in-flight, skip
       - else                                   →  submit
  5. Submit one  `sbatch stage_sample.sh ...`  per sample. Throttled so we
     never have more than --max-jobs stage-* jobs in queue at once.
  6. Poll until all jobs we submitted have left the queue.
  7. Generate the Enhancerflow samplesheet from the staging dir.
  8. Print a one-screen summary: done / in-flight / submitted / failed.

Why this design:
  - ROSE has no replicate concept. The Whyte 2013 / dbSUPER convention pools
    replicate IP reads (and input reads) into ONE BAM before SE calling.
    Per-replicate SE calling gives shallow signal and an unstable elbow.
  - Single entry point = single command for the user. The script handles
    queue throttling, idempotent restart, _se/_pe merging, and samplesheet
    generation without manual intervention.

Outputs written:
  - <staging_dir>/<SAMPLE>/<SAMPLE>_case.bam (+.bai)
  - <staging_dir>/<SAMPLE>/<SAMPLE>_control.bam (+.bai)
  - <staging_dir>/<SAMPLE>/<SAMPLE>_peaks.bed
  - <staging_dir>/<SAMPLE>/<SAMPLE>.bigWig
  - <output_samplesheet> with columns: sample, condition, timepoint, bam,
    peaks, control_bam — Enhancerflow ingests it as-is.

Logs:
  - logs/<SAMPLE>.log per sample (combined stdout+stderr — both streams go
    to the same file via sbatch -o)
  - the orchestrator itself writes to whatever it's pointed at by the caller
    (typically logs/staging-driver.log via the sbatch wrapper)
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


CHIPSEQ_COLS = ("sample", "fastq_1", "fastq_2", "replicate",
                "antibody", "control", "control_replicate")
ENHANCERFLOW_COLS = ("sample", "condition", "timepoint",
                     "bam", "peaks", "control_bam")

# Strip trailing _se or _pe to derive the biological-sample base name.
# ENCSR597UDW_se → ENCSR597UDW.  ENCSR593INW → ENCSR593INW (unchanged).
_END_TYPE_RE = re.compile(r"_(se|pe)$")


# ── helpers ────────────────────────────────────────────────────────────────────

def ts() -> str:
    return time.strftime("%F %T")


def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", flush=True)


def base_sample(name: str) -> str:
    """E.g. ENCSR597UDW_se → ENCSR597UDW, ENCSR593INW → ENCSR593INW."""
    return _END_TYPE_RE.sub("", name)


def bam_path(merged_lib: Path, full_sample: str, rep: str) -> Path:
    return merged_lib / f"{full_sample}_REP{rep}.mLb.clN.sorted.bam"


def peak_path(merged_lib: Path, full_sample: str, rep: str, mode: str) -> Path:
    suffix = "narrowPeak" if mode == "narrow" else "broadPeak"
    return (merged_lib / "macs3" / f"{mode}_peak"
            / f"{full_sample}_REP{rep}_peaks.{suffix}")


def stage_jobs_in_queue() -> dict[str, str]:
    """Return {sample: jobid} for all stage-* jobs in this user's queue."""
    user = os.environ.get("USER", "")
    try:
        res = subprocess.run(
            ["squeue", "-u", user, "-h", "-o", "%i %j"],
            capture_output=True, text=True, timeout=30)
    except (subprocess.SubprocessError, FileNotFoundError):
        return {}
    if res.returncode != 0:
        return {}
    out: dict[str, str] = {}
    for line in res.stdout.splitlines():
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        jobid, name = parts
        if name.startswith("stage-"):
            out[name[len("stage-"):]] = jobid
    return out


# ── chipseq samplesheet → biological-sample groups ────────────────────────────

def build_sample_groups(chipseq_csv: Path) -> dict[str, list[tuple]]:
    """
    Return  base_sample → [(full_sample, rep, control, control_rep), ...]

    full_sample is what appears in nf-core BAM filenames (may have _se/_pe).
    base_sample is the biological-sample name we use as the staging dir name.
    Replicate rows are deduplicated; FASTQ-per-replicate multiplicity is dropped.
    """
    seen_rows: set[tuple[str, str]] = set()
    groups: dict[str, list[tuple]] = defaultdict(list)
    with chipseq_csv.open() as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in CHIPSEQ_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: chipseq samplesheet missing columns: {missing}")
        for row in reader:
            antibody = (row["antibody"] or "").strip()
            if not antibody:
                continue  # control-only row — already represented as someone's control
            full_sample = row["sample"].strip()
            rep = row["replicate"].strip()
            ctrl = (row["control"] or "").strip()
            ctrl_rep = (row["control_replicate"] or "").strip()
            if not (full_sample and rep and ctrl and ctrl_rep):
                log(f"[warn] skipping incomplete row: sample={full_sample!r} "
                    f"rep={rep!r} control={ctrl!r}")
                continue
            key = (full_sample, rep)
            if key in seen_rows:
                continue
            seen_rows.add(key)
            groups[base_sample(full_sample)].append(
                (full_sample, rep, ctrl, ctrl_rep))
    return groups


# ── samplesheet generation ────────────────────────────────────────────────────

def write_enhancerflow_samplesheet(
    staging_dir: Path,
    output: Path,
    condition: str,
    min_peaks: int,
    absolute: bool,
) -> tuple[int, int]:
    """Walk staging/<SAMPLE>/ subdirs; write one row per complete sample.
    Returns (rows_written, rows_skipped)."""
    root = staging_dir.resolve() if absolute else staging_dir
    output.parent.mkdir(parents=True, exist_ok=True)
    rows_written, rows_skipped = 0, 0
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ENHANCERFLOW_COLS)
        writer.writeheader()
        for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            sample = sample_dir.name
            bam     = sample_dir / f"{sample}_case.bam"
            bai     = sample_dir / f"{sample}_case.bam.bai"
            ctl     = sample_dir / f"{sample}_control.bam"
            ctl_bai = sample_dir / f"{sample}_control.bam.bai"
            peaks   = sample_dir / f"{sample}_peaks.bed"
            bigwig  = sample_dir / f"{sample}.bigWig"
            missing = [p.name for p in (bam, bai, ctl, ctl_bai, peaks, bigwig)
                       if not p.exists()]
            if missing:
                log(f"[samplesheet skip] {sample}: missing {missing}")
                rows_skipped += 1
                continue
            with peaks.open() as pf:
                peak_count = sum(1 for _ in pf)
            if peak_count < min_peaks:
                log(f"[samplesheet skip] {sample}: only {peak_count} peaks")
                rows_skipped += 1
                continue
            writer.writerow({
                "sample":      sample,
                "condition":   condition,
                "timepoint":   "",
                "bam":         str(bam),
                "peaks":       str(peaks),
                "control_bam": str(ctl),
            })
            rows_written += 1
    return rows_written, rows_skipped


# ── orchestration core ────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chipseq-samplesheet", required=True, type=Path)
    p.add_argument("--chipseq-outdir",      required=True, type=Path)
    p.add_argument("--staging-dir",         required=True, type=Path)
    p.add_argument("--logs-dir",            required=True, type=Path,
                   help="directory for per-sample sbatch logs")
    p.add_argument("--output-samplesheet",  required=True, type=Path,
                   help="path to write the Enhancerflow samplesheet")
    p.add_argument("--stage-script", type=Path,
                   default=Path(__file__).resolve().parent / "stage_sample.sh")
    p.add_argument("--condition", default="H3K27ac",
                   help="ChIP target written into every samplesheet row")
    p.add_argument("--peak-mode", choices=("narrow", "broad"), default="narrow")
    p.add_argument("--max-jobs", type=int, default=50,
                   help="cap on concurrent stage-* jobs in queue (default: 50)")
    p.add_argument("--poll-seconds", type=int, default=30,
                   help="queue-poll interval while waiting (default: 30)")
    p.add_argument("--min-peaks", type=int, default=1)
    p.add_argument("--dry-run", action="store_true",
                   help="print the sbatch command per sample without submitting")
    p.add_argument("--force", action="store_true",
                   help="re-stage even if <SAMPLE>.bigWig already exists")
    p.add_argument("--skip-stage", action="store_true",
                   help="don't submit staging jobs; only (re)generate the "
                        "Enhancerflow samplesheet from the existing staging dir")
    p.add_argument("--no-wait", action="store_true",
                   help="submit jobs and exit immediately; don't wait for them "
                        "to finish (samplesheet won't be generated in this mode)")
    p.add_argument("--clean-stale", action="store_true",
                   help="DELETE staging dirs whose name is not a valid base "
                        "sample (e.g. leftover _se / _pe split dirs from a "
                        "previous run). USE WITH CARE.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.chipseq_samplesheet.is_file():
        sys.exit(f"ERROR: chipseq samplesheet not found: {args.chipseq_samplesheet}")
    if not args.chipseq_outdir.is_dir():
        sys.exit(f"ERROR: chipseq outdir not found: {args.chipseq_outdir}")
    if not args.skip_stage:
        if not args.stage_script.is_file():
            sys.exit(f"ERROR: stage script not found: {args.stage_script}")
        if not args.dry_run and shutil.which("sbatch") is None:
            sys.exit("ERROR: sbatch not found on PATH (run on a SLURM submit host)")

    merged_lib = (args.chipseq_outdir / "bwa" / "merged_library").resolve()
    if not merged_lib.is_dir():
        sys.exit(f"ERROR: expected {merged_lib} (chipseq output layout)")

    args.staging_dir.mkdir(parents=True, exist_ok=True)
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    staging_dir = args.staging_dir.resolve()
    logs_dir = args.logs_dir.resolve()

    # ── group chipseq rows by biological sample (_se/_pe folded) ─────────────
    groups = build_sample_groups(args.chipseq_samplesheet)
    log(f"chipseq samplesheet → {len(groups)} biological samples after "
        f"_se/_pe folding")

    # ── optionally delete stale staging dirs ──────────────────────────────────
    if args.clean_stale:
        expected = set(groups.keys())
        for d in sorted(staging_dir.iterdir()):
            if d.is_dir() and d.name not in expected:
                log(f"[clean-stale] removing {d}")
                shutil.rmtree(d)

    # If we're only regenerating the samplesheet, jump straight there.
    if args.skip_stage:
        rw, rs = write_enhancerflow_samplesheet(
            staging_dir, args.output_samplesheet,
            args.condition, args.min_peaks, absolute=True)
        log(f"[samplesheet] wrote {rw} rows  skipped {rs}  "
            f"→ {args.output_samplesheet}")
        return 0

    # ── decide which samples to submit ────────────────────────────────────────
    in_flight = stage_jobs_in_queue()
    if in_flight:
        log(f"{len(in_flight)} stage-* jobs already in queue; those samples "
            f"will be left alone")

    to_submit: list[tuple[str, list[Path], list[Path], list[Path]]] = []
    n_done = n_inflight = n_missing = n_empty = 0

    for sample in sorted(groups.keys()):
        reps = groups[sample]
        out_dir = staging_dir / sample
        bigwig = out_dir / f"{sample}.bigWig"
        if bigwig.exists() and not args.force:
            n_done += 1
            continue
        if sample in in_flight and not args.force:
            log(f"[skip-inflight] {sample}: stage-{sample} jobid={in_flight[sample]}")
            n_inflight += 1
            continue

        case_bams = [bam_path(merged_lib, f, r)        for f, r, _, _ in reps]
        ctrl_bams = [bam_path(merged_lib, c, cr)       for _, _, c, cr in reps]
        peaks     = [peak_path(merged_lib, f, r, args.peak_mode)
                     for f, r, _, _ in reps]

        missing = [str(p) for p in (*case_bams, *ctrl_bams, *peaks) if not p.is_file()]
        if missing:
            log(f"[skip-missing] {sample}: input missing — e.g. {missing[0]}")
            n_missing += 1
            continue

        total_peak_lines = sum(
            sum(1 for _ in p.open()) for p in peaks)
        if total_peak_lines == 0:
            log(f"[skip-empty] {sample}: all replicate narrowPeaks are empty")
            n_empty += 1
            continue

        to_submit.append((sample, case_bams, ctrl_bams, peaks))

    log(f"plan: done={n_done}  in_flight={n_inflight}  to_submit={len(to_submit)}  "
        f"missing={n_missing}  empty={n_empty}")

    # ── submit, throttled ────────────────────────────────────────────────────
    submitted_ids: list[str] = []
    submitted_samples: list[str] = []
    n_submitted = n_submit_failed = 0

    for sample, case_bams, ctrl_bams, peaks in to_submit:
        log_file = logs_dir / f"{sample}.log"
        # Combine stdout+stderr into a single log file by aiming -o and -e at it.
        cmd = [
            "sbatch",
            f"--job-name=stage-{sample}",
            f"--output={log_file}",
            f"--error={log_file}",
            str(args.stage_script),
            sample,
            ",".join(str(p) for p in case_bams),
            ",".join(str(p) for p in ctrl_bams),
            ",".join(str(p) for p in peaks),
            str(staging_dir),
        ]
        if args.dry_run:
            log(f"[dry-run] {sample}: " + " ".join(cmd))
            n_submitted += 1
            continue

        # throttle: wait until there's room
        if args.max_jobs > 0:
            while True:
                current = stage_jobs_in_queue()
                if len(current) < args.max_jobs:
                    break
                log(f"[throttle] {len(current)} stage-* jobs in queue "
                    f"(cap {args.max_jobs}); sleeping {args.poll_seconds}s")
                time.sleep(args.poll_seconds)

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            log(f"[submit-fail] {sample}: sbatch returned {res.returncode}  "
                f"stderr={res.stderr.strip()}")
            n_submit_failed += 1
            continue
        jobid = res.stdout.strip().split()[-1] if res.stdout.strip() else "?"
        submitted_ids.append(jobid)
        submitted_samples.append(sample)
        n_submitted += 1
        log(f"[submit] {sample:30s} reps={len(case_bams)}  jobid={jobid}  "
            f"-> {log_file.name}")

    log(f"submission complete: submitted={n_submitted}  failed={n_submit_failed}")

    if args.dry_run or args.no_wait:
        log("skipping wait + samplesheet generation (--dry-run or --no-wait)")
        return 0

    # ── wait for everything we submitted ──────────────────────────────────────
    if submitted_ids:
        log(f"waiting for {len(submitted_ids)} jobs to finish "
            f"(poll every {args.poll_seconds}s)")
        last_in_queue = -1
        while True:
            current = stage_jobs_in_queue()
            # Only count the samples WE submitted in this run.
            still = [s for s in submitted_samples if s in current]
            if not still:
                break
            if len(still) != last_in_queue:
                log(f"  ... still running/pending: {len(still)}")
                last_in_queue = len(still)
            time.sleep(args.poll_seconds)
        log("all submitted jobs have left the queue")

    # ── generate samplesheet from final staging state ────────────────────────
    rw, rs = write_enhancerflow_samplesheet(
        staging_dir, args.output_samplesheet,
        args.condition, args.min_peaks, absolute=True)
    log(f"[samplesheet] wrote {rw} rows  skipped {rs}  "
        f"→ {args.output_samplesheet}")

    # Final per-sample status check based on filesystem (catches silent failures).
    n_final_ok = n_final_fail = 0
    for s in submitted_samples:
        if (staging_dir / s / f"{s}.bigWig").exists():
            n_final_ok += 1
        else:
            n_final_fail += 1
            log(f"[final-fail] {s}: no bigWig produced  "
                f"(check {logs_dir / (s + '.log')})")

    log("")
    log("=" * 60)
    log("SUMMARY")
    log(f"  biological samples seen          : {len(groups)}")
    log(f"  already done (skipped)           : {n_done}")
    log(f"  already in flight (skipped)      : {n_inflight}")
    log(f"  input missing (skipped)          : {n_missing}")
    log(f"  empty narrowPeaks (skipped)      : {n_empty}")
    log(f"  submitted this run               : {n_submitted}")
    log(f"    → finished with bigWig          : {n_final_ok}")
    log(f"    → failed (no bigWig produced)   : {n_final_fail}")
    log(f"  Enhancerflow samplesheet rows    : {rw}  → {args.output_samplesheet}")
    log("=" * 60)
    return 0 if n_final_fail == 0 and n_submit_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
