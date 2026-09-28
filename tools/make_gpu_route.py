#!/usr/bin/env python3
"""
Build / extend the CPU-vs-GPU route file for the hybrid nf-core/chipseq
pipeline (--bwa_route). One line per sample-replicate ("<sample>_REP<n>", i.e.
the pipeline's meta.id without its _T<run> suffix):  <key>\t<cpu|gpu>

Rules
  * A key with ANY run already aligned or failed on CPU (BWA_MEM COMPLETED /
    CACHED / FAILED) is pinned to 'cpu': its cached BAMs are kept and a
    multi-run sample never mixes aligners.
  * Existing assignments in --route are never changed (append-only): moving an
    aligned key to the other lane would re-align it and re-run its downstream.
    Exception: --requeue-failed moves 'gpu' keys whose GPU task failed/was
    ignored and that have no GPU alignment back to 'cpu'.
  * Remaining keys are ranked by input size (largest first: the GPU's gain per
    sample grows with size, its ~40 s start-up cost dominates small inputs) and
    given to the GPU until the two lanes are predicted to finish together.

Task status comes from Nextflow trace files and/or `nextflow log` output
(columns: name, status), e.g. for every run of the session:
  for r in $(awk -F'\\t' -v s=<SESSION_UUID> '$6==s{print $3}' .nextflow/history); do
      nextflow log "$r" -f name,status; done > logs/reports/session_task_status.tsv
Keep the old trace files: each launch overwrites logs/reports/nfcore_trace.txt.

Usage
  python3 tools/make_gpu_route.py --samplesheet raw/geo_results/human_chipseq_samplesheet.local.csv \\
      --status logs/reports/nfcore_trace.amazing_elion.txt [--status more.tsv ...] \\
      --route raw/geo_results/bwa_route_human.tsv [--cpu-slots 10] [--layouts SE,PE] [--dry-run]
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

RUN_SUFFIX = re.compile(r"_T\d+$")
DONE = {"COMPLETED", "CACHED"}


def read_status(paths):
    """{(process, meta.id): set(statuses)} from trace TSVs or `nextflow log` output."""
    status = defaultdict(set)
    for path in paths:
        with open(path) as fh:
            first = fh.readline()
            fh.seek(0)
            if "\t" in first and "status" in first.split("\t") and "name" in first.split("\t"):
                rows = ((r["name"], r["status"]) for r in csv.DictReader(fh, delimiter="\t"))
            else:
                rows = (tuple(l.rstrip("\n").split("\t")[:2]) for l in fh if "\t" in l)
            for name, st in rows:
                m = re.match(r"^(\S+) \((.+)\)$", name.strip())
                if not m:
                    continue
                status[(m.group(1).split(":")[-1], m.group(2))].add(st.strip())
    return status


def read_samplesheet(path):
    """{key: {'layout': SE|PE, 'bytes': total fastq bytes, 'runs': n}} (pipeline's ID rules)."""
    keys = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            key = f"{r['sample'].strip().replace(' ', '_')}_REP{int(r['replicate'])}"
            fq = [r["fastq_1"].strip()] + ([r["fastq_2"].strip()] if r.get("fastq_2", "").strip() else [])
            size = sum(os.path.getsize(f) for f in fq if os.path.exists(f))
            k = keys.setdefault(key, {"layout": "PE" if len(fq) == 2 else "SE", "bytes": 0, "runs": 0})
            k["bytes"] += size
            k["runs"] += 1
    return keys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--status", action="append", default=[], help="trace TSV or `nextflow log -f name,status` output (repeatable)")
    ap.add_argument("--route", required=True, help="route TSV to create or extend (append-only)")
    ap.add_argument("--cpu-slots", type=float, default=10, help="concurrent 32-core BWA_MEM tasks the node sustains (default 10)")
    ap.add_argument("--layouts", default="SE,PE", help="layouts allowed on the GPU lane (default SE,PE)")
    # Cost models, seconds per run: fixed + per GB of gzipped FASTQ.
    # CPU: 32-thread BWA_MEM realtime fitted on 1,945 GEO human tasks (2026-09-28);
    # GPU: Parabricks fq2bam wall time from the 2026-09-28 pilot.
    ap.add_argument("--cpu-se", default="71,150", help="CPU SE model 'fixed_s,s_per_gzGB'")
    ap.add_argument("--cpu-pe", default="221,194", help="CPU PE model 'fixed_s,s_per_gzGB'")
    ap.add_argument("--gpu-se", default="35,10", help="GPU SE model 'fixed_s,s_per_gzGB'")
    ap.add_argument("--gpu-pe", default="35,14", help="GPU PE model 'fixed_s,s_per_gzGB'")
    ap.add_argument("--requeue-failed", action="store_true", help="move failed, unaligned 'gpu' keys back to 'cpu'")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    def model(s):
        f, g = (float(x) for x in s.split(","))
        return lambda b, runs: runs * f + g * b / 1e9
    cost = {("cpu", "SE"): model(a.cpu_se), ("cpu", "PE"): model(a.cpu_pe),
            ("gpu", "SE"): model(a.gpu_se), ("gpu", "PE"): model(a.gpu_pe)}
    allowed = set(a.layouts.split(","))

    keys = read_samplesheet(a.samplesheet)
    status = read_status(a.status)

    def key_status(proc):
        out = defaultdict(set)
        for (p, mid), sts in status.items():
            if p == proc:
                out[RUN_SUFFIX.sub("", mid)] |= sts
        return out
    bwa, pb = key_status("BWA_MEM"), key_status("PARABRICKS_FQ2BAM")

    route = {}
    if os.path.exists(a.route):
        for line in open(a.route):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not f[0].startswith("#"):
                route[f[0]] = f[1]
    existing = dict(route)

    requeued = []
    if a.requeue_failed:
        for k, lane in existing.items():
            if lane == "gpu" and (pb.get(k, set()) & {"FAILED"}) and not (pb.get(k, set()) & DONE):
                route[k] = "cpu"
                requeued.append(k)

    pinned = 0
    for k in keys:
        if k not in route and bwa.get(k, set()) & (DONE | {"FAILED"}):
            route[k] = "cpu"
            pinned += 1

    free = [k for k in keys if k not in route]
    # lane loads (seconds) of work that is still to do
    cpu_load = sum(cost[("cpu", keys[k]["layout"])](keys[k]["bytes"], keys[k]["runs"]) for k in free)
    gpu_load = sum(cost[("gpu", keys[k]["layout"])](keys[k]["bytes"], keys[k]["runs"])
                   for k, lane in route.items() if lane == "gpu" and k in keys and not (pb.get(k, set()) & DONE))
    to_gpu = 0
    for k in sorted(free, key=lambda k: keys[k]["bytes"], reverse=True):
        lay = keys[k]["layout"]
        c = cost[("cpu", lay)](keys[k]["bytes"], keys[k]["runs"])
        g = cost[("gpu", lay)](keys[k]["bytes"], keys[k]["runs"])
        if lay in allowed and gpu_load + g <= (cpu_load - c) / a.cpu_slots:
            route[k] = "gpu"
            gpu_load += g
            cpu_load -= c
            to_gpu += 1
    for k in free:
        route.setdefault(k, "cpu")

    n_gpu = sum(1 for k in free if route[k] == "gpu")
    print(f"sample-replicates in samplesheet : {len(keys)}")
    print(f"already in route file            : {len(existing)}  (requeued gpu->cpu: {len(requeued)})")
    print(f"pinned to cpu (aligned/failed)   : {pinned}")
    print(f"remaining, newly routed          : {len(free)}  -> gpu {n_gpu}, cpu {len(free) - n_gpu}")
    print(f"predicted lane time              : CPU {cpu_load / a.cpu_slots / 3600:.1f} h ({a.cpu_slots:g} slots)"
          f"  |  GPU {gpu_load / 3600:.1f} h")
    if a.dry_run:
        print("(dry run: route file not written)")
        return
    tmp = a.route + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("# sample_replicate\tlane   (append-only; see tools/make_gpu_route.py)\n")
        for k in sorted(route):
            fh.write(f"{k}\t{route[k]}\n")
    os.replace(tmp, a.route)
    print(f"wrote {a.route}")


if __name__ == "__main__":
    sys.exit(main())
