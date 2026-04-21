# nf-core/chipseq Performance Optimization Plan

**For:** dbsuper_pipeline ENCODE H3K27ac cohort (913 experiments / 5 batches)
**Pipeline:** nf-core/chipseq v2.1.0
**Cluster:** CIAI (SLURM, partition `cscc-cpu-p`, QoS `cscc-cpu-qos`)
**Target config:** [../../../member4/dbsuper_pipeline/nfcore_chipseq.config](../../../member4/dbsuper_pipeline/nfcore_chipseq.config)

---

## TL;DR — the three killer lines

Your current config serializes the entire pipeline. Just **three lines** are responsible for the weeks-vs-days difference:

```groovy
process.cpus     = 20          // every task grabs 20 CPUs, even tiny ones
process.maxForks = 1           // only ONE task pipeline-wide at a time
executor.queueSize = 1         // only ONE task queued
```

Remove/replace those + switch `executor` from `'local'` to `'slurm'` and your member4 batch goes from **~3–5 weeks → ~1–2 days**.

GPU is not an option for this pipeline (alignment + peak-calling + ROSE are all CPU-only; see §5).

---

## Table of contents

1. [What the current config is doing wrong](#1-what-the-current-config-is-doing-wrong)
2. [Answers to the four questions](#2-answers-to-the-four-questions)
3. [Step-by-step fix — with code](#3-step-by-step-fix--with-code)
4. [Expected wall-clock speedup](#4-expected-wall-clock-speedup)
5. [Why GPU doesn't help here](#5-why-gpu-doesnt-help-here)
6. [SLURM submission pattern](#6-slurm-submission-pattern)
7. [Rollout / verification plan](#7-rollout--verification-plan)
8. [Troubleshooting the loop-device error](#8-troubleshooting-the-loop-device-error)

---

## 1. What the current config is doing wrong

The config at `member4/dbsuper_pipeline/nfcore_chipseq.config` sets:

```groovy
process {
  executor = 'local'        // runs everything inside ONE SLURM allocation
  cpus    = 20              // blanket default for EVERY process step
  memory  = '64 GB'         // blanket default for EVERY process step
  time    = '72h'
  maxForks = 1              // ← pipeline-wide: only 1 task runs at a time
}
executor { queueSize = 1 }  // ← pipeline-wide: only 1 task queued
```

Four problems:

1. **`maxForks = 1`** — Only one task in the entire DAG executes at a time. A 180-sample batch has ~2,000 tasks (15 steps × 180 samples + joint steps). Serial execution of 2,000 × ~10 min avg = **~333 hours ≈ 14 days** bare-minimum; realistically 3–5 weeks with retries.

2. **`queueSize = 1`** — Same effect, from the executor side.

3. **`executor = 'local'`** — Nextflow runs every task inside the one SLURM node you allocated (cn-17, 24 CPUs). No other node ever gets used, no matter how many you have. All 2,000 tasks funnel through 24 cores.

4. **Blanket `cpus=20, memory=64 GB`** — Overrides nf-core/chipseq's per-process labels. A FastQC that needs 1 CPU + 4 GB holds 20 CPUs for 10 minutes. Even if `maxForks` were higher, this alone blocks packing.

---

## 2. Answers to the four questions

### Q: Is the memory and CPU allocation enough?

**Partially yes, but blanket-sized wrong.** nf-core/chipseq tags each process with a **label** (`process_low`, `process_medium`, `process_high`, etc.) and the base.config sizes each label appropriately. Your config throws that away with a global `cpus=20, memory=64 GB`. Consequences:

| Step | Needs (nf-core default label) | Your config gives | Result |
|---|---|---|---|
| FastQC | process_low: 2 CPU, 12 GB | 20 CPU, 64 GB | 10× wasted CPU |
| Trim Galore | process_high: 12 CPU, 72 GB | 20 CPU, 64 GB | **Memory under-spec by 8 GB** |
| BWA-MEM / Bowtie2 | process_high: 12 CPU, 72 GB | 20 CPU, 64 GB | Memory under-spec; CPU fine |
| Picard MarkDuplicates | process_high_memory: 4 CPU, 200 GB | 20 CPU, 64 GB | **Memory under-spec by 136 GB** — will OOM on deep samples |
| MACS2 callpeak | process_medium: 6 CPU, 36 GB | 20 CPU, 64 GB | Over-spec, blocks packing |
| deepTools plotProfile | process_medium: 6 CPU, 36 GB | 20 CPU, 64 GB | Over-spec |
| Any `_single` step | process_single: 1 CPU, 6 GB | 20 CPU, 64 GB | 20× wasted CPU |

**Fix:** delete the blanket defaults, use `withLabel` overrides. See §3 Step 2.

### Q: Does it parallelize?

**Currently no.** With correct settings nf-core/chipseq parallelizes at **two levels**:

- **Workflow-level (task fan-out):** all 180 FastQCs run at once, all alignments at once, etc. Controlled by `executor.queueSize` and per-process `maxForks`. Currently both set to 1 → no fan-out.
- **Tool-level (multithreading):** BWA/Bowtie2/samtools/deepTools use the `cpus` directive value as their thread count. Already wired up by nf-core.

**Fix:** delete `maxForks=1` and `queueSize=1`, set `queueSize=50`. See §3 Step 3.

### Q: Can we use GPU?

**No.** See §5 for the per-tool breakdown. Short version: every hot tool (BWA, Bowtie2, samtools, Picard, MACS2, deepTools, HOMER, ROSE) is CPU-only. The one tool with experimental GPU support (fastp CUDA) isn't used in the default chipseq 2.1.0 path and accounts for ~5% of runtime anyway.

### Q: Local executor vs SLURM executor — which is faster?

**SLURM, by a huge margin.** Concrete numbers for the member4 batch (~180 samples):

| Configuration | Max parallel tasks | Est. wall time |
|---|---|---|
| Current (`local`, `maxForks=1`, 24 CPU node) | **1** | ~3–5 weeks |
| `local`, remove maxForks, add labels, same node | ~4–8 | ~5–10 days |
| **`slurm`, `queueSize=50`, labels, cluster-wide** | **50+** | **~1–2 days** |
| Same + all 5 batches launched in parallel | 250+ | ~1–2 days for whole 913-sample cohort |

Why SLURM wins: with `executor='slurm'`, Nextflow submits each task as its own `sbatch` job. Alignment tasks run on node A, MACS2 on node B, FastQC on node C, simultaneously. The driver job itself becomes tiny (2 CPU / 8 GB, just the orchestrator).

---

## 3. Step-by-step fix — with code

### Step 1 — Switch executor to SLURM

**In `nfcore_chipseq.config`, change:**

```groovy
// BEFORE
process {
  executor = 'local'
  ...
}
```

**To:**

```groovy
// AFTER
process {
  executor = 'slurm'
  queue = 'cscc-cpu-p'
  clusterOptions = '--qos=cscc-cpu-qos'
  ...
}

executor {
  name = 'slurm'
  queueSize = 50              // max concurrent SLURM child jobs
  submitRateLimit = '10 sec'  // throttle scheduler submissions
  pollInterval = '30 sec'
}
```

**Why it's faster:** Nextflow stops cramming everything into one node and instead sprays tasks across the cluster. `queueSize=50` is a safe starting point — raise it if your `cscc-cpu-qos` allows more concurrent jobs (check with `sacctmgr show qos cscc-cpu-qos format=MaxJobs,MaxSubmit`).

### Step 2 — Replace blanket defaults with nf-core labels

**DELETE these three lines:**

```groovy
cpus   = 20
memory = '64 GB'
time   = '72h'
```

**REPLACE with the label-based block** (right-sized for cn-17-class nodes with 24 CPU / 256 GB RAM — adjust `process_high_memory` if your partition's biggest node has less):

```groovy
process {
  executor = 'slurm'
  queue = 'cscc-cpu-p'
  clusterOptions = '--qos=cscc-cpu-qos'

  // nf-core/chipseq built-in labels — sized per cluster node class
  withLabel:process_single      { cpus = 1;  memory = 6.GB;   time = 4.h  }
  withLabel:process_low         { cpus = 2;  memory = 12.GB;  time = 4.h  }
  withLabel:process_medium      { cpus = 6;  memory = 36.GB;  time = 8.h  }
  withLabel:process_high        { cpus = 12; memory = 72.GB;  time = 16.h }
  withLabel:process_high_memory { cpus = 4;  memory = 200.GB; time = 16.h }
  withLabel:process_long        { time = 20.h }

  errorStrategy = { task.exitStatus in [143,137,104,134,139,140] ? 'retry' : 'finish' }
  maxRetries = 2
  maxErrors  = '-1'
}
```

**Why it's faster:** lightweight steps (FastQC, HOMER annotate) only claim 1–2 CPUs, so dozens can pack into the same cluster capacity while heavy alignments run in parallel elsewhere.

### Step 3 — Remove the two `=1` parallelism killers

**DELETE:**

```groovy
process.maxForks = 1   // line 26 of current config
executor.queueSize = 1 // line 32 of current config
```

Leave `maxForks` **unset** (Nextflow default = unlimited per-process) and set `queueSize = 50` as shown in Step 1. If the loop-device error returns (see §8), re-introduce **per-label** caps only, never pipeline-wide:

```groovy
// fallback ONLY if loop-device errors return
withLabel:process_high { cpus = 12; memory = 72.GB; time = 16.h; maxForks = 8 }
```

### Step 4 — Add nf-core safety caps via `params`

Add these inside the existing `params {}` block so any mis-labelled task can't ask for more than your node can supply:

```groovy
params {
  igenomes_base = "${launchDir}/ref/igenomes"
  macs_gsize    = 2700000000

  // NEW — cluster safety caps
  max_cpus   = 24
  max_memory = '200.GB'
  max_time   = '72.h'
}
```

### Step 5 — Keep these (they're already right)

Don't change the following — they're all correct:

```groovy
workDir = "/scratch/${System.getenv('USER')}/nf_work/nfcore_chipseq"

singularity {
  enabled = true
  autoMounts = true
  cacheDir = "/scratch/${System.getenv('USER')}/singularity_cache"
  pullTimeout = '12h'
}

env.NXF_HOME            = "/scratch/${System.getenv('USER')}/.nextflow"
env.TMPDIR              = "${workDir}/tmp"
env.SINGULARITY_CACHEDIR = "/scratch/${System.getenv('USER')}/singularity_cache"
env.SINGULARITY_TMPDIR   = "${workDir}/singularity_tmp"

trace   { enabled = true; file = "logs/nfcore_trace.txt";    overwrite = true }
report  { enabled = true; file = "logs/nfcore_report.html";  overwrite = true }
timeline{ enabled = true; file = "logs/nfcore_timeline.html"; overwrite = true }
dag     { enabled = true; file = "logs/nfcore_dag.png";       overwrite = true }
```

### Step 6 — Stop pinning the driver to `cn-17`

Your current driver submit has `#SBATCH --nodelist=cn-17`. Once you're on `executor='slurm'`, the driver no longer needs to be on any particular node (it just orchestrates) — and pinning it means if cn-17 is busy you wait. **Remove `--nodelist=cn-17`** and shrink the driver allocation to match its actual workload (Nextflow JVM + job polling = ~2 CPU, 8 GB):

```bash
# OLD
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --nodelist=cn-17

# NEW
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# no --nodelist
```

Full driver script in §6.

---

## 4. Expected wall-clock speedup

For the member4 batch (~180 samples, ~2,000 tasks). Alignment dominates at ~60% of total compute.

| Change layered in | Config state | Parallelism | Est. wall time | vs baseline |
|---|---|---|---|---|
| Baseline (today) | `local`, maxForks=1, 24 CPU node | 1 | 3–5 weeks | 1× |
| + remove maxForks/queueSize | `local`, labels, 24 CPU node | 4–8 | 5–10 days | ~3× |
| + switch to `executor=slurm` | `slurm`, labels, queueSize=50 | 50 | **1–2 days** | **~15×** |
| + launch all 5 batches in parallel drivers | same, × 5 drivers | 250 (cluster-shared) | 1–2 days for full 913 | ~100× for whole cohort |

Why alignment dominates: each H3K27ac ChIP fastq is ~30–80M reads, single-ended, ~500 MB–2 GB gz. BWA-MEM at 12 CPU on a node takes 15–45 min per sample. 180 samples serial = ~60 hours just for alignment; 180 samples at 20-way parallelism = ~3 hours.

---

## 5. Why GPU doesn't help here

| Pipeline step | Tool | Hardware | Why no GPU |
|---|---|---|---|
| QC | FastQC | CPU | Java histogram scans, no GPU port |
| Trimming | Trim Galore | CPU | Perl wrapper around cutadapt — no GPU |
| Trimming (alt) | fastp | CPU (exp. GPU) | CUDA quality scoring exists but nf-core/chipseq 2.1.0 default is Trim Galore; ~5% of runtime anyway |
| Alignment | BWA-MEM / Bowtie2 | **CPU only** | BWT + hash-based string matching — not a GPU-amenable workload. No GPU port exists. **60% of total runtime.** |
| BAM ops | samtools / Picard | CPU (I/O-bound) | Bottleneck is disk, not compute. NVMe scratch helps here, GPU doesn't. |
| Peak calling | MACS2 / MACS3 | CPU | Python statistics (Poisson models). No GPU path. Fast anyway (~5 min/sample). |
| Super-enhancer | ROSE | CPU | Sorting + signal summation. Minutes per sample. Not a bottleneck. |
| BigWig | deepTools bamCoverage | CPU | Multi-threaded (`-p` flag). No GPU. |
| Annotation | HOMER | CPU | String matching against genome DB. Already fast. |

**Conclusion: don't pursue GPU.** Spend effort on parallelism (`executor=slurm`) and I/O (keep workDir on `/scratch`, which you already do). These give orders of magnitude more speedup than any GPU switch would.

---

## 6. SLURM submission pattern

See the companion file [submit_driver.sh.template](submit_driver.sh.template) for a ready-to-use driver script.

Key pattern: **driver submits children, children do the work.**

```
┌──────────────────────────────────────────────────────────┐
│ sbatch submit_driver.sh      (2 CPU / 8 GB / 7 day job)  │
│                                                           │
│   └─ nextflow run nf-core/chipseq ...                    │
│         └─ executor='slurm' submits N child jobs:        │
│                                                           │
│   ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐           │
│   │ FastQC │ │ BWA    │ │ MACS2  │ │ ROSE   │  ... × 50 │
│   │ 1 CPU  │ │ 12 CPU │ │ 6 CPU  │ │ 1 CPU  │           │
│   │ nodeA  │ │ nodeB  │ │ nodeC  │ │ nodeD  │           │
│   └────────┘ └────────┘ └────────┘ └────────┘           │
└──────────────────────────────────────────────────────────┘
```

**Driver script skeleton** (full version in the template file):

```bash
#!/bin/bash
#SBATCH -J nf-chipseq-m4
#SBATCH -p cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=07-00:00:00
#SBATCH -o logs/driver-%j.log
#SBATCH -e logs/driver-%j.err

module load nextflow singularity   # or equivalent activation
cd /scratch/$USER/work/member4/dbsuper_pipeline

nextflow run nf-core/chipseq -r 2.1.0 \
  -profile singularity \
  -c nfcore_chipseq.config \
  -resume \
  --input  raw/encode_results/nfcore_chipseq_samplesheet.final.csv \
  --outdir out/member4 \
  --genome GRCh38 \
  --read_length 50 \
  --narrow_peak
```

---

## 7. Rollout / verification plan

Don't just apply all changes and re-launch 180 samples. Staged rollout:

### Stage A — Dry run (0 cost)

```bash
cd /scratch/$USER/work/member4/dbsuper_pipeline
nextflow run nf-core/chipseq -r 2.1.0 \
  -c nfcore_chipseq.config \
  -profile singularity \
  --input raw/encode_results/nfcore_chipseq_samplesheet.final.csv \
  --outdir out/member4 \
  --genome GRCh38 \
  -preview
```

`-preview` parses the config, resolves labels, and prints the DAG without launching anything. Checks your syntax.

### Stage B — 1-sample smoke test (~1 hour)

Make a mini samplesheet with just the first ChIP+control pair:

```bash
head -n 1 raw/encode_results/nfcore_chipseq_samplesheet.final.csv  > /tmp/smoke.csv
grep -m 2 ENCSR raw/encode_results/nfcore_chipseq_samplesheet.final.csv >> /tmp/smoke.csv
```

Launch as a SLURM driver (see §6). Then in another shell:

```bash
watch -n 5 'squeue -u $USER'
```

You should see the driver job PLUS several child jobs appearing/completing. If you only see the driver and tasks run serially inside it, `executor='slurm'` isn't taking effect — check the config was loaded (`nextflow log`).

### Stage C — Inspect the trace

After the smoke test:

```bash
column -t logs/nfcore_trace.txt | less -S
```

Verify different processes got different `cpus`/`memory` (e.g. FASTQC=1 CPU 6 GB, BWA_MEM=12 CPU 72 GB). If everything is still 20 CPU / 64 GB, your label overrides aren't being applied — check `withLabel:` syntax.

### Stage D — Full member4 batch

```bash
sbatch optimization_plan/submit_driver.sh
```

Monitor with `squeue`, `logs/nfcore_report.html` (open after finish), and `logs/nfcore_timeline.html`.

### Stage E — Launch batches 1, 2, 3, 5 in parallel

Once member4 finishes cleanly, the other 4 batches can launch simultaneously — each with its own driver + its own `workDir` + its own `--outdir`. They won't collide.

```bash
for b in member1 member2 member3 member5; do
  sbatch --job-name=nf-$b optimization_plan/submit_driver.sh "$b"
done
```

(Assumes the driver script accepts a batch name — see the template.)

---

## 8. Troubleshooting the loop-device error

The `maxForks=1` workaround was set after a "Singularity: no loop devices available" error you don't remember precisely. That error happens when **many containers mount `.sif` images simultaneously on one node**, exhausting the kernel's `/dev/loop*` devices (usually 8–32 available).

### Why it probably won't return with `executor='slurm'`

With SLURM executor, containers spread across many nodes. On any single node you'll have at most a handful of concurrent containers — well below the loop-device limit. You already have the other two common causes addressed:

- `singularity.cacheDir` on `/scratch` (not NFS `$HOME`) ✓
- `SINGULARITY_TMPDIR` on local scratch ✓

### If it returns anyway

1. Pre-pull all images before launching, so `.sif` files exist:

```bash
nextflow pull nf-core/chipseq -r 2.1.0
nextflow run nf-core/chipseq -r 2.1.0 -c nfcore_chipseq.config -preview
```

2. Cap concurrency per heavy label only (not pipeline-wide):

```groovy
withLabel:process_high { cpus = 12; memory = 72.GB; maxForks = 8 }
```

3. As a last resort, add runtime options to skip loop mounting:

```groovy
singularity {
  enabled = true
  autoMounts = true
  cacheDir = "/scratch/${System.getenv('USER')}/singularity_cache"
  pullTimeout = '12h'
  runOptions = '--no-mount tmp --cleanenv'
}
```

**Do NOT** revert to `maxForks=1` pipeline-wide. The cost is too high.

---

## Checklist — what you actually need to change

In [nfcore_chipseq.config](../../../member4/dbsuper_pipeline/nfcore_chipseq.config):

- [ ] Change `process.executor = 'local'` → `process.executor = 'slurm'`
- [ ] Add `process.queue = 'cscc-cpu-p'`
- [ ] Add `process.clusterOptions = '--qos=cscc-cpu-qos'`
- [ ] **Delete** `process.cpus = 20`
- [ ] **Delete** `process.memory = '64 GB'`
- [ ] **Delete** `process.maxForks = 1`
- [ ] Add `withLabel:process_*` blocks (Step 2)
- [ ] Change `executor.queueSize = 1` → `executor.queueSize = 50`
- [ ] Add `executor.submitRateLimit = '10 sec'`
- [ ] Add `max_cpus / max_memory / max_time` to `params{}`

For the driver submit script:

- [ ] Drop `--cpus-per-task=24` → `--cpus-per-task=2`
- [ ] Drop `--mem=64G` → `--mem=8G`
- [ ] Remove `--nodelist=cn-17`
- [ ] Add `-resume` to the `nextflow run` command

---

## Companion files in this folder

- [recommended_nfcore_chipseq.config](recommended_nfcore_chipseq.config) — full new config, ready to diff against the current one
- [submit_driver.sh.template](submit_driver.sh.template) — SLURM driver submit script template
