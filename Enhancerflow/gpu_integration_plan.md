# Plan: run the CPUs and the 1 GPU together in dbsuper_pipeline (nf-core/chipseq)

## Context

The user asked whether Karakoram's single GPU can speed up the Phase-1 nf-core/chipseq runs. This was prompted by
benchmark: Parabricks `fq2bam` was 5.9–7.2× faster than 32-core BWA-MEM, with identical peaks.

They want to know:
1. Is running CPU + GPU together faster?
2. How do we implement it in `dbsuper_pipeline`?

The work goes on branch **`slrum_karakoram`** (spelled this way locally and on origin). It is pushed only after a
CPU-vs-GPU similarity test passes.

Remaining human work: 4,456 alignments. Human downloads are essentially complete (6,031 of 6,158 GSMs in the sheet). Mouse: 134 done; ~3.2k more are still to download.

User decisions (2026-09-25):
- Mixing aligners is OK: keep the 2,013 finished CPU BWA alignments.
- The GPU is ours to use freely.
- The human run stays stopped until the mouse Enhancerflow job finishes.
- Don't touch the uncommitted `Enhancerflow/*` edits.

Research came from 7 read-only agents, with adversarial checks on the "worth it" and "safe" claims. The numbers below
are the corrected ones.

**What Parabricks is, in brief.** It is NVIDIA's GPU re-implementation of BWA-MEM. Thousands of GPU cores align reads
in parallel, and it produces the same alignments. It is free, runs in a container, and only speeds up the alignment step.

## Answer

**Possible: yes.**
- Parabricks 4.7.1-1 has already run successfully on this exact node, driver 615.71 and 48 GB RTX 6000 Ada, in runs on 2026-09-23.
- The SIF is at `/storage/software/containers/nvcr.io-nvidia-clara-clara-parabricks-4.7.1-1.img`.

**Faster together: yes, as an extra lane next to the CPUs. Never as GPU-only.**
- One GPU handles one sample at a time.
- On our libraries it is ≈3.5 CPU BWA slots (range 2.6–4.9). Per-sample speedup is SE ≈2.7–3.9×, PE ≈4.2–6.5×, and small libraries ≈1×, because of ~50 s fixed start-up.

Remaining-human **alignment phase**. It ends at a hard barrier (`groupTuple`) before any downstream step, so savings here
count 1:1 toward the total:

| setup | alignment phase |
|---|---|
| today's config (≈4.8 BWA slots actually used) | ~150 h |
| today's config + GPU lane | ~70–86 h (saves ~3 days) |
| right-sized requests (≈9–11 slots, SMT-limited) | ~66–82 h |
| **right-sized + GPU lane** | **~45–60 h (saves ~1 day more)** |
| everything on the GPU | 148–277 h, ~560 h with per-task queue waits (slower) |

**The biggest lever is not the GPU.**
- The last run reserved 306/384 cores and 1.2 TB but used ~118 cores, and was memory-blocked ≥35% of the time.
- TrimGalore reserves 32 c / 128 GB and uses ≤6 c / 6 GB; MACS3 and MarkDup use ~1 of 16 cores.
- After the barrier, ~125k downstream tasks run. The GPU can't help with them; they are limited by `queueSize=30` and a
  `submitRateLimit` of 1 job/s (≥35 h of submissions alone).
- Nextflow 25.10.4 does not hash cpus, memory or executor settings (verified in its bytecode), so fixing all of this keeps the `-resume` cache.

## Design

**What we reuse from ** (`/storage/projects/yusuf.maulana/parabricks-benchmark/pipeline-checkouts-backup/chipseq-parabricks`, commit `52bd58a`):
- His fq2bam module and flags: `--no-markdups --filter-flag 256 --bwa-options='-M'`, the read group after the FASTQs, and `--in-se-fq` for SE.
  Together these equal chipseq's `bwa mem -M -R | samtools view -F 0x100`.
- **Not** his `--aligner parabricks` switch, because it:
  - orphans 2,780 finished tasks (2,013 BWA + 627 sort + 140 index);
  - publishes to `…/parabricks/`, which breaks Enhancerflow;
  - uses a symlinked `--ref` that fails on our layout (exit 255).

**Hybrid split inside `FASTQ_ALIGN_BWA`.** `aligner` stays `bwa`, so outputs stay in `out/<org>/bwa/…` and Enhancerflow is unchanged:
```
TRIMGALORE reads ── branch on route[meta.id minus _T<n>] ──┬─ cpu → BWA_MEM (untouched → 2,013 stay cached)
                                                            └─ gpu → PARABRICKS_FQ2BAM (→ <meta.id>.Lb.bam)
                        mix → existing BAM_SORT_STATS_SAMTOOLS → merge → MarkDup → … (all unchanged)
```
- **Routing is per sample-replicate**, so a multi-run GSM never mixes aligners. Keys missing from the route file default to `cpu`.
- The route file is **append-only**. The one exception: a GPU key that failed and has no alignment can be moved to CPU.
- **The GPU lane runs with Nextflow's `local` executor, inside the driver job.** The driver is submitted with
  `sbatch -p gpu --gres=gpu:1 -c 36 --mem=128G`.
  - This removes the 4–12 min SLURM wait per GPU task seen in runs, which would otherwise erase the gain.
  - It holds the GPU for the run; OK per the user.
  - After the alignment barrier the driver can be restarted without `--gres`: everything is cached, and the GPU is freed.

## Implementation (repo `dbsuper_pipeline`, branch `slrum_karakoram`)

**Commit 1 — right-size resources** in `nfcore_chipseq.config`. It helps even without the GPU; values come from trace p99 RSS / p90 %cpu:
- TRIMGALORE: 12 c / 16 GB. Still `--cores 8`, so output is unchanged.
- FASTQC: 4 c / 24 GB
- SAMTOOLS_SORT: 8 c / 24 GB
- **BWA_MEM: keep 32 c** (`-t` sets the PE batch size), memory → 64 GB
- PICARD_MARKDUPLICATES / PICARD_MERGESAMFILES: 4 c / 48 GB. Their Java heap is 0.8 × memory, so keep it generous.
- MACS3_CALLPEAK: 2 c / 16 GB
- BAMTOOLS_FILTER, BEDTOOLS_GENOMECOV: 4 c / 32 GB
- DEEPTOOLS_COMPUTEMATRIX: 16 c / 16 GB
- DEEPTOOLS_PLOTFINGERPRINT: 12 c / 16 GB
- `queueSize`: 30 → 120
- `submitRateLimit`: `'1 sec'` → `'5 sec'`

**Commit 2 — pipeline patch plus installer:**
- `pipelines/chipseq-2.1.0-hybrid.patch`, against nf-core/chipseq tag 2.1.0 (`76e2382`).
- `tools/install_chipseq_hybrid.sh`: clones 2.1.0 from `~/.nextflow/assets/nf-core/chipseq`, applies the patch, and deploys it
  **once** to a stable directory outside the git working tree, `geo_data/pipelines/chipseq-2.1.0-hybrid/`.
  It refuses to overwrite without `--force`.

  A stable location matters because a local-path launch hashes `bin/` and `assets/` files by path + mtime. Touching that
  directory mid-session would re-run the downstream steps.

Patch contents:
- `modules/local/parabricks_fq2bam/main.nf`, adapted from module:
  - input: reads + a Parabricks reference dir;
  - `--ref <that dir>/genome.fa`, `--num-gpus 1`;
  - `--tmp-dir` via `mktemp -d` under `$TMPDIR` (NVMe `/scratch`), removed on exit;
  - `--memory-limit` ≈ 0.66 × task.memory (the default is half the node, and there are no cgroups);
  - emits `<meta.id>.Lb.bam` + `versions.yml`.
- `subworkflows/nf-core/fastq_align_bwa/main.nf`: the `branch` + `mix` above. The `BWA_MEM(...)` call is left byte-identical.
- `conf/modules.config`: PARABRICKS ext.args, defined unconditionally in the bwa block (in  patch they only exist under
  `aligner=='parabricks'`), and `ext.prefix = "${meta.id}.Lb"`.
- `nextflow.config`: new params `bwa_route`, `parabricks_ref` and `pb_pe_bwa_opts`, all null by default. The pilot decides
  PE `-K 320000000`, which mimics our CPU's 320 Mbp batches. With the params unset, the pipeline behaves exactly like stock.

**Commit 3 — GPU wiring and tooling:**
- `nfcore_chipseq.config`, the GPU lane block:
  ```groovy
  withName: '.*:FASTQ_ALIGN_BWA:PARABRICKS_FQ2BAM' {
    executor = 'local'; maxForks = 1        // runs inside the GPU-holding driver; one GPU
    cpus = 32; memory = 96.GB; time = 2.h   // 16 vs 32 threads: pilot decides
    containerOptions = '--nv'               // only here (the GPU is not cgroup-fenced)
    errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries = 2; maxErrors = 25          // circuit breaker: a systematic GPU failure stops the run
  }
  ```
  Also copy the SIF into `~/.singularity_cache/` under the same filename, so there is no pull and no dependency on file.
- `tools/make_parabricks_ref.sh <GRCh38|mm10>` builds gitignored `ref/parabricks/<g>/` with **real copies** of
  `genome.fa`, `.fai`, and `genome.fa.{amb,ann,bwt,pac,sa}` (≈9 GB each).
  **It never writes into `/storage/datasets/igenomes`.** That directory is group-writable, and one new file in
  `BWAIndex/version0.6.0` would invalidate all 2,013 BWA_MEM tasks.
- `tools/make_gpu_route.py` (stdlib):
  1. Pin to `cpu` every key with a CPU run COMPLETED, CACHED or FAILED, read from `nextflow log <each run of the session>`.
     The trace alone is not enough: it is overwritten on every launch.
  2. Rank the remaining keys by FASTQ Gbp, largest first, optionally SE-only or PE-only per the pilot.
  3. Fill the GPU lane until predicted lane hours balance. Cost models: GPU ≈ 50 s + 26 s/Gbp; CPU ≈ 86 s + 107 s/Gbp (SE),
     303 s + 106 s/Gbp (PE), ÷ `--cpu-slots` (default 10).
  4. Write `metadata/bwa_route_human.tsv` and print the predicted hours.
  5. `--requeue-failed` moves GPU keys that failed and have no alignment to `cpu`.
- `submit_driver.sh`, new flags:
  - `--pipeline <dir>`: drops `-r 2.1.0` for a local path.
  - `--route <tsv>`
  - `--parabricks-ref <dir>`
  - **`--resume-session <uuid>`**, which emits `-resume <uuid>`.

  A bare `-resume` takes the last history entry. That has already gone wrong twice: a `-preview` run supplied today's
  human UUID, and the Aug-2 human run adopted the mouse session.
  Human = `23fad154-a2e9-431c-b992-309fa2afc9bc`.

  The driver also gets a guard: if the route has any `gpu` key but the job holds no GPU (`SLURM_JOB_GPUS` unset), it
  aborts with the correct `sbatch` line.
- README:
  - document the GPU lane and the resume rules;
  - remove the wrong advice that `rm -rf .nextflow/` is safe.

**Standing rules for the human session once switched:**
- Always launch with the deployed hybrid directory, and never modify it.
- Never change `--fasta` / `--bwa_index` / `--gtf`.
- No `-preview`/test runs from `dbsuper_pipeline/`.
- Going back to stock `nf-core/chipseq -r 2.1.0` would re-align the GPU samples.
- At the first switch only, SAMPLESHEET_CHECK, GTF2BED and GENOME_BLACKLIST_REGIONS re-run (cheap, nothing cached depends on them).

**Mouse** stays on the stock CPU pipeline, session `1871e38d…`. It has the smallest GPU gain (all SE, small), and switching
it would re-run its finished downstream steps.

## Validation — gates the production switch and the push

Everything here runs from `/scratch/$USER/hybrid_test/`, with its own launch dir and `-w`. It can run now, while mouse
Enhancerflow finishes.

1. **Standalone GPU smoke test.** `sbatch -p gpu --gres=gpu:1`, `apptainer exec --nv`, `pbrun fq2bam` with our ref dir, on already-CPU-aligned GEO libraries:
   - SE ≤50 bp
   - SE ~64–75 bp
   - a large SE
   - 2 PE, each with and without `-K 320000000`
   - the bad-mates PE GSM5172182

   Record:
   - wall time (fixed vs per-Gbp);
   - 16 vs 32 threads;
   - GPU memory;
   - that `/storage/datasets` mtimes are unchanged;
   - whether mismatched mates are rejected;
   - whether an option such as `--fix-mate` restores the MC/MQ tags. These are the only difference  BAMs showed, and they shift Picard PCT_CHIMERAS.
2. **Record-level BAM comparison** against our existing CPU `.Lb.bam` files (read-only from `nf_work`):
   - Fields compared: flag, chrom, pos, MAPQ, CIGAR, NM, MD, AS, and mate fields.
   - Pass: SE identical except tags/header; PE ≥99.9% identical with the chosen `-K`.
   - If PE fails, route **SE only** to the GPU.
3. **Mini pipeline test.** 3 case/control pairs through the deployed hybrid pipeline with the route set to all-cpu vs all-gpu. Compare:
   - narrowPeak: identical, or Jaccard ≥ 0.99;
   - FRiP;
   - MarkDup metrics;
   - Enhancerflow ROSE2 super-enhancer calls for one pair.

   Also run a mixed-route test to confirm outputs land in `bwa/merged_library/` and stage cleanly.
4. **Production canary**, after mouse Enhancerflow finishes:
   1. Run the route helper first.
   2. Relaunch human with the hybrid directory, the right-sized config, `--resume-session 23fad154…` and an **all-cpu** route.
   3. Within ~10 min the log must show FASTQC 5,066 / TRIMGALORE 4,874 / BWA_MEM 2,013 / SORT 627 / INDEX 140 **cached**,
      with only the 3 genome-prep tasks re-running. Otherwise `scancel`.
   4. Then stop it, regenerate the route (append), and relaunch as the GPU driver.
5. **Monitoring.** Watch GPU busy time (`nvidia-smi`, trace) and the balance between the CPU and GPU lanes; rebalance only keys
   that have not started yet. After the barrier, optionally restart without `--gres`. Push `slrum_karakoram` once 1–3 pass and the canary is clean.

## Risks / not verified

- SE on Parabricks and ≤50 bp reads have never been run; tests 1–2 settle this.
- PE equivalence with our no-`-K` CPU BAMs is unproven.
- The ~50 s fixed start-up cost is fitted from PE/T2T runs.
- The GPU lane's gain is conditional on low per-task overhead (local executor) and on the CPU side being right-sized. The whole-pipeline downstream estimate is ±50%.
- There is no GPU/memory fencing (no cgroups). Parabricks pre-allocates ~38 GB of the 48 GB, so no second GPU process may run.
- Alternative not in scope: **bwa-mem2** (CPU; claims ~2× on every slot with output identical to bwa). The same route hook could host it later, after its own pilot.
