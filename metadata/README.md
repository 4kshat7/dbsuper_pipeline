# ENCODE H3K27ac Organism | *Homo sapiens* Cohort — Summary

**Source query:** `encodefetch --assay-title "Histone ChIP-seq" --target-label H3K27ac --organism "Homo sapiens" --file-type fastq --status released`

**Inputs:**
- [raw/encode_results/manifest.tsv](raw/encode_results/manifest.tsv) — one row per fastq file (4,277 rows)
- [raw/encode_results/nfcore_chipseq_samplesheet.csv](raw/encode_results/nfcore_chipseq_samplesheet.csv) — nf-core/chipseq samplesheet (4,277 rows)

---

## 1. File / Experiment Counts

| Metric | Value |
|---|---|
| Total fastq files | **4,277** |
| Total unique experiments | **913** |
| ChIP (H3K27ac) experiments | **458** |
| Control (input) experiments | **455** |
| ChIP fastq files | 2,013 |
| Control fastq files | 2,264 |
| Total data volume | **3.41 TB** |
| └─ ChIP | 1.47 TB |
| └─ Control | 1.95 TB |
| Organism | *Homo sapiens* (100%) |
| File format | fastq (100%) |
| File status | released (100%) |

---

## 2. Single-end vs Paired-end

| Run type | Files | Experiments |
|---|---|---|
| **Single-ended** | 4,094 (95.7%) | 834 |
| **Paired-ended** | 183 (4.3%) | 79 |

Breakdown by role:

| | Single | Paired |
|---|---|---|
| ChIP files | 1,916 | 97 |
| Control files | 2,178 | 86 |

> Paired-ended fastq files come as R1+R2 pairs, so 183 files ≈ 91 paired pairs across 79 experiments.

---

## 3. Replicate Structure (ChIP only, N=458 experiments)

| Bio replicates | # experiments |
|---|---|
| 1 | 301 (65.7%) |
| 2 | 131 (28.6%) |
| 3 | 26 (5.7%) |

Mean biological replicates per ChIP experiment: **1.40**

---

## 4. Biology

- **Unique biosamples (cell types / tissues):** 173
- **Top biosamples:** dorsolateral prefrontal cortex (1,164 files), heart right ventricle (197), heart left ventricle (186), spleen (172), skin epidermis (142), A549 (119), transverse colon (108), pancreas (76), adrenal gland (62)

**ChIP experiments by biosample classification:**

| Classification | Count |
|---|---|
| tissue | 296 |
| cell line | 86 |
| primary cell | 49 |
| in vitro differentiated cells | 27 |

**Donor life stage (ChIP):** adult 334 · embryonic 48 · child 21 · newborn 12 · unknown 23 · NA 20
**Donor sex:** female 2,264 files · male 1,791 files · unknown 222 files

---

## 5. Sequencing Platforms

| Platform | Files |
|---|---|
| Illumina NextSeq 500 | 2,091 |
| Illumina HiSeq 2500 | 1,094 |
| Illumina HiSeq 2000 | 853 |
| Illumina HiSeq 4000 | 73 |
| Illumina Genome Analyzer IIx | 69 |
| Illumina Genome Analyzer (others) | 79 |
| Illumina NovaSeq 6000 | 18 |

---

## 6. Metadata Tables Written

- [experiment_metadata.tsv](experiment_metadata.tsv) — 913 rows × 12 cols (one row per experiment — compact)
- [experiment_metadata_full.tsv](experiment_metadata_full.tsv) — 913 rows × 30 cols (full biology + donor + provenance)

Columns in the full table: experiment_accession, is_control, matched_control, target_label, assay_title, description, lab, biosample_term_id, biosample_term_name, classification, organ_slims, cell_slims, developmental_slims, system_slims, organism, donor_accession, donor_sex, donor_life_stage, donor_age, donor_age_units, donor_ethnicity, run_type, platform, bio_replicate_count, tech_replicate_count, replication_type, date_released, status_exp, n_fastq_files, total_bytes.

---

## 7. Caveat on the nf-core samplesheet

The `nfcore_chipseq_samplesheet.csv` has **one row per fastq file (4,277 rows)** and its `sample` column covers **both ChIP and control experiments as separate samples (913 uniques)**. In the standard nf-core/chipseq template, only ChIP samples are rows and controls are referenced via the `control` column — double-check whether your pipeline run expects the standard layout or this expanded one before Stage 1.
