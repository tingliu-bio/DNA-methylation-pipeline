# EM-seq Pipeline (bwa-mem3 / hg38)

This repo documents a Snakemake-based pipeline for processing Enzymatic Methyl-seq (EM-seq) whole-genome bisulfite sequencing data, using `bwa-mem3` with methylation-aware alignment. The pipeline handles multi-lane paired-end samples on a SLURM HPC cluster.

---

## Overview

The full pipeline flow:

<pre>
FASTQ (per lane)
   │
   ▼
1. Quality Trimming        TrimGalore --paired
   │
   ▼
2. Alignment (per lane)    bwa-mem3 --meth → samtools sort → .bam
   │
   ▼
3. Merge Lanes             samtools merge → per-sample .bam
   │
   ▼
4. Filter Non-converted    mark-nonconverted-reads.py → samtools view -F 512
   Reads
   │
   ▼
5. Mark Duplicates         sambamba markdup
   │
   ▼
6. BAM QC                  samtools flagstat + stats
   │
   ▼
7. Methylation Calling     MethylDackel mbias → extract → _CpG.methylKit
</pre>

---

## 0. Build Reference Index

Before the first run, build the methylation-aware bwa-mem3 index. This step is only needed once per reference genome and requires a high-memory node (~150 GB) with AVX2 support.

The script runs:
```bash
bwa-mem3 index --meth --max-memory 150G genome.fa
```

> **Note:** The `--constraint=avx2` SLURM flag is required — bwa-mem3 will fail on nodes without AVX2 support.

---

## 1. Input: Sample Table

The pipeline takes a CSV with one row per sequencing lane:

| Column     | Description                          |
|------------|--------------------------------------|
| `sample`   | Sample ID                            |
| `lane`     | Lane identifier (e.g. `L001`)        |
| `flowcell` | Flowcell ID                          |
| `R1`       | Path to R1 FASTQ (.fq.gz)           |
| `R2`       | Path to R2 FASTQ (.fq.gz)           |

Multiple lanes per sample are merged automatically at step 3.

---

## 2. Configuration

Edit `config.yaml` before running:

```yaml
samples_csv:      "/path/to/master_samples_input.csv"
work_dir:         "/path/to/output_directory"              # all output written here
ref_genome:       "/path/to/bwa-mem3_hg38/genome.fa"       # bwa-mem3 index
ref_bundle:       "/path/to/bwa-mem3_hg38/"                # contains mark-nonconverted-reads.py
ref_methyldackel: "/path/to/hg38/genome.fa"                # FASTA for MethylDackel
```

---

## 3. Running the Pipeline

```bash
sbatch submit_smk_sbatch.sh
```

This submits the Snakemake master job, which then dispatches up to 120 concurrent SLURM jobs via `cluster.yaml`.

---

## 4. Key Design Decisions

**Non-converted read filtering (step 4)**
`mark-nonconverted-reads.py` flags reads that failed bisulfite conversion; `samtools view -F 512` then removes them. This step is placed after lane merging to maximize the read depth available for filtering decisions.

**Intermediate files as `temp()`**
Per-lane BAMs, merged BAMs, and filtered BAMs are all marked as `temp()` in Snakemake and automatically deleted after downstream rules complete, to keep disk usage manageable.

**MethylDackel: mbias before extraction**
MethylDackel runs an M-bias analysis first to identify per-read-position methylation artifacts, then uses the trimming parameters from that run for the final extraction step. A per-run `mbias_trimming_summary.tsv` is generated for QC review.

**Read group tags**
Each per-lane BAM carries a full `@RG` tag (`ID`, `SM`, `CN`, `PL`, `LB`, `PU`) set at alignment time, so sample identity is preserved through merging.

---

## 5. Output

All output is written under `work_dir/`:

```
{work_dir}/
├── BWA/
│   └── {sample}.markdup.bam       # Final deduplicated BAM (retained)
├── MethylDackel/
│   ├── {sample}_CpG.methylKit     # CpG methylation (methylKit format)
│   ├── {sample}_CHG.methylKit
│   ├── {sample}_CHH.methylKit
│   └── mbias_trimming_summary.tsv
├── QC/
│   ├── {sample}.flagstat
│   └── {sample}.stats
└── logs/                          # Per-rule logs (trim, bwa, merge, etc.)
```

> Intermediate BAMs (per-lane, merged, filtered) are automatically removed by Snakemake after the pipeline completes.

---

## Tools and Versions

| Tool                        | Version  | Purpose                              |
|-----------------------------|----------|--------------------------------------|
| Snakemake                   | ≥ 7.0    | Workflow management                  |
| TrimGalore                  | 0.6.6    | Adapter trimming                     |
| bwa-mem3                    | —        | Methylation-aware alignment          |
| samtools                    | 1.20     | BAM processing and QC                |
| sambamba                    | —        | Duplicate marking                    |
| mark-nonconverted-reads.py  | —        | Non-converted read filtering         |
| MethylDackel                | 0.6.1    | Methylation extraction               |

---

## Notes

- Intermediate BAMs are marked as `temp()` in Snakemake and automatically removed after downstream rules complete, to save disk space.
- Non-converted read filtering uses `samtools view -F 512`, which removes reads flagged as failing QC by `mark-nonconverted-reads.py`.
- MethylDackel runs mbias analysis first to determine optimal trimming positions before extraction.
- Pipeline is under active development; WGBS workflow coming soon.
