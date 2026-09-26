#!/bin/bash
#SBATCH --job-name=em-seq_smk
#SBATCH -p himem
#SBATCH --mem=4G
#SBATCH --time=5-00:00:00
#SBATCH -c 1
#SBATCH -o logs/slurm/snakemake_master_%j.out
#SBATCH -e logs/slurm/snakemake_master_%j.err


work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

module load snakemake

snakemake --snakefile EM-seq_bwa-mem3_workflow.smk \
    --configfile config_ref_boundle.yaml \
    --cluster-config cluster_config.yaml \
    --jobs 120 \
    --cluster "sbatch -p {cluster.partition} -c {threads} --mem={cluster.mem} \
               -t {cluster.time} {cluster.extra} \
               -o logs/slurm/%x_%j.out -e logs/slurm/%x_%j.err" \
    --latency-wait 60 \
    --rerun-incomplete
