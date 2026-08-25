#!/bin/bash
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=46000M
#SBATCH -t 3-10:00 # Runtime in D-HH:MM
#SBATCH -J WGBS_trim
#SBATCH --array=0-64 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

ref_genome="/seqs_for_alignment_pipelines/genome.fa"

module load bwa/0.7.15;
module load samtools/1.20;
module load bwa_meth/0.2.7;
module load python3/3.10.9
module load trim_galore/0.6.6

file_path="${work_path}/raw_fq"

if [ ! -d "${work_path}/fastq_trimmed" ]; then
        mkdir -p "${work_path}/fastq_trimmed"
fi

if [ ! -d "${work_path}/BWA" ]; then
        mkdir -p "${work_path}/BWA"
fi

##assign each sample into different ARRAYID
samples=${work_path}/fq1_file_list.txt

prefixs=($(cat $samples))
echo  "prefixs="$prefixs
input_fq1=${prefixs[${SLURM_ARRAY_TASK_ID}]}

##get the prefix name
if [[ "$input_fq1" =~ ${file_path}/(.*)-R1.fastq.gz ]]
 then
	echo "input_fq1=$input_fq1"
	prefix="${BASH_REMATCH[1]}"
	echo "prefix=$prefix"

trim_galore \
  --paired \
  --gzip \
  --cores 4 \
  --output_dir "${work_path}/fastq_trimmed" \
  "${input_fq1}" "${file_path}/${prefix}-R2.fastq.gz" && echo "** $prefix trim_galore done  **"

bwameth.py -t4 --read-group  "@RG\tID:"${prefix}"\tSM:${prefix}\tCN:PMGC\tPL:Illumina\tLB:${prefix}\tPU:lane_name" \
  --reference ${ref_genome} \
  ${work_path}/fastq_trimmed/${prefix}-R1_val_1.fq.gz \
  ${work_path}/fastq_trimmed/${prefix}-R2_val_2.fq.gz > ${work_path}/BWA/${prefix}.sam && echo "** $prefix bwameth done  **"

samtools sort -@4 -m 2G -O bam -o ${work_path}/BWA/${prefix}.bam \
  -T ${prefix} \
  ${work_path}/BWA/${work_path}/BWA/${prefix}.sam && echo "** $prefix sam2bam done  **"

samtools index ${work_path}/BWA/${prefix}.bam && echo "** $prefix bam index done  **"

/mark-nonconverted-reads.py \
  --reference ${ref_genome} \
  --bam ${work_path}/BWA/${prefix}.bam | samtools view -b > ${work_path}/BWA/${prefix}_filtered.bam && echo "** $prefix mark-nonconverted done  **"
samtools index ${work_path}/BWA/${prefix}_filtered.bam && echo "** $prefix filtered bam index done  **"



fi
