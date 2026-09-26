#!/bin/bash
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=46000M
#SBATCH -t 3-10:00 # Runtime in D-HH:MM
#SBATCH -J WGBS_MethylDackel
#SBATCH --array=0-12 # job array index - 13 samples, index 0-12

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

ref_genome="/seqs_for_alignment_pipelines/genome.fa"

module load samtools/1.20
module load MethylDackel/0.6.1

if [ ! -d "${work_path}/MethylDackel" ]; then
        mkdir -p "${work_path}/MethylDackel"
fi

##assign each sample into different ARRAYID
completed_txt="${work_path}/completed_samples.txt"
mapfile -t complete_samples < <(awk -F'\t' 'NR>1 && $4=="COMPLETE" {print $1}' "$completed_txt")

sample="${complete_samples[${SLURM_ARRAY_TASK_ID}]}"
echo "sample=${sample}"

input_bam="${work_path}/merged_bam/${sample}.markdup.bam"
echo "input_bam=${input_bam}"

##run MethylDackel mbias, capture stdout for maxval extraction
##stdout contains the per-position table (Strand Read Position nMethylated nUnmethylated)
##which is used to determine the actual max read length
mbias_output=$(MethylDackel mbias \
    --txt \
    -@ 4 \
    ${ref_genome} \
    ${input_bam} \
    ${work_path}/MethylDackel/${sample} 2>&1)

echo "${mbias_output}"
echo "** ${sample} mbias done **"

##get actual max read length from captured mbias stdout (column 3 = Position)
maxval=$(echo "${mbias_output}" | awk 'NR>1 && $3~/^[0-9]+$/ { if ($3+0 > max) max=$3+0 } END { print max }')
echo "maxval=${maxval}"

##parse --OT and --OB directly from SVG files
OT_svg="${work_path}/MethylDackel/${sample}_OT.svg"
OB_svg="${work_path}/MethylDackel/${sample}_OB.svg"

otsvg=$(grep -oE '>--OT [0-9]+,[0-9]+,[0-9]+,[0-9]+' "${OT_svg}" | awk '{print $2}')
obsvg=$(grep -oE '>--OB [0-9]+,[0-9]+,[0-9]+,[0-9]+' "${OB_svg}" | awk '{print $2}')

echo "otsvg=${otsvg}"
echo "obsvg=${obsvg}"

##validate parsed values
if [[ -z "${otsvg}" || -z "${obsvg}" || -z "${maxval}" ]]; then
    echo "ERROR: Failed to parse OT/OB from SVG or maxval from mbias output. Exiting."
    exit 1
fi

##split into arrays
IFS=',' read -ra otv <<< "${otsvg}"
IFS=',' read -ra obv <<< "${obsvg}"

##fix 0 values: MethylDackel 0.6.1 bug - 0 causes "Invalid bounds string"
##positions A,C (index 0,2) = start positions: replace 0 with 1
##positions B,D (index 1,3) = end positions:   replace 0 with maxval
for i in 0 2; do
    if [[ "${otv[i]}" -eq 0 ]]; then otv[i]=1; fi
    if [[ "${obv[i]}" -eq 0 ]]; then obv[i]=1; fi
done
for i in 1 3; do
    if [[ "${otv[i]}" -eq 0 ]]; then otv[i]=${maxval}; fi
    if [[ "${obv[i]}" -eq 0 ]]; then obv[i]=${maxval}; fi
done

OT="--OT ${otv[0]},${otv[1]},${otv[2]},${otv[3]}"
OB="--OB ${obv[0]},${obv[1]},${obv[2]},${obv[3]}"

echo "OT=${OT}"
echo "OB=${OB}"

##append trimming parameters to summary file for QC review
summary_file="${work_path}/MethylDackel/mbias_trimming_summary.tsv"
if [ ! -f "${summary_file}" ]; then
    echo -e "sample\tOT_original\tOB_original\tmaxval\tOT_fixed\tOB_fixed" > "${summary_file}"
fi
echo -e "${sample}\t${otsvg}\t${obsvg}\t${maxval}\t${OT#--OT }\t${OB#--OB }" >> "${summary_file}"
echo "** ${sample} trimming parameters saved to ${summary_file} **"

##run MethylDackel extract with bounds, variant site exclusions and min depth filter
#MethylDackel extract \
#    --CHH \
#    --CHG \
#    $OT \
#    $OB \
#    --minOppositeDepth 5 \
#    --maxVariantFrac 0.1 \
#    --minDepth 5 \
#    --methylKit \
#    -@ 4 \
#    -o ${work_path}/MethylDackel/${sample} \
#    ${ref_genome} \
#    ${input_bam} && echo "** ${sample} MethylDackel extract done **"

##run MethylDackel extract with bounds, variant site exclusions and min depth filter
##capture stderr to extract number of positions excluded as variants
extract_output=$(MethylDackel extract \
    --CHH \
    --CHG \
    $OT \
    $OB \
    --minOppositeDepth 5 \
    --maxVariantFrac 0.1 \
    --minDepth 5 \
    --methylKit \
    -@ 4 \
    -o ${work_path}/MethylDackel/${sample} \
    ${ref_genome} \
    ${input_bam} 2>&1)

echo "${extract_output}"
echo "** ${sample} MethylDackel extract done **"

##parse number of positions excluded due to variants
positions_excluded=$(echo "${extract_output}" | grep -oP '[0-9]+ positions were excluded' | grep -oP '^[0-9]+')
[[ -z "${positions_excluded}" ]] && positions_excluded="0"
echo "positions_excluded=${positions_excluded}"

##append the number of excluded likely being variants positions for QC review
excluded_position_file="${work_path}/MethylDackel/excluded_positions_summary.tsv"
if [ ! -f "${excluded_position_file}" ]; then
   echo -e "sample\tpositions_excluded_likely_variants" > "${excluded_position_file}"
fi

echo -e "${sample}\t${positions_excluded}" >> "${excluded_position_file}"
echo "** ${sample} QC metrics saved to ${excluded_position_file} **"
