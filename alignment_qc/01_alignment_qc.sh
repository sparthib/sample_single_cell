#!/bin/bash

#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 10
#SBATCH --job-name=highqualbam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/minimap2/genome/primary_over_30.%a.txt
#SBATCH -e logs/minimap2/genome/primary_over_30.%a.txt
#SBATCH --array=1-8
#SBATCH -t 4-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER= ## path_to_alignment_stats/genome/primary_over_30
CONFIG=$path_to_config/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)


echo $sample
alignment_dir= ## path_to_alignment_outputs
input_dir=$alignment_dir/genome/bams
output_dir=$input_dir/primary_over_30_chr_only
mkdir -p $output_dir
bam=$input_dir/${sample}_sorted.bam

ml load samtools 

samtools view -h -b -q 30 -F 0x800 $bam > $output_dir/${sample}_primary_over_30.bam
echo "finished filtering bam by mapping quality and removing duplicates"

samtools sort $output_dir/${sample}_primary_over_30.bam -o $output_dir/${sample}_primary_over_30_sorted.bam
samtools index $output_dir/${sample}_primary_over_30_sorted.bam $output_dir/${sample}_primary_over_30_sorted.bam.bai
echo "finished sorting and indexing bam"

samtools view -h -b $output_dir/${sample}_primary_over_30_sorted.bam chr1 chr2 chr3 chr4 \
chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 \
chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > $output_dir/${sample}_primary_over_30_chr_only.bam
echo "finished subsetting bam by chromosome"

samtools sort $output_dir/${sample}_primary_over_30_chr_only.bam -o $output_dir/${sample}_primary_over_30_chr_only_sorted.bam
samtools index $output_dir/${sample}_primary_over_30_chr_only_sorted.bam $output_dir/${sample}_primary_over_30_chr_only_sorted.bam.bai
echo "finished sorting and indexing bam"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
samtools flagstat $output_dir/${sample}_primary_over_30_chr_only_sorted.bam >> ${LOGS_FOLDER}/${sample}_bam_flagstat.txt

echo "finished computing flag stats"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
