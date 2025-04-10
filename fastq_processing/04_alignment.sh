#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH -c 20
#SBATCH --job-name=fastq2bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/minimap2/genome.%a.txt
#SBATCH -e logs/minimap2/genome.%a.txt
#SBATCH --array=1-num_samples


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER= ## path_to_alignment_stats/genome
CONFIG=$path_to_config/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample


input_fastq=$fastq_dir/03_blaze_processed/${sample}_matched_reads.fastq.gz

LOGS_FOLDER=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/bam_stats_genome_gencode_splice
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz

SAM_FOLDER=04_minimap2_output_dir/genome/sams
BAM_FOLDER=04_minimap2_output_dir/genome/bams

mkdir -p $SAM_FOLDER
mkdir -p $BAM_FOLDER
# 
cd ~/minimap2

#remove older sam version 
rm ${SAM_FOLDER}/${sample}_stranded.sam
./minimap2 -ax splice -y --secondary=no -t ${SLURM_CPUS_PER_TASK} $REFERENCE_FASTA $input_fastq > ${SAM_FOLDER}/${sample}.sam


# # remove older bam version
rm ${BAM_FOLDER}/${sample}_stranded_sorted.bam
rm ${BAM_FOLDER}/${sample}_stranded_sorted.bam.bai

ml load samtools

samtools view -bS ${SAM_FOLDER}/${sample}.sam -o ${BAM_FOLDER}/${sample}.bam
# rm ${SAM_FOLDER}/${sample}.sam if you don't need the sam file anymore

samtools sort ${BAM_FOLDER}/${sample}.bam -o ${BAM_FOLDER}/${sample}_sorted.bam
# rm ${BAM_FOLDER}/${sample}.bam if you don't need the bam file anymore
samtools index ${BAM_FOLDER}/${sample}_sorted.bam ${BAM_FOLDER}/${sample}_sorted.bam.bai

echo "finished indexing stranded bam"
index stats ${sample}_index_stats.txt
samtools idxstats ${BAM_FOLDER}/${sample}_sorted.bam > ${LOGS_FOLDER}/${sample}_index_stats.txt


echo "finished computing stats"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
samtools flagstat ${BAM_FOLDER}/${sample}_sorted.bam >> ${LOGS_FOLDER}/${sample}_bam_flagstat.txt

## flag stats will give you a number of reads that are mapped and unmapped, secondary, and supplementary
echo "finished computing flagstats: contains percentage mapped reads"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"