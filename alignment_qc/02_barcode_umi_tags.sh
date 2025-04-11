#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=add_barcode
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/add_barcode.%a.txt
#SBATCH -e logs/add_barcode.%a.txt
#SBATCH --array=1-8
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

CONFIG= #INSERT PATH
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
num_cells=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $5}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

## input_bam_dir= #INSERT PATH
## output_bam_dir=  #INSERT PATH
input_bam=$input_bam_dir/${sample}_primary_over_30_chr_only.bam
output_bam=$output_bam_dir/${sample}_tagged.bam
rm $output_bam
ml load python/3.10.13
python3 02_barcode_umi_tags.py $input_bam $output_bam


#sort tagged file 
input_bam=$output_bam_dir/${sample}_tagged.bam
output_bam=$output_bam_dir/sorted/${sample}_tagged.bam
output_bai=$output_bam_dir/sorted/${sample}_tagged.bam.bai
rm $output_bam
rm $output_bai

ml load samtools

samtools sort $input_bam -o $output_bam
# rm ${BAM_FOLDER}/${sample}.bam
samtools index $output_bam $output_bai



echo "**** Job ends ****"
date +"%Y-%m-%d %T"