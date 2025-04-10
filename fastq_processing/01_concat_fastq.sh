#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=concat_fqs
#SBATCH -c 2
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/concat/concat_fqs.%a.txt
#SBATCH -e logs/concat/concat_fqs.%a.txt
#SBATCH --array=1-number_of_samples

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


path_to_config = ## INSERT PATH 
CONFIG=$path_to_config/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample
output_fastq_dir = ##INSERT PATH
OUTPUT_FOLDER=$output_fastq_dir/01_input_fastqs

cd $INPUT_FOLDER/fastq_pass/
  cat *.fastq.gz > $OUTPUT_FOLDER/${sample}.fastq.gz
  
  
  
  

echo "**** Job ends ****"
date +"%Y-%m-%d %T"