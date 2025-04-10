#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=blaze
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/blaze.txt
#SBATCH -e logs/blaze.txt
#SBATCH --array=1-number_of_samples
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


path_to_config = ## INSERT PATH 
CONFIG=$path_to_config/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

source activate blaze_env

fastq_dir= ## INSERT PATH
input_fastq=$fastq_dir/02_nanofilt_processed/${sample}.fastq.gz

output_dir= $fastq_dir/03_blaze_processed/
mkdir -p $output_dir
output_prefix=${output_dir}/${sample}_

num_cells= ## INSERT NUM CELLS YOU EXPECT FOR THE SAMPLE. CAN BE ENCODED ON CONFIG FILE


echo "processing input fastq"
blaze --expect-cells $num_cells --output-prefix $output_prefix \
--threads $SLURM_CPUS_PER_TASK  $input_fastq  
## there is a high sensitivity option if we aren't getting enough cells 
## it is a less conservative option and will call more cells. 
## if running this mode, BLAZE suggests empty droplet removal which is 
## on by default in the later versions. 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
