#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medpri
#SBATCH --job-name=RGI_WGS_sequences
#SBATCH --array=1-167%100
echo $SLURM_ARRAY_TASK_ID
echo $SLURM_JOB_ID
module load apps/anaconda3
source activate /mnt/scratch/igfs-anaconda/conda-envs/rgi
cd /mnt/scratch2/users/********
#file_list=$(ls *.fna | sed -n ${SLURM_ARRAY_TASK_ID}p)
for file_list in `ls *.fna | sed -n $(expr $(expr ${SLURM_ARRAY_TASK_ID} \* 100) - 99),$(expr ${SLURM_ARRAY_TASK_ID} \* 100)p`;
do rgi main --input_sequence /mnt/scratch2/users/40309916/$file_list \--output_file /mnt/scratch2/users/40309916/RGI_output/$file_list --input_type contig --local --clean;
done
