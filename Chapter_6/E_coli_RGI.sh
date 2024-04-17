#!/bin/bash
#SBATCH --time=11-23:59:59
#SBATCH --partition=bio-compute
#SBATCH --job-name=RGI_E_coli
#SBATCH	--mem=50G
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=RGI_E_coli-%A-%a.err
#SBATCH --array=1-96%100
echo $SLURM_JOB_ID
module load apps/anaconda3
source activate /mnt/scratch/igfs-anaconda/conda-envs/rgi
cd /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/prokka_analysis/prokka_output/protein_fasta_files

file_list=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)
for file_list in `ls *.faa | sed -n $(expr $(expr ${SLURM_ARRAY_TASK_ID} \* 100) - 99),$(expr ${SLURM_ARRAY_TASK_ID} \* 100)p`; do
rgi main --input_sequence /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/prokka_analysis/prokka_output/protein_fasta_files/$file_list \--output_file /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/prokka_analysis/prokka_output/roary_analysis/RGI_output/$file_list --input_type protein --local --clean;
done
