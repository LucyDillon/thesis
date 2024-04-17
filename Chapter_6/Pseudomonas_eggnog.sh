#!/bin/sh
#SBATCH --time=13-23:59:00
#SBATCH --partition=k2-lowpri,bio-compute
#SBATCH --mem=60G
#SBATCH --cpus-per-task=24
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=pa_eggnog-%A-%a.err
#SBATCH --array=1-71%100
#SBATCH --job-name=pa_eggnog
cd /mnt/scratch2/users/40309916/Pseudomonas_genomes/genomes/prokka_analysis/prokka_output/protein_fasta_file

module load apps/anaconda3/5.2.0
module add diamond/0.9
source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog/

file_list=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)
for file_list in `ls *.faa | sed -n $(expr $(expr ${SLURM_ARRAY_TASK_ID} \* 100) - 99),$(expr ${SLURM_ARRAY_TASK_ID} \* 100)p`; 
  do emapper.py --data_dir /mnt/scratch2/igfs-anaconda/conda-envs/eggnog/lib/python3.9/site-packages/data -i $file_list --output $file_list.output -m diamond  --dbmem --dmnd_iterate no --cpu 24; 
done
