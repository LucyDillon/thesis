#!/bin/sh
#SBATCH --time=7-23:59:00
#SBATCH --partition=bio-compute
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --job-name=eggnog
#SBATCH --error=diamond_dn-%A-%a.err
echo $SLURM_JOB_ID

module load apps/anaconda3
module add  diamond/0.9
cd /mnt/scratch2/users/********/PATRIC_Genomes/gff_files

diamond blastp -d /mnt/scratch2/igfs-anaconda/conda-envs/eggnog/lib/python3.9/site-packages/data/eggnog_proteins.dmnd -q /mnt/scratch2/users/40309916/PATRIC_Genomes/gff_files/split_prodigal_proteinsdn --more-sensitive --threads 2 -e 0.001000 -o /mnt/scratch2/users/40309916/PATRIC_Genomes/gff_files/split_prodigal_proteinsdn_sequences.txt --top 3 --query-cover 0 --subject-cover 0
# Version: diamond v0.9.24.125
