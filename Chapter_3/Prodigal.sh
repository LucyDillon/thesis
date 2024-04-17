#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --partition=lowpri
#SBATCH --job-name=prodigal_run
cd /mnt/scratch2/users/********/PATRIC_Genomes
module load apps/anaconda3
source activate /mnt/scratch/igfs-anaconda/conda-envs/prokka
for i in *.fna; do prodigal -i $i -o $i.genes -f gff -a $i.proteins.gff; done
