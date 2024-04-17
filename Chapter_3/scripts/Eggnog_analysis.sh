#!/bin/sh
#SBATCH --time=6-23:59:00
#SBATCH --partition=lowpri
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=eggnog
echo $SLURM_JOB_ID

cd /mnt/scratch2/users/********/PATRIC_Genomes/gff_files
module load apps/anaconda3
source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog

emapper.py -m no_search --annotate_hits_table diamond_modified_files_3.txt --override -o diamond_modified_3_eggnog.txt --cpu 16 --dbmem
#Eggnog version emapper-2.1.6
