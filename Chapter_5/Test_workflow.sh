#!/bin/sh
#SBATCH --time=1-23:59:00
#SBATCH --partition=k2-lowpri,bio-compute
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=SM_test-%A-%a.err
#SBATCH --job-name=snakemake_test

module load apps/anaconda3
module load snakemake/V5.31.1_Python3.8.5
module add  diamond/0.9
source activate /mnt/scratch2/igfs-anaconda/conda-envs/prokka
source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog

snakemake --cores 8 --snakefile Snakefile.smk Eggnog/Eggnog_output.txt
