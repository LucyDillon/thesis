#!/bin/sh
#SBATCH --time=06-23:59:59
#SBATCH --partition=k2-lowpri
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --job-name=sourmash
cd /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/Genomes
module load apps/anaconda3
source activate /mnt/scratch2/users/40309916/panaroo

sourmash sketch dna *.fna
