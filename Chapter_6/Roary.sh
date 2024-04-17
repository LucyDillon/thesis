#!/bin/sh
#SBATCH --time=6-23:59:59
#SBATCH --partition=k2-lowpri,lowpri,bio-compute
#SBATCH --mem=350G
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=pa_roary_core-%A-%a.err
#SBATCH --job-name=core_roary
cd /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/prokka_analysis/prokka_output/roary_analysis
module load apps/anaconda3
source activate /mnt/scratch2/users/40309916/Roary

roary -e --mafft -p 32 *.gff -f roary_core_gene -g 200000
