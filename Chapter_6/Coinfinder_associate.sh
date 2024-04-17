#!/bin/sh
#SBATCH --time=15-23:59:59
#SBATCH --partition=bio-himem
#SBATCH --mem=350G
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=pa_coinfinder-%A-%a.err
#SBATCH --job-name=PA_coinfinder
cd /mnt/scratch2/users/40309916/Pseudomonas_genomes/genomes/prokka_analysis/prokka_output/roary_analysis/roary_core_gene
module load apps/anaconda3
source activate /mnt/scratch2/users/40309916/coinfinder

coinfinder -i gene_presence_absence1.csv -I -p tree_1e-06.nwk -o Pseudomonas_coinfinder --associate -x 32
