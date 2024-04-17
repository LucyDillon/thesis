#!/bin/sh
#SBATCH --time=23:59:59
#SBATCH --partition=bio-compute
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name=EC_plot_csv
cd /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/Genomes
module load apps/anaconda3
source activate /mnt/scratch2/users/40309916/panaroo

sourmash plot distances.cmp --csv plot_output.csv
