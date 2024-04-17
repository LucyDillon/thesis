#!/bin/sh
#SBATCH --time=1-23:59:00
#SBATCH --partition=bio-compute
#SBATCH --mem=50G
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --error=pa_eggnog-%A-%a.err
#SBATCH --job-name=pa_eggnog
cd /mnt/scratch2/users/40309916/Pseudomonas_genomes/genomes/prokka_analysis/prokka_output/roary_analysis/roary_core_gene

module load apps/anaconda3/5.2.0
module add diamond/0.9
source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog/

emapper.py --data_dir /mnt/scratch2/igfs-anaconda/conda-envs/eggnog/lib/python3.9/site-packages/data -i Key_sequences_updated.faa --output EggnogOutput -m diamond --cpu 80 --dbmem --dmnd_iterate no
