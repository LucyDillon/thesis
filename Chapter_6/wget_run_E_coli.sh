#!/bin/bash

#SBATCH --time=06-23:59:59
#SBATCH --partition=lowpri
#SBATCH --job-name=E_coli_WGET_run
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
cd /mnt/scratch2/users/40309916/E_coli_genomes
for i in `cat E_coli_genome_ids_final.txt`; do wget -qN "ftp://ftp.patricbrc.org/genomes/$i/$i.*fna";
done
