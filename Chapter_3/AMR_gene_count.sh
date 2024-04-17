#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medpri
#SBATCH --job-name=gene_count_AMRFinder
#SBATCH --array=1-167%100
echo $SLURM_ARRAY_TASK_ID
echo $SLURM_JOB_ID
cd /mnt/scratch2/users/********
for i in *.txt.txt.txt;
do cat AMRFinder_genes.txt | xargs -I % grep -c "%$" $i > $i.AMRFindercount.txt;
done
