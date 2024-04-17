#!/bin/bash
#SBATCH --time=13-23:59:00
#SBATCH --partition=lowpri
#SBATCH --job-name=weka
#SBATCH --mem=50G
cd /mnt/scratch2/users/********/eggnog_arff

module load weka/3.8.6/java-8.0.151

for i in *.arff; do java weka.classifiers.trees.J48 -t $i -C 0.25 -M 2 -g > $i.dot; done
