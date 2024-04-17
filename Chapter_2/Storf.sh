#!/bin/sh
#SBATCH --time=10:59:59
#SBATCH --partition=k2-medpri,bio-compute
#SBATCH --mem=20G
#SBATCH --mail-user=ldillon05@qub.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL

module load apps/pip_python310/22.3.1/python3-3.10.5

StORF-Reporter -anno Prokka Multiple_Out_Dirs -p ./ -overwrite True -sout True
