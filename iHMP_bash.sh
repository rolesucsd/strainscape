#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=500GB

source /home/roles/anaconda3/bin/activate

python iHMP_main.py 