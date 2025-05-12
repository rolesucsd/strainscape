#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roles@ucsd.edu

export TMPDIR=/panfs/roles/tmp

source /home/roles/anaconda3/bin/activate

conda activate sepsis

dRep dereplicate ehormaechi/drep2 -sa 0.98 --ignoreGenomeQuality