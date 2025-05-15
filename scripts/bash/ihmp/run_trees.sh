#!/bin/bash
#SBATCH --job-name=chromosome_trees
#SBATCH --output=chromosome_trees.log
#SBATCH --error=chromosome_trees.err
#SBATCH --time=60:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8

source /home/roles/anaconda3/bin/activate
python create_trees.py --genes_file ../reference/bakta/wolr2_reference.tsv --output_dir ../reference/trees --threads 8
