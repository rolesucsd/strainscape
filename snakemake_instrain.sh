#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1       # match your heaviest rule (megahit)
#SBATCH --mem=5G                # assemblies often need tens of GB
#SBATCH --time=48:00:00          # give plenty of time for co‐assembly + binning
#SBATCH --job-name=snakemake_iHMP

source /home/roles/anaconda3/bin/activate
#mamba env remove -n bakta

#mamba create -c bioconda -n bakta bakta=1.10.3

#conda activate bakta

#rm -r .snakemake

#rm -r /ddn_scratch/roles/Panpiper/panpiper/databases/bakta/db

#bakta_db download --output /ddn_scratch/roles/strain_analysis/iHMP/bakta

#conda deactivate

conda activate panpiper

export TMPDIR=/ddn_scratch/roles/strain_analysis/tmp

snakemake -s instrain_mags.smk --unlock

#snakemake -s instrain.smk --touch -c 1

snakemake \
  -s instrain_mags.smk \
  -j 30 \
  --use-conda \
  --keep-going \
  --stats stats.json \
  --profile /ddn_scratch/roles/Panpiper/profile \
  --touch
