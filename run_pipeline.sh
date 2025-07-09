#!/bin/bash

# iHMP Analysis Pipeline
# This script runs the complete iHMP analysis pipeline from raw data to final results

set -e  # Exit on error

# Configuration
CONFIG_FILE="config.yaml"
SNAKEMAKE_ENV="snakemake/envs/python.yaml"
R_ENV="snakemake/envs/r.yaml"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status messages
print_status() {
    echo -e "\033[1;34m==>\033[0m $1"
}

# Function to print error messages
print_error() {
    echo -e "\033[1;31mError:\033[0m $1"
    exit 1
}

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    print_error "Config file not found: $CONFIG_FILE"
fi

# Run Snakemake workflow
print_status "Running Snakemake workflow..."
snakemake --use-conda --conda-frontend conda \
    --configfile "$CONFIG_FILE" \
    --conda-prefix "$SNAKEMAKE_ENV" \
    --cores 4 || print_error "Snakemake workflow failed"

# Run Python analysis
print_status "Running Python analysis..."
python src/python/ihmp/main.py \
    --config "$CONFIG_FILE" || print_error "Python analysis failed"