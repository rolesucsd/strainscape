#!/bin/bash -l
#SBATCH --partition=short
#SBATCH --job-name=process-fastq
#SBATCH --time=2:00:00
#SBATCH --mem=30G

coverage_dir="/ddn_scratch/roles/strain_analysis/iHMP/coverages/BU/coverage_results"
output_coverage_summary="/ddn_scratch/roles/strain_analysis/iHMP/coverages/BU/output/coverage_summary.csv"

# Create or overwrite the summary files with headers
echo "Sample,Average_Coverage,Coverage_Percent" > "$output_coverage_summary"

# Loop through all coverage files
for coverage_file in "$coverage_dir"/*_coverage.txt; do
    # Extract the sample name
    sample=$(basename "$coverage_file" _coverage.txt)

    # Calculate average coverage and coverage percent using awk
    # coverage file format: contig position coverage
    read avg_cov cov_pct < <(awk '
    {
        sum+=$3; count++;
        if($3>0) covered++
    }
    END {
        if(count>0) {
            avg = sum/count
            pct = (covered/count)*100
        } else {
            avg = 0
            pct = 0
        }
        print avg, pct
    }' "$coverage_file")

    # Append to coverage_summary.csv
    echo "${sample},${avg_cov},${cov_pct}" >> "$output_coverage_summary"
done

echo "Processed coverage files. Results are in $output_coverage_summary
