# iHMP Analysis Pipeline

A comprehensive pipeline for analyzing strain-level evolution in the iHMP (Integrative Human Microbiome Project) dataset, with a focus on reproducibility, statistical rigor, and clear documentation.

## Pipeline Overview

The pipeline consists of several interconnected modules:

1. **Assembly** (`assembly.smk`):
   - Co-assembly of metagenomic reads
   - Contig filtering and quality control
   - Generation of STB files for strain tracking

2. **Mapping** (`mapping.smk`):
   - Read mapping to assembled contigs
   - BAM file processing and sorting
   - Mapping statistics generation

3. **MAGs** (`mags.smk`):
   - Metagenome-assembled genome binning
   - Quality assessment of bins
   - Taxonomic classification

4. **InStrain** (`instrain.smk`):
   - Strain-level profiling
   - SNV calling and filtering
   - Scaffold information processing
   - Combined analysis across samples

5. **SNP Tracking** (`snp_tracking.smk`):
   - Mutation trajectory analysis
   - Frequency change tracking
   - Cluster analysis of mutations

## Statistical Methodology and Parameter Justification

### Assembly Parameters
- **Minimum Contig Length (1000bp)**: Balances assembly quality with strain-level resolution
- **Minimum Coverage (10x)**: Ensures reliable assembly of true genomic regions
- **Maximum Contamination (10%)**: Based on CheckM2 recommendations for high-quality MAGs

### Mapping Parameters
- **Minimum Mapping Quality (1)**: Ensures reliable read placement
- **Minimum Coverage (10x)**: Provides sufficient depth for variant calling
- **Maximum Insert Size**: Based on library preparation protocol

### InStrain Parameters
- **Minimum Coverage (5x)**: Ensures reliable variant detection
- **Minimum Frequency (0.05)**: Filters out sequencing errors while retaining true variants
- **P-value Threshold (0.05)**: Standard statistical significance level

### SNP Tracking Parameters
- **Minimum Slope (0.01/week)**: Captures significant frequency changes over time
- **Minimum Absolute Change (0.1)**: Ensures biological relevance of detected changes
- **Cluster Distance Threshold**: Based on hierarchical clustering optimization

## Directory Structure

```
.
├── snakemake/
│   ├── Snakefile           # Main pipeline file
│   ├── rules/             # Individual rule files
│   │   ├── assembly.smk   # Assembly rules
│   │   ├── mapping.smk    # Mapping rules
│   │   ├── mags.smk       # MAG binning rules
│   │   ├── instrain.smk   # InStrain analysis rules
│   │   └── snp_tracking.smk # SNP tracking rules
│   ├── config/            # Configuration files
│   │   ├── config.yaml    # Main configuration
│   │   ├── patients.yaml  # Patient information
│   │   └── samples.yaml   # Sample information
│   ├── envs/             # Conda environment files
│   └── scripts/          # Pipeline scripts
├── data/
    ├── assembly/         # Assembly outputs
    ├── mapping/          # Mapping outputs
    ├── bins/            # MAG bins
    ├── instrain/        # InStrain outputs
    └── bakta/           # Annotation output

```

## Installation

### Using Docker (Recommended)
```bash
# Pull the container
docker pull roles/ihmp-pipeline:latest

# Run the pipeline
docker run -v $(pwd)/data:/data roles/ihmp-pipeline
```

### Manual Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/roles/ihmp.git
   cd ihmp
   ```

2. Install dependencies:
   ```bash
   # Create conda environment
   conda env create -f snakemake/envs/instrain.yml
   ```

## Usage

### Running the Pipeline

1. **Configure the Pipeline**
   ```yaml
   # config.yaml
   patients:
     - "patient1"
     - "patient2"
   
   params:
     min_coverage: 10
     min_frequency: 0.05
     min_slope: 0.01
   ```

2. **Run the Pipeline**
   ```bash
   # Using Snakemake
   snakemake --use-conda --configfile snakemake/config/config.yaml
   
   # Using Docker
   docker run -v $(pwd)/data:/data roles/ihmp-pipeline
   ```

### Output Files

The pipeline generates the following key outputs:

1. **Assembly Outputs**
   - `data/assembly/{patient}/final_contigs.fa`: Co-assembled contigs
   - `data/assembly/{patient}/contigs_gt1kb.fa`: Filtered contigs
   - `data/assembly/{patient}/assembly.stb`: Strain tracking file

2. **Mapping Outputs**
   - `data/mapping/{patient}/{sample}.filtered.sorted.bam`: Processed BAM files
   - `data/mapping/{patient}/{sample}.filtered.flagstat.txt`: Mapping statistics

3. **InStrain Outputs**
   - `data/instrain/{patient}/each/{sample}/output/{sample}_SNVs.tsv`: SNV information
   - `data/instrain/{patient}/each/{sample}/output/{sample}_scaffold_info.tsv`: Scaffold information
   - `data/instrain/{patient}/combined/snv_info.tsv`: Combined SNV data
   - `data/instrain/{patient}/combined/scaffold_info.tsv`: Combined scaffold data

4. **SNP Tracking Outputs**
   - `data/instrain/{patient}/SNV_filtered.txt`: Filtered variants
   - `data/instrain/{patient}/SNV_filtered_trend.txt`: Temporal trends
   - `data/instrain/{patient}/SNV_filtered_trend_trajectory.txt`: Trajectory analysis
   - `data/instrain/{patient}/SNV_filtered_trend_trajectory_cluster_mapped.txt`: Clustered trajectories
   - `data/instrain/{patient}/SNV_filtered_trend_trajectory_cluster_mapped_mutation.txt`: Mutation types

## Development

### Adding New Rules
1. Create a new rule file in `snakemake/rules/`
2. Add the rule file to `snakemake/Snakefile`
3. Update wildcard functions in `snakemake/scripts/wildcards.py`

### Testing
```bash
# Run Snakemake dry-run
snakemake -n --configfile snakemake/config/config.yaml

# Run specific rule
snakemake -R rule_name --configfile snakemake/config/config.yaml
```

## Citation

If you use this pipeline in your research, please cite:

```
@software{ihmp2024,
  author = {Renee Oles},
  title = {iHMP Analysis Pipeline},
  year = {2024},
  url = {https://github.com/roles/ihmp}
}
```

## License

MIT 