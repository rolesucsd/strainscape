# iHMP Analysis Pipeline

A comprehensive pipeline for analyzing strain-level evolution in the iHMP (Integrative Human Microbiome Project) dataset.

## Overview

This pipeline processes metagenomic data through several stages:
1. **Data Processing** (Snakemake): Assembly, mapping, and initial variant calling
2. **Python Analysis**: SNP analysis, mutation typing, and trajectory analysis
3. **R Analysis**: Statistical analysis and visualization

## Directory Structure

```
.
├── config.yaml              # Main configuration file
├── run_pipeline.sh         # Top-level script to run everything
├── src/
│   ├── python/            # Python analysis code
│   │   ├── analysis/     # Core analysis functions
│   │   ├── visualization/# Plotting functions
│   │   └── utils/        # Utility functions
│   └── r/                # R analysis code
│       ├── analysis/     # Statistical analysis functions
│       ├── visualization/# Plotting functions
│       └── utils/        # Utility functions
├── scripts/
│   └── bash/            # Utility bash scripts
├── snakemake/
│   ├── rules/           # Snakemake rule files
│   ├── config/          # Snakemake configuration
│   └── envs/            # Conda environment files
├── data/
│   ├── raw/            # Raw input data
│   └── processed/      # Processed intermediate files
├── results/
│   ├── intermediate/   # Intermediate results
│   └── final/         # Final results and figures
└── docs/
    ├── api/           # API documentation
    └── examples/      # Usage examples
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/roles/ihmp.git
   cd ihmp
   ```

2. Install dependencies:
   ```bash
   # Create conda environment
   conda env create -f snakemake/envs/instrain.yml
   
   # Install R package
   Rscript -e "renv::restore()"
   ```

## Usage

### Quick Start

Run the complete pipeline:
```bash
./run_pipeline.sh
```

### Step-by-Step

1. Run Snakemake workflow:
   ```bash
   snakemake --use-conda --configfile config.yaml
   ```

2. Run Python analysis:
   ```bash
   python src/python/ihmp/iHMP_main.py --config config.yaml
   ```

3. Run R analysis:
   ```bash
   Rscript src/r/run_analysis.R --config config.yaml
   ```

## Configuration

Edit `config.yaml` to configure:
- Input/output paths
- Analysis parameters
- Resource limits

## Output

Results are organized in the `results/` directory:
- `intermediate/`: Intermediate files from each step
- `final/`: Final results and visualizations

## Development

### Adding New Analysis

1. Add Python functions to `src/python/analysis/`
2. Add R functions to `src/r/analysis/`
3. Update `run_pipeline.sh` and `run_analysis.R`

### Testing

```bash
# Run Python tests
pytest tests/python/

# Run R tests
Rscript -e "testthat::test_dir('tests/r/')"
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