# StrainScape

A comprehensive toolkit for analyzing strain evolution in microbiome data, with a focus on iHMP (Integrative Human Microbiome Project) data.

## Repository Structure

```
strainscape/
├── snakemake/                 # Snakemake workflow for initial data processing
│   ├── instrain_mags.smk     # Main workflow file
│   └── rules/                # Individual workflow rules
│
├── src/                      # Source code for analysis
│   ├── python/               # Python analysis modules
│   │   ├── trees/           # Interval tree creation and management
│   │   ├── coverage/        # Coverage analysis
│   │   ├── filtering/       # BAM and pangenome filtering
│   │   └── mapping/         # Mapping and alignment utilities
│   │
│   └── r/                   # R visualization and analysis
│       ├── figures/         # Figure generation scripts
│       └── analysis/        # Statistical analysis scripts
│
├── scripts/                  # Utility scripts
│   ├── bash/                # Shell scripts for setup and utilities
│   └── python/              # Python utility scripts
│
├── tests/                   # Test suite
│   ├── python/             # Python tests
│   └── r/                  # R tests
│
└── docs/                    # Documentation
    ├── usage/              # Usage guides
    └── api/                # API documentation
```

## Installation

```bash
# Clone the repository
git clone https://github.com/rolesucsd/strainscape.git
cd strainscape

# Install Python package
pip install -e .

# Install R dependencies
R -e "install.packages(c('tidyverse', 'ggplot2', 'devtools'))"
```

## Usage

### Running the Snakemake Workflow

```bash
snakemake -s snakemake/instrain_mags.smk --configfile config.yaml
```

### Python Analysis

```python
from strainscape import trees, coverage

# Create interval trees
trees.compute_and_save_trees(genes_file, output_dir)

# Analyze coverage
coverage.summarize_coverage(coverage_dir)
```

### R Visualization

```r
# Load the package
library(strainscape)

# Generate figures
plot_coverage_summary(coverage_data)
plot_strain_evolution(evolution_data)
```

## Development

This project uses a development workflow with two main branches:
- `main`: Production-ready code
- `develop`: Development branch for new features

### Contributing

1. Create a feature branch from `develop`
2. Make your changes
3. Add tests
4. Submit a pull request

## License

[Add your license here] 