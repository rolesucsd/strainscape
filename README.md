# iHMP Analysis Tools

A collection of tools for analyzing iHMP (Integrative Human Microbiome Project) data.

## Installation

```bash
pip install -e .
```

## Usage

### Creating Interval Trees

```bash
create-trees --genes_file path/to/genes.txt --output_dir path/to/output
```

### Summarizing Coverage

```bash
summarize-coverage --coverage_dir path/to/coverage_results
```

## Development

This package is organized as follows:

- `ihmp_analysis/`: Main package directory
  - `trees.py`: Functions for creating and managing interval trees
  - `coverage.py`: Functions for analyzing coverage data
  - `utils.py`: Common utility functions

## License

[Add your license here] 