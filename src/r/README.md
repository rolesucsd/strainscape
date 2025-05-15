# iHMP Analysis Tools

An R package for analyzing iHMP (Integrative Human Microbiome Project) data, focusing on strain-level evolution in the human microbiome.

## Installation

```r
# Install from GitHub
devtools::install_github("roles/ihmp")
```

## Features

- Mutation accrual analysis
- Gene-level metrics calculation
- Species filtering and visualization
- Consistent visualization theme and color palette

## Usage

### Mutation Accrual Analysis

```r
# Analyze mutation accrual patterns
results <- analyze_mutation_accrual(
  combined_df = mutation_data,
  metadata = sample_metadata,
  outdir = "output"
)
```

### Gene Metrics

```r
# Calculate gene-level metrics
gene_metrics <- calculate_gene_metrics(
  snp_metrics = snp_data,
  outdir = "output"
)
```

### Species Filtering

```r
# Filter and analyze species with multiple strains
filtered_data <- filter_species_multiple_strains(
  nucleotide_div_filter = diversity_data,
  scaffold_data_edit = scaffold_data,
  outdir = "output"
)
```

## Visualization

The package provides a consistent theme and color palette for all visualizations:

```r
# Use the iHMP theme
ggplot(data, aes(x, y)) +
  geom_point() +
  theme_ihmp()

# Use the iHMP color palette
colors <- ihmp_colors()
```

## Development

### Running Tests

```r
devtools::test()
```

### Building Documentation

```r
devtools::document()
```

## License

MIT

## Citation

If you use this package in your research, please cite:

```
@software{ihmp2024,
  author = {Renee Oles},
  title = {iHMP Analysis Tools},
  year = {2024},
  url = {https://github.com/roles/ihmp}
}
``` 