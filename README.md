<!-- README.md – StrainScape -->

<h1 align="center">StrainScape</h1>
<p align="center">
  <em>End-to-end Snakemake workflow for strain-level analysis of the iHMP microbiome cohort</em>
</p>

<p align="center">
  <a href="https://github.com/rolesucsd/strainscape/actions"><img alt="CI" src="https://github.com/rolesucsd/strainscape/actions/workflows/ci.yml/badge.svg"></a>
  <a href="LICENSE"><img alt="MIT license" src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
  <a href="https://doi.org/10.5281/zenodo.12345678"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.12345678.svg" alt="DOI"></a>
</p>

---

Below is a proposed overhaul of your `README.md`. It’s organized into clear sections, includes a concise summary up front, and reflects the core components of the Strainscape pipeline.

> **Note:** Be sure to update any placeholder URLs, author names, or contact info as appropriate for your project.

---

## 🚀 Summary

Strainscape is an end‐to‐end Snakemake-powered workflow for strain-level analysis of metagenomic datasets, originally built around the iHMP cohort. It automates everything from raw FASTQ processing through assembly, mapping, and InStrain profiling to downstream mutation tracking, trend analysis, and visualization. The pipeline leverages Conda environments for reproducible software stacks, Python for data wrangling and mutation analysis, and R for statistical summaries and plots.

---

## 📦 Features

* **Modular Snakemake rules** for assembly, mapping, binning, InStrain profiling, and mutation tracking
* **Conda-managed environments** for each major toolset (e.g., InStrain, Bakta)
* Automatic **merging** of per-sample scaffolds and SNV tables
* **Trend and trajectory** computation of mutation frequencies over time
* **Gene mapping** and **mutation-type classification** (synonymous, non-synonymous, etc.)
* **Configurable** via a single YAML file with patient IDs, paths, and resource parameters
* **Logging** and **flagstat** summaries for mapping QC

---

## 📝 Requirements

* Linux / macOS
* [Snakemake ≥ 8.4.8](https://github.com/snakemake/snakemake)
* [Conda](https://docs.conda.io) (Mamba recommended)
* Python 3.8+ (with pandas, numpy, etc.)
* R 4.0+ (for optional downstream plotting)

---

## 🛠️ Installation

1. **Clone the repo**

   ```bash
   git clone https://github.com/rolesucsd/strainscape.git
   cd strainscape
   ```

   ([GitHub][1])

2. **Create a top-level Conda env** (optional)

   ```bash
   mamba create -n strainscape python=3.10 snakemake=8.4.8 -c conda-forge -c bioconda
   conda activate strainscape
   ```

   ([GitHub][1])

3. **Install per-rule environments**
   Snakemake will automatically build/use the environments defined under `snakemake/envs/` during the run (via `--use-conda`).

---

## ⚙️ Configuration

All settings live in **`config/config.yaml`**:

* `patient_ids`: list of patient identifiers
* `paths`: base directories for raw data, logs, and outputs
* `metadata`: path to sample metadata file
* Conda environment mappings for each rule

*Edit `config/config.yaml` to point at your data directories and metadata.*

---

## ▶️ Usage

### 1. Dry-run

```bash
snakemake --use-conda --cores 1 --configfile config/config.yaml --dry-run
```

### 2. Full execution

```bash
snakemake --use-conda --cores 16 --configfile config/config.yaml
```

*Sample log messages* will appear in `logs/` and a summary report in `pipeline.log`. ([GitHub][1])

---

## 📂 Directory Structure

```
.
├── config/
│   └── config.yaml          # Pipeline settings
├── docs/                    # Documentation and examples
├── snakemake/
│   ├── rules/               # Individual rule files (.smk)
│   ├── wildcards.smk        # Path–wildcard definitions
│   └── utils.py             # Shared Python helpers
├── snakemake/envs/          # Conda env YAMLs for each rule
├── run_pipeline.sh          # Convenience wrapper script
├── src/                     # Python scripts for merging, analysis
├── results/                 # Final outputs (InStrain tables, trends)
└── logs/                    # Snakemake logs and flagstat outputs
```

---

## 🔄 Pipeline Overview

1. **Assembly** (`assembly.smk`): de novo contig generation with MEGAHIT
2. **Mapping** (`mapping.smk`): BWA indexing, read alignment, BAM sorting, and flagstat QC
3. **MAGS** (`mags.smk`): binning with MetaBAT2 and JGI summarization
4. **InStrain** (`instrain.smk`): per‐sample SNV/scaffold profiling
5. **SNP Tracking** (`snp_tracking.smk`): combine tables, filter, compute trends, map to genes, classify mutation types

---

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add/modify Snakemake rules or Python scripts
4. Update `config/config.yaml` examples as needed
5. Open a pull request

---

## 📜 License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details. ([GitHub][1])

---

## ✉️ Contact

For questions or issues, please open an issue on GitHub or contact **[roles@ucsd.edu](mailto:roles@ucsd.edu)**.

[1]: https://github.com/rolesucsd/strainscape "GitHub - rolesucsd/strainscape"

---

## 🙏 Acknowledgements

* README template inspired by the [https://snakemake.bio](https://snakemake.bio) community.
* Uses open-source tools: **Snakemake**, **inStrain**, **MEGAHIT**, **MetaBAT2**, **Bakta**, **pandas**, **ggplot2**.

Feel free to open an issue if you hit a bug or have a feature request.
