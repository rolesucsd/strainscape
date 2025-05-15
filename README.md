<!-- README.md – StrainScape -->

<p align="center">
  <img src="docs/img/strainscape_logo.svg" alt="StrainScape logo" height="120">
</p>

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

## ✨ What is StrainScape?

StrainScape is a reproducible pipeline that:

1. **Assembles metagenomes** per participant (MEGAHIT)  
2. **Maps reads** and filters alignments (BWA + SAMtools)  
3. **Profiles strains** (inStrain)  
4. **Bins genomes** (MetaBAT2)  
5. **Annotates genes** (Bakta)  
6. **Performs downstream population genetics** in Python/R  
7. Generates publication-ready figures and summary tables.

It was built to explore evolutionary dynamics in the **Integrative Human
Microbiome Project (iHMP)** but is easily adapted to any time-series
metagenome dataset.

---

## 📂 Repository layout

```

.
├── snakemake/           # workflow rules, envs/, config/
├── docs/                # user guide, figures, logo
├── scripts/             # helper bash / python utilities
├── src/
│   ├── python/          # analysis & plotting modules
│   └── r/               # statistical models & ggplot themes
├── config.yaml          # top-level user settings
└── run\_pipeline.sh      # one-liner launcher

````

---

## ⚡ Quick start

> Tested on Linux / macOS, Conda ≥ 23.1, Snakemake ≥ 7.30

```bash
# 1. clone
git clone https://github.com/rolesucsd/strainscape.git
cd strainscape

# 2. set paths & parameters
nano config.yaml         # edit input FASTQ locations, sample sheet…

# 3. create the base env with Snakemake + Mamba
conda env create -f snakemake/envs/core.yml
conda activate strainscape

# 4. run the whole workflow on 32 cores
./run_pipeline.sh --cores 32
````

Outputs land in `results/`:

```
results/
├── assembly/
├── mapping/
├── instrain/
├── bins/
├── figures/         # PNG / PDF plots
└── tables/          # TSV / CSV summary files
```

---

## 🛠️ Advanced usage

| Task                        | Command                                           |
| --------------------------- | ------------------------------------------------- |
| Dry-run the DAG             | `snakemake -n`                                    |
| Resume failed jobs          | `snakemake --keep-going`                          |
| Clean intermediates         | `snakemake --delete-temp-output`                  |
| Override a rule’s resources | `snakemake map_reads_bwa --resources mem_mb=8000` |
| Render an HTML report       | `snakemake --report report.html`                  |

---

## 🔬 Workflow diagram

```
FASTQ → assembly → mapping ┐
                           ├─► inStrain profile ─► population tables
bins  ◄────────────────────┘
```

*(see `docs/workflow.svg` for the full DAG)*

---

## 📑 Configuration keys (excerpt)

| Key                | Description                       | Example                           |
| ------------------ | --------------------------------- | --------------------------------- |
| `samples:`         | TSV mapping `sample_id  fastq_R1` | `CSM67UB9  /path/sample.fastq.gz` |
| `threads:`         | Default threads per rule          | `8`                               |
| `instrain.min_cov` | Minimum read depth                | `10`                              |
| `filters.length`   | Keep contigs ≥ bp                 | `1000`                            |

Full reference in [`docs/config_reference.md`](docs/config_reference.md).

---

## 🚀 Development

```bash
# create dev env with pre-commit, pytest, coverage…
conda env create -f environment-dev.yml
conda activate strainscape-dev
pre-commit install
```

* **Unit tests:** `pytest -q`
* **Linting & black:** automatically run by pre-commit
* **Docs:** `mkdocs serve`

Contributions via pull request are welcome – please open an issue first
if you plan a large change.

---

## 📖 Citing StrainScape

If StrainScape helps your research, please cite:

```bibtex
@software{olles_scape_2025,
  author       = {Renee A. Oles and contributors},
  title        = {StrainScape: a reproducible pipeline for strain-level
                  microbiome analysis},
  year         = {2025},
  publisher    = {GitHub},
  journal      = {Zenodo},
  doi          = {10.5281/zenodo.12345678},
  url          = {https://github.com/rolesucsd/strainscape}
}
```

---

## 📝 License

This project is released under the MIT License – see [`LICENSE`](LICENSE) for details.

---

## 🙏 Acknowledgements

* README template inspired by the [https://snakemake.bio](https://snakemake.bio) community.
* Uses open-source tools: **Snakemake**, **inStrain**, **MEGAHIT**, **MetaBAT2**, **Bakta**, **pandas**, **ggplot2**.

Feel free to open an issue if you hit a bug or have a feature request.