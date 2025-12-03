# StrainScape pipeline schematic

The diagram below summarizes the end-to-end workflow used by StrainScape, from raw FASTQ inputs through mutation-type reporting. Each stage is annotated with the core tools and outputs so you can trace how files move through the Snakemake rules.

![StrainScape pipeline schematic](pipeline_schematic.svg)

## Stage-by-stage notes

- **Input discovery** – `config/config.yaml` provides `patient_ids` and paths, while per-patient `*_full_path.txt` manifests list the FASTQs consumed by downstream rules. Sample names are parsed in `snakemake/Snakefile` via `read_samples`. 
- **Co-assembly** – `megahit_coassembly` pools each patient's reads and runs MEGAHIT with a 27–77 k-mer sweep, emitting `final_contigs.fa` under `data_dir/assembly/{patient}`.【F:snakemake/rules/assembly.smk†L16-L44】
- **Contig filtering and annotation** – `filter_contigs` trims assemblies to >1 kb, `make_stb` creates scaffold-to-bin mappings for InStrain, and `bakta_coassembly` annotates genes to produce patient-specific TSV/FNA outputs.【F:snakemake/rules/assembly.smk†L46-L103】
- **Mapping** – `bwa_index` prepares indices on filtered contigs; `map_reads_bwa` aligns reads, writes raw and filtered/sorted BAMs, and produces per-sample flagstat QC files for each patient/sample combination.【F:snakemake/rules/mapping.smk†L15-L70】
- **Depth and binning** – `jgi_summarize_depths` aggregates coverage across sorted BAMs, then `metabat2_binning` and `checkm2` generate quality-scored MAG bins alongside `bins.txt` scaffold-bin maps.【F:snakemake/rules/mags.smk†L15-L89】
- **InStrain profiling** – `inStrain_profile` runs per-sample strain profiling on sorted BAMs using filtered contigs + STB, and the combine rules merge scaffold and SNV tables per patient while tolerating empty inputs.【F:snakemake/rules/instrain.smk†L16-L66】
- **SNP tracking** – `process_scaffolds` filters combined scaffold info with bin/metadata joins; `merge_snv_info` limits SNVs to passing scaffolds; `calculate_trends`, `map_genes`, and `analyze_mutation_types` derive time trends, gene mappings, and mutation classifications; final rules merge patient-level outputs into cohort-wide tables.【F:snakemake/rules/snp_tracking.smk†L6-L122】【F:snakemake/rules/snp_tracking.smk†L124-L158】

## How to use the schematic

1. Confirm paths and patients in `config/config.yaml` before running.
2. Follow the top-to-bottom flow when debugging: check assembly outputs, then mapping/BAM quality, then binning depth files, InStrain tables, and finally SNP tracking derivatives.
3. The arrows show data passing between modules; reuse of patient-level artifacts (bins, annotations, scaffold filters) feeds the downstream per-sample or per-patient rules highlighted in the diagram.
