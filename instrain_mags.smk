#!/usr/bin/env python
"""
Co‑assembly → annotation → mapping → binning → inStrain workflow
for single‑end (R1‑only) FASTQs.
"""

import os, re, pandas as pd
configfile: "patient_ids.yaml"

# ────────── CONSTANTS ──────────
BASE_DIR      = "/ddn_scratch/roles/strain_analysis/iHMP"
PATIENT_IDS   = config["patient_ids"]

FASTQ_LIST_DIR = os.path.join(BASE_DIR, "input", "samples")
ASSEMBLY_DIR   = os.path.join(BASE_DIR, "assembly")
REF_BAKTA_DIR  = os.path.join(BASE_DIR, "reference", "bakta")
DEPTH_DIR      = os.path.join(BASE_DIR, "depth")
BIN_DIR        = os.path.join(BASE_DIR, "bins")

COASSEMBLY       = lambda p: os.path.join(ASSEMBLY_DIR, p, "final_contigs.fa")
FILTERED_CONTIGS = lambda p: os.path.join(ASSEMBLY_DIR, p, "contigs_gt1kb.fa")
STB_FILE         = lambda p: os.path.join(ASSEMBLY_DIR, p, "assembly.stb")
DEPTH_FILE       = lambda p: os.path.join(DEPTH_DIR, f"{p}.depth.txt")
PATIENT_BIN_DIR  = lambda p: os.path.join(BIN_DIR, p)
INSTR_DIR        = lambda p: os.path.join(BASE_DIR, "instrain", p)

# ────────── CONDA ENVS ─────────
CONDA_MEGAHIT   = "/ddn_scratch/roles/strain_analysis/iHMP/config/megahit.yml"
CONDA_MAPPING   = "/ddn_scratch/roles/strain_analysis/iHMP/config/mapping.yml"
CONDA_INSTRAIN  = "/ddn_scratch/roles/strain_analysis/iHMP/config/instrain.yml"
CONDA_BAKTA     = "/ddn_scratch/roles/Panpiper/panpiper/workflow/envs/bakta.yml"
CONDA_METABAT   = "/ddn_scratch/roles/strain_analysis/iHMP/config/metabat2.yml"
CONDA_SEQKIT   = "/ddn_scratch/roles/strain_analysis/iHMP/config/seqkit.yml"

# ────────── HELPERS ────────────
def get_samples(patient: str):
    """
    Return dict {sample_id : R1_path} for every line in <patient>_full_path.txt>.
    """
    fp_file = os.path.join(FASTQ_LIST_DIR, f"{patient}_full_path.txt")
    df      = pd.read_csv(fp_file, header=None, names=["R1"])
    df["R1"] = df["R1"].str.strip()
    df["fname"] = df["R1"].apply(os.path.basename)
    df["sample"] = df["fname"].str.replace(
        r"(\.|_)R1\.fastq(\.gz)?$", "", regex=True
    )
    return dict(zip(df["sample"], df["R1"]))

# per‑sample convenience paths (wildcards left intact)
MAP_BAM   = lambda p, s: os.path.join(BASE_DIR, "mapping", p, f"{s}.all.bam")
SORT_BAM  = lambda p, s: os.path.join(BASE_DIR, "mapping", p, f"{s}.filtered.sorted.bam")
FLAGSTAT  = lambda p, s: os.path.join(BASE_DIR, "mapping", p, f"{s}.filtered.flagstat.txt")
SNV_FILE  = lambda p, s: os.path.join(INSTR_DIR(p), "each", s, "SNVs.tsv")
SCF_FILE  = lambda p, s: os.path.join(INSTR_DIR(p), "each", s, "scaffold_info.tsv")

def all_flagstats():
    return [
        FLAGSTAT(p, s) for p in PATIENT_IDS for s in get_samples(p)
    ]

# ────────── RULE ALL ───────────
rule all:
    input:
        # per‑patient artefacts
        expand(COASSEMBLY("{patient}"),            patient=PATIENT_IDS),
        expand(FILTERED_CONTIGS("{patient}"),      patient=PATIENT_IDS),
        expand(STB_FILE("{patient}"),              patient=PATIENT_IDS),
        expand(os.path.join(REF_BAKTA_DIR, "{patient}", "bakta.done"),
               patient=PATIENT_IDS),
        expand(DEPTH_FILE("{patient}"),            patient=PATIENT_IDS),
        expand(PATIENT_BIN_DIR("{patient}"),       patient=PATIENT_IDS),

        # per‑sample artefacts
        lambda wc: all_flagstats(),
        lambda wc: [SNV_FILE(p, s) for p in PATIENT_IDS for s in get_samples(p)],
        lambda wc: [SCF_FILE(p, s) for p in PATIENT_IDS for s in get_samples(p)],

        # merged tables
        expand(os.path.join(INSTR_DIR("{patient}"), "combined_scaffold_info.tsv"),
               patient=PATIENT_IDS),
        expand(os.path.join(INSTR_DIR("{patient}"), "combined_SNV_info.tsv"),
               patient=PATIENT_IDS)

rule *:
    log:
        lambda wildcards: f"logs/{wildcards.__rule__}.{getattr(wildcards, 'sample', 'na')}.{getattr(wildcards, 'patient', 'na')}.log"
    benchmark:
        lambda wildcards: f"benchmarks/{wildcards.__rule__}.{getattr(wildcards, 'sample', 'na')}.{getattr(wildcards, 'patient', 'na')}.txt"

# ────────── ASSEMBLY / FILTER / BAKTA ─────────────────────────────────────────
rule megahit_coassembly:
    input:
        lambda wc: list(get_samples(wc.patient).values())
    output: fa = COASSEMBLY("{patient}")
    params:
        outdir=os.path.join(ASSEMBLY_DIR, "{patient}"), kmin=27, kmax=77,
        kstep=10, merge="20,0.95", minlen=200
    conda:  CONDA_MEGAHIT
    threads: 8
    resources: mem="200G"
    shell:
        r"""
        mkdir -p {params.outdir}
        megahit \
          --force \
          -r $(printf "%s," {input}) \
          --min-contig-len {params.minlen} \
          --k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} \
          --merge-level {params.merge} \
          -o {params.outdir}
        mv {params.outdir}/final.contigs.fa {output.fa}
        """

rule filter_contigs:
    input:  fa = COASSEMBLY("{patient}")
    output: fa = FILTERED_CONTIGS("{patient}")
    conda:  CONDA_SEQKIT
    shell:  "seqkit seq -m 1000 {input.fa} -o {output.fa}"

rule make_stb:
    input:  fa = FILTERED_CONTIGS("{patient}")
    output: stb = STB_FILE("{patient}")
    run:
        with open(input.fa) as fin, open(output.stb, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(f"{line[1:].split()[0]}\t{wildcards.patient}\n")

rule bakta_coassembly:
    input:  fa = FILTERED_CONTIGS("{patient}")
    output: touch(os.path.join(REF_BAKTA_DIR, "{patient}", "bakta.done"))
    params:
        db="/ddn_scratch/roles/Panpiper/panpiper/databases/bakta2/db",
        outdir=os.path.join(REF_BAKTA_DIR, "{patient}"),
        prefix="{patient}"
    conda:  CONDA_BAKTA
    resources: mem="50G"
    shell:
        "mkdir -p {params.outdir} && bakta --db {params.db} --output {params.outdir} "
        "--prefix {params.prefix} --keep-contig-headers --force {input.fa}"

# ────────── BWA INDEX, MAP, FILTER ───────────────────────────────────────────
rule bwa_index:
    input:  fa = FILTERED_CONTIGS("{patient}")
    output: amb = FILTERED_CONTIGS("{patient}") + ".amb",
            ann = FILTERED_CONTIGS("{patient}") + ".ann",
            bwt = FILTERED_CONTIGS("{patient}") + ".bwt",
            pac = FILTERED_CONTIGS("{patient}") + ".pac",
            sa  = FILTERED_CONTIGS("{patient}") + ".sa"
    conda:  CONDA_MAPPING
    shell:  "bwa index {input.fa}"

rule map_reads_bwa:
    input:
        ref   = FILTERED_CONTIGS("{patient}"),
        amb   = FILTERED_CONTIGS("{patient}") + ".amb",
        ann   = FILTERED_CONTIGS("{patient}") + ".ann",
        bwt   = FILTERED_CONTIGS("{patient}") + ".bwt",
        pac   = FILTERED_CONTIGS("{patient}") + ".pac",
        sa    = FILTERED_CONTIGS("{patient}") + ".sa",
        r1    = lambda wc: get_samples(wc.patient)[wc.sample]
    output:
        bam_all = MAP_BAM("{patient}", "{sample}"),
        sorted   = SORT_BAM("{patient}", "{sample}"),
        flagstat = FLAGSTAT("{patient}", "{sample}"),
        index    = SORT_BAM("{patient}", "{sample}") + ".bai"
    threads: 8
    conda:  CONDA_MAPPING
    shell:
        r"""
        mkdir -p $(dirname {output.bam_all})
        bwa mem -t {threads} {input.ref} {input.r1} \
          | samtools view -b - > {output.bam_all}
        samtools flagstat {output.bam_all} > {output.flagstat}
        samtools view -b -F 4 -q 1 {output.bam_all} | samtools sort -o {output.sorted}
        samtools index {output.sorted}
        """

# ────────── DEPTH & BINNING ──────────────────────────────────────────────────
rule jgi_summarize_depths:
    input:
        bams=lambda wc: [SORT_BAM(wc.patient, s) for s in get_samples(wc.patient)]
    output:
        depth = DEPTH_FILE("{patient}")
    conda: CONDA_METABAT
    shell:
        """
        mkdir -p {DEPTH_DIR}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams}
        """

rule metabat2_binning:
    input:  contig = FILTERED_CONTIGS("{patient}"), depth = DEPTH_FILE("{patient}")
    output: directory(PATIENT_BIN_DIR("{patient}"))
    params: outdir = PATIENT_BIN_DIR("{patient}")
    conda:  CONDA_METABAT
    resources: mem="50G"
    shell:
        """
        mkdir -p {params.outdir}
        metabat2 -i {input.contig} -a {input.depth} -o {params.outdir}/bin \
                 -s 1500 -m 1500 --maxP 95 --minS 60 --maxEdges 200 \
                 --seed 1 --saveCls
        """

# ────────── INSTRAIN ────────────────────────────────────────────────────────
rule inStrain_profile:
    input:  bam_sorted = SORT_BAM("{patient}", "{sample}")
    output: snv = SNV_FILE("{patient}", "{sample}"), scaffold = SCF_FILE("{patient}", "{sample}")
    params: outdir = os.path.join(INSTR_DIR("{patient}"), "each", "{sample}", "output"),
           genome = FILTERED_CONTIGS("{patient}"),
           stb    = STB_FILE("{patient}")
    conda:  CONDA_INSTRAIN
    resources: mem="100G"
    threads: 8
    shell:
        r"""
        mkdir -p {params.outdir}
        inStrain profile {input.bam_sorted} {params.genome} \
            -s {params.stb} -o {params.outdir} \
            --min_cov 10 -p {threads} --database_mode -d \
            --pairing_filter all_reads \
            --skip_plot_generation
        """

# ────────── MERGE TABLES ────────────────────────────────────────────────────
rule combine_scaffold_info:
    input:
        lambda wc: [SCF_FILE(wc.patient, s) for s in get_samples(wc.patient)]
    output:
        combined = os.path.join(INSTR_DIR("{patient}"), "combined_scaffold_info.tsv")
    run:
        import pandas as pd, os
        dfs = [
            pd.read_csv(f, sep="\t").assign(Sample=os.path.basename(os.path.dirname(f)))
            for f in input if os.path.getsize(f)
        ]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined, sep="\t", index=False)

rule combine_SNV_info:
    input:
        lambda wc: [SNV_FILE(wc.patient, s) for s in get_samples(wc.patient)]
    output:
        combined = os.path.join(INSTR_DIR("{patient}"), "combined_SNV_info.tsv")
    run:
        import pandas as pd, os
        dfs = [
            pd.read_csv(f, sep="\t").assign(Sample=os.path.basename(os.path.dirname(f)))
            for f in input if os.path.getsize(f)
        ]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined, sep="\t", index=False)
