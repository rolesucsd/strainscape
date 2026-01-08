"""
MAGs rules for the iHMP pipeline.
Handles depth calculation and binning of metagenomic data.

This module contains rules for processing metagenomic assembled genomes (MAGs), including:
- Depth calculation using jgi_summarize_bam_contig_depths
- Binning using MetaBAT2
- Quality assessment of bins

Dependencies:
- Sorted BAM files from mapping
- Filtered contigs from assembly
- MetaBAT2 for binning
"""

import os

rule semibin2:
    """
    SemiBin2 in co-assembly mode:
    - Use `single_easy_bin` with MULTIPLE BAMs mapped to the SAME co-assembly.
      (Correct per SemiBin2 docs. Pretrained `--environment` is NOT used for co-assembly.)
    - Self-supervised training happens internally; we force reproducibility with --random-seed.
    - We set `-m/--min-len` to MIN_CONTIG_LEN to keep contig cutoff consistent with assembly filter.
    - `--sequencing-type short_read` to select short-read model.

    Parameters:
    - Minimum contig length: 1500bp
    
    Input:
        contig: Filtered contigs file
        depth: Depth file from jgi_summarize_depths
    Output:
        directory: Directory containing binned contigs
    Resources:
        mem: 50GB

    """
    input:
        contig = FILTERED_CONTIGS("{patient}"),
        bams   = lambda wc: [SORT_BAM(wc.patient, s) for s in sorted(read_samples(wc.patient))]
    output:
        directory(PATIENT_BIN_OUT("{patient}"))
    params:
        outdir = PATIENT_BIN_DIR("{patient}"),
        minlen=1500,
        engine="auto",          # auto-detect GPU, fall back to CPU
        seed=1
    threads: 16
    conda:
        config['conda_envs']['semibin']
    shell:
        r"""
        mkdir -p "{params.outdir}"
        SemiBin2 single_easy_bin \
          --sequencing-type short_read \
          --self-supervised \
          --random-seed {params.seed} \
          -m {params.minlen} \
          -t {threads} \
          -i {input.contig} \
          -b {input.bams} \
          -o {params.outdir} \
          --engine {params.engine}
        """

rule checkm2:
    """
    Classify and assess MAGs using CheckM2.

    This rule runs CheckM2 to estimate the taxonomy and quality of bins.

    Input:
        bin_dir: Directory with binned MAGs (FASTA files)
    Output:
        results: Directory with CheckM2 output including taxonomy and QC
    """
    input:
        bin_dir = PATIENT_BIN_OUT("{patient}")
    output:
        results = directory(CHECKM2_OUT("{patient}")),
        tsv = os.path.join(CHECKM2_OUT("{patient}"), "quality_report.tsv")
    params:
        db=config['reference']['checkm2_db'],
        threads=8
    conda:
        config['conda_envs']['checkm2']
    resources:
        mem = "32G"
    shell:
        """
        checkm2 predict --threads {params.threads} -x fa.gz --database_path {params.db} -i {input.bin_dir} -o {output.results}
        """

rule make_bin_txt:
    """
    Generate a bin.txt file (scaffold→bin mapping) and merge
    CheckM2 metrics (Completeness, Contamination, Genome_Size).
    """
    input:
        bin_dir    = PATIENT_BIN_OUT("{patient}"),
        checkm2_tsv = os.path.join(CHECKM2_OUT("{patient}"), "quality_report.tsv")
    output:
        bin_file = BIN_TXT("{patient}")
    resources:
        mem = "5G"
    log:
        os.path.join(config['paths']['log_dir'], "make_bin_txt_{patient}.log")
    shell:
        """
        python ../strainscape/make_bin_txt.py \
          --bin_dir {input.bin_dir} \
          --checkm2-tsv {input.checkm2_tsv} \
          --output_file {output.bin_file} \
          --log_file {log} 2>&1
        """

rule gtdbtk_classify:
    """
    Assign taxonomy to MAGs with GTDB-Tk directly from the bin directory.
    """
    input:
        bin_dir = PATIENT_BIN_OUT("{patient}"),
        db_dir  = config['reference']['gtdbtk_db']
    output:
        out_dir = directory(GTDBTK_OUT("{patient}"))
    params:
        scratch = os.path.join(GTDBTK_OUT("{patient}"), "tmp/")
    threads: 1            # scheduler can downscale; GTDB-Tk uses this for HMM/align
    resources:
        mem = "150G"       # raise only if logs show OOM
    conda:
        config['conda_envs']['gtdbtk']
    shell:
        r"""
        mkdir -p "{params.scratch}/tmp"
        export GTDBTK_DATA_PATH="{input.db_dir}"
        gtdbtk classify_wf \
          --genome_dir "{input.bin_dir}" \
          --out_dir "{output.out_dir}" \
          --cpus {threads} \
          --pplacer_cpus 1 \
          --scratch_dir "{params.scratch}" \
          -x fa.gz
        """